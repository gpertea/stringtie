/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <vector>
#include <memory>
#include <limits>
#include <cmath>
#include <cstdlib>

#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/jellyfish.hpp>

#include <jflib/multiplexed_io.hpp>
#include <gzip_stream.hpp>

#include <src/mer_database.hpp>
#include <src/error_correct_reads.hpp>
#include <src/error_correct_reads_cmdline.hpp>
#include <src/verbose_log.hpp>

namespace err = jellyfish::err;

using jellyfish::mer_dna;
typedef std::vector<const char*> file_vector;
typedef jellyfish::stream_manager<file_vector::const_iterator> stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> read_parser;


typedef uint64_t hkey_t;
typedef uint64_t hval_t;

static args_t args;

// Poisson term. Computes a term of the Poisson distribution: e^{-\lambda} \lambda^i / i!.
double poisson_term(double lambda, unsigned int i) {
  static const double facts[11] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 };
  static const double tau       = 6.283185307179583;

  if(i < 11)
    return exp(-lambda) * pow(lambda, i) / facts[i];
  else
    return exp(-lambda + i) * pow(lambda / i, i) / sqrt(tau * i);
}

// Contaminant database. If a Jellyfish database is given, return true
// iff the k-mer is in the database. With no database, it always
// returns false.
class contaminant_check {
public:
  virtual ~contaminant_check() { }

  virtual bool is_contaminant(const mer_dna& m) = 0;
  virtual bool is_contaminant(const mer_dna& m, mer_dna& tmp) = 0;
  virtual void debug(const char*msg) = 0;
};

class contaminant_no_database : public contaminant_check {
public:
  virtual ~contaminant_no_database() { }
  virtual bool is_contaminant(const mer_dna& m) { return false; }
  virtual bool is_contaminant(const mer_dna& m, mer_dna& tmp) { return false; }
  virtual void debug(const char* msg) { std::cerr << msg << " no database" << std::endl; }
};

class contaminant_database : public contaminant_check {
  mer_array ary_;
public:
  contaminant_database(binary_reader& reader, size_t size) : ary_(size, mer_dna::k() * 2, 0, 126) {
    while(reader.next()) {
      if(!ary_.set(reader.key()))
        throw std::runtime_error("Size of hash for contaminant too small");
    }
  }
  virtual ~contaminant_database() { }
  virtual void debug(const char* msg) { std::cerr << msg << std::endl; }
  virtual bool is_contaminant(const mer_dna& m) { return ary_.has_key(m); }
  virtual bool is_contaminant(const mer_dna& m, mer_dna& tmp) {
    size_t id;
    return ary_.get_key_id(m, &id, tmp);
  }
};

template<class instance_t>
class error_correct_t : public jellyfish::thread_exec {
  read_parser            _parser;
  int                    _skip;
  int                    _good;
  int                    _anchor;
  std::string            _prefix;
  int                    _min_count;
  int			 _cutoff;
  char                   _qual_cutoff;
  int                    _window;
  int                    _error;
  bool                   _gzip;
  const database_query*  _mer_database;
  contaminant_check*     _contaminant;
  bool                   _trim_contaminant;
  int                    _homo_trim;
  double                 _collision_prob; // collision probability = a priori error rate / 3
  double                 _poisson_threshold;
  bool                   _no_discard;

  jflib::o_multiplexer * _output;
  jflib::o_multiplexer * _log;
public:
  error_correct_t(int nb_threads, stream_manager& streams) :
    _parser(4 * nb_threads, 100, 1, streams),
    _skip(0), _good(1), _min_count(1), _cutoff(4), _window(0), _error(0), _gzip(false),
    _mer_database(0), _contaminant(0), _trim_contaminant(false),
    _homo_trim(std::numeric_limits<int>::min()), _no_discard(false) { }

private:
  // Open the data (error corrected reads) and log files. Default to
  // STDOUT and STDERR if none specified.
  std::ostream* open_file(const std::string prefix, const char* suffix,
                          const std::string def) {
    std::ostream* res;
    std::string file;

    if(prefix.empty())
      file = def;
    else {
      file = prefix;
      file += suffix;
    }
    if(_gzip) {
      if(!prefix.empty())
        file += ".gz";
      res = new gzipstream(file.c_str());
    } else
      res = new std::ofstream(file.c_str());
    if(!res->good())
      throw std::runtime_error(err::msg() << "Failed to open file '" << file << "'" << err::no);
    res->exceptions(std::ios::eofbit|std::ios::failbit|std::ios::badbit);
    return res;
  }

public:
  void do_it(int nb_threads) {
    // Make sure they are deleted when done
    std::unique_ptr<std::ostream> details(open_file(_prefix, ".log", "/dev/fd/2"));
    std::unique_ptr<std::ostream> output(open_file(_prefix, ".fa", "/dev/fd/1"));
    // Multiplexers, same thing
    std::unique_ptr<jflib::o_multiplexer>
      log_m(new jflib::o_multiplexer(details.get(), 3 * nb_threads, 1024));
    std::unique_ptr<jflib::o_multiplexer>
      output_m(new jflib::o_multiplexer(output.get(), 3 * nb_threads, 1024));
    _log    = log_m.get();
    _output = output_m.get();

    exec_join(nb_threads);
  }

  virtual void start(int id) {
    instance_t(*this, id).start();
  }

  error_correct_t& skip(int s) { _skip = s; return *this; }
  error_correct_t& good(int g) { _good = g; return *this; }
  error_correct_t& anchor(int a) { _anchor = a; return *this; }
  error_correct_t& prefix(const char *s) { _prefix = s; return *this; }
  error_correct_t& prefix(const std::string s) { _prefix = s; return *this; }
  error_correct_t& min_count(int m) { _min_count = m; return *this; }
  error_correct_t& cutoff(int p) { _cutoff = p; return *this; }
  error_correct_t& qual_cutoff(char c) { _qual_cutoff = c; return *this; }
  //  error_corret_t & advance(int a) { _advance = a; return *this; }
  error_correct_t& window(int w) { _window = w; return *this; }
  error_correct_t& error(int e) { _error = e; return *this; }
  error_correct_t& gzip(bool g) { _gzip = g; return *this; }
  error_correct_t& mer_database(database_query* q) { _mer_database = q; return *this; }
  error_correct_t& contaminant(contaminant_check* c) { _contaminant = c; return *this; }
  error_correct_t& trim_contaminant(bool t) { _trim_contaminant = t; return *this; }
  error_correct_t& homo_trim(int t) { _homo_trim = t; return *this; }
  error_correct_t& collision_prob(double cp) { _collision_prob = cp; return *this; }
  error_correct_t& poisson_threshold(double t) { _poisson_threshold = t; return *this; }
  error_correct_t& no_discard(bool d) { _no_discard = d; return *this; }

  read_parser& parser() { return _parser; }
  int skip() const { return _skip; }
  int good() const { return _good; }
  int anchor() const { return _anchor; }
  const std::string & prefix() const { return _prefix; }
  int min_count() const { return _min_count; }
  int cutoff() const { return _cutoff; }
  char qual_cutoff() const { return _qual_cutoff; }
  //  int advance() const { return _advance; }
  int window() const { return _window ? _window : mer_dna::k(); }
  int error() const { return _error ? _error : mer_dna::k() / 2; }
  bool gzip() const { return _gzip; }
  const database_query* mer_database() const { return _mer_database; }
  contaminant_check* contaminant() const { return _contaminant; }
  bool trim_contaminant() const { return _trim_contaminant; }
  bool do_homo_trim() const { return _homo_trim != std::numeric_limits<int>::min(); }
  int homo_trim() const { return _homo_trim; }
  double collision_prob() const { return _collision_prob; }
  double poisson_threshold() const { return _poisson_threshold; }
  bool no_discard() const { return _no_discard; }

  jflib::o_multiplexer& output() { return *_output; }
  jflib::o_multiplexer& log() { return *_log; }
};

class error_correct_instance {
public:
  typedef error_correct_t<error_correct_instance> ec_t ;

private:
  ec_t&   _ec;
  //  int     _id;
  size_t  _buff_size;
  char*   _buffer;
  kmer_t  _tmp_mer;
  mer_dna _tmp_mer_dna;

  static const char* error_contaminant;
  static const char* error_no_starting_mer;
  static const char* error_homopolymer;

public:
  error_correct_instance(ec_t& ec, int id) :
    _ec(ec), _buff_size(0), _buffer(0) { }
    //    _ec(ec), _id(id), _buff_size(0), _buffer(0) { }
  ~error_correct_instance() {
    free(_buffer);
  }

  void start() {
    jflib::omstream output(_ec.output());
    jflib::omstream details(_ec.log());
    kmer_t          mer, tmer;

    uint64_t nb_reads = 0;
    while(true) {
      read_parser::job job(_ec.parser());
      if(job.is_empty()) break;
      for(size_t i = 0; i < job->nb_filled; ++i) {
        if(i % 2 == 0)
          output << jflib::endr;
        const std::string& header   = job->data[i].header;
        const std::string& sequence = job->data[i].seq;
        const char* const  seq_s    = sequence.c_str();
        const char* const  seq_e    = seq_s + sequence.size();
        const char* const  qual_s   = job->data[i].qual.c_str();

        nb_reads++;
        insure_length_buffer(sequence.size());

        const char* error = "";
        const char *input = seq_s + _ec.skip();
        char       *out   = _buffer + _ec.skip();
        //Prime system. Find and write starting k-mer
        if(!find_starting_mer(mer, input, seq_e, out, &error)) {
          details << "Skipped " << header << ": " << error << "\n";
          details << jflib::endr;
          if(_ec.no_discard())
            output << ">" << header << "\nN\n";
          continue;
        }
        // Extend forward and backward
        tmer = mer;
        const ssize_t start_off = input - seq_s;
        forward_log fwd_log(_ec.window(), _ec.error());
        char *end_out =
          extend(forward_mer(tmer),
                 forward_ptr<const char>(input),
                 forward_ptr<const char>(qual_s + start_off),
                 forward_counter(start_off),
                 forward_ptr<const char>(seq_e),
                 forward_ptr<char>(out), fwd_log,
                 &error);
        if(!end_out) {
          details << "Skipped " << header << ": " << error << "\n";
          details << jflib::endr;
          if(_ec.no_discard())
            output << ">" << header << "\nN\n";
          continue;
        }
        assert(input > seq_s + mer_dna::k());
        assert(out > _buffer + mer_dna::k());
        assert(input - seq_s == out - _buffer);
        tmer = mer;
        backward_log bwd_log(_ec.window(), _ec.error());
        char *start_out =
          extend(backward_mer(tmer),
                 backward_ptr<const char>(input - mer_dna::k() - 1),
                 backward_ptr<const char>(qual_s + start_off - mer_dna::k() - 1),
                 backward_counter(start_off - mer_dna::k() - 1),
                 backward_ptr<const char>(seq_s - 1),
                 backward_ptr<char>(out - mer_dna::k() - 1), bwd_log,
                 &error);
        if(!start_out) {
          details << "Skipped " << header << ": " << error << "\n";
          details << jflib::endr;
          if(_ec.no_discard())
            output << ">" << header << "\nN\n";
          continue;
        }
        start_out++;
        assert(start_out >= _buffer);
        assert(_buffer + _buff_size >= end_out);

        if(_ec.do_homo_trim()) {
          end_out = homo_trim(_buffer, start_out, end_out, fwd_log, bwd_log, &error);
          if(!end_out) {
            details << "Skipped " << header << ": " << error << "\n";
            details << jflib::endr;
            if(_ec.no_discard())
              output << ">" << header << "\nN\n";
            continue;
          }
        }
        assert(end_out >= _buffer);
        assert(_buffer + _buff_size >= end_out);

        output << ">" << header
               << " " << fwd_log << " " << bwd_log << "\n"
               << substr(start_out, end_out) << "\n";
      } // for(size_t i...  Loop over reads in job
    } // while(true)... loop over all jobs
    details.close();
    output.close();
  }

private:
  enum log_code { OK, TRUNCATE, ERROR };

  template<typename dir_mer, typename elog, typename counter>
  log_code check_contaminant(dir_mer& mer, elog& log, const counter& cpos, const char**error) {
    if(_ec.contaminant()->is_contaminant(mer.canonical(), _tmp_mer_dna)) {
      if(_ec.trim_contaminant()) {
        log.truncation(cpos);
        return TRUNCATE;
      }
      *error = error_contaminant;
      return ERROR;
    }
    return OK;
  }

  // Log a substitution.
  template<typename dir_mer, typename out_dir_ptr, typename elog, typename counter>
  log_code log_substitution(dir_mer& mer, out_dir_ptr& out, elog& log, const counter& cpos, int from, int to,
                            const char**error) {
    if(from == to)
      return OK;

    mer.replace(0, to);
    switch(log_code c = check_contaminant(mer, log, cpos, error)) {
    case OK: break;
    default: return c;
    }

    if(log.substitution(cpos, from >= 0 ? mer_dna::rev_code(from) : 'N', to >= 0 ? mer_dna::rev_code(to) : 'N')) {
      int diff = log.remove_last_window();
      out = out - diff;
      log.truncation(cpos - diff);
      return TRUNCATE;
    }
    return OK;
  }

  // Extend and correct read. Copy from input to out. mer should be
  // represent a "good" starting k-mer at the input position.
  // out point to the next character to be written.
  template<typename dir_mer, typename in_dir_ptr, typename out_dir_ptr,
           typename counter, typename elog>
  char * extend(dir_mer mer, in_dir_ptr input, in_dir_ptr qual,
                counter pos, in_dir_ptr end,
                out_dir_ptr out, elog &log, const char** error) {
    counter  cpos       = pos;
    uint32_t prev_count = _ec.mer_database()->get_val(mer.canonical());

    for( ; input < end; ++input, ++qual) {
      const char base = *input;

      if(base == '\n')
        continue;
      cpos = pos;
      ++pos;

      const int ori_code = mer_dna::code(base);
      mer.shift(ori_code >= 0 ? ori_code : 0);
      if(ori_code >= 0) {
        switch(check_contaminant(mer, log, cpos, error)) {
        case OK: break;
        case TRUNCATE: goto done;
        case ERROR: return 0;
        }
      }
      uint64_t counts[4];
      int      ucode = 0;
      int      level;

      const int count = _ec.mer_database()->get_best_alternatives(mer, counts, ucode, level);

      // No coninuation whatsoever, trim.
      if(count == 0) {
        log.truncation(cpos);
        goto done;
      }

      if(count == 1) { // One continuation. Is it an error?
        prev_count = counts[ucode];
        switch(log_substitution(mer, out, log, cpos, ori_code, ucode, error)) {
        case OK: break;
        case TRUNCATE: goto done;
        case ERROR: return 0;
        }
        *out++ = mer.base(0);
        continue;
      }
      // We get here if there is more than one alternative base to try
      // at some level. Note that if the current base is low quality
      // and all alternatives are higher quality, then the current
      // base will have a count of zero. If the current base is non N
      // and has a high count or previous count is low enough that
      // continuity does not apply, output current base. But if the current
      // base has count of zero, all alternatives are low quality and prev_count is low
      // then trim
      if(ori_code >= 0){ //if the current base is valid base (non N)
	if(counts[ori_code] > (uint64_t)_ec.min_count()) {
          if(counts[ori_code] >= (uint32_t)_ec.cutoff() || *qual >= _ec.qual_cutoff()) {
            *out++ = mer.base(0);
            continue;
          }
          // Now we ask for a probability of getting
          // counts[ori_code] errors with p=1/300 in sum_counts
          // trials.  If this probability is < 10e-6, do not correct
          double p = (double)(counts[0] + counts[1] + counts[2] + counts[3]) * _ec.collision_prob();
          const double prob = poisson_term(p, counts[ori_code]);
          if(prob < _ec.poisson_threshold()) {
            *out++ = mer.base(0);
            continue;
          }
	} else if(level == 0 && counts[ori_code] == 0) {
          // definitely an error and all alternatives are low quality
          log.truncation(cpos);
          goto done;
	}
      } else if(level == 0) { //if all alternatives are low quality
	log.truncation(cpos);
	goto done;
      }

      // We get here if there are multiple possible substitutions, the
      // original count is low enough and the previous count is high (good) or
      // the current base is an N
      // We find out all alternative bases
      // that have a continuation at the same or better level.  Then
      // if the current base is N, pick the one with the highest
      // count that is the most similar to the prev_count,
      // otherwise pick the one with the most similar count.
      // If no alternative continues, leave the base alone.
      int          check_code               = ori_code;
      bool         success                  = false;
      uint64_t     cont_counts[4]; //here we record the counts for the continuations
      bool         continue_with_correct_base[4];
      int          read_nbase_code          = -1;
      bool         candidate_continuations[4];
      unsigned int ncandidate_continuations = 0;

      //here we determine what the next base in the read is
      if(input + 1 < end)
        read_nbase_code = mer_dna::code(*(input + 1));

      for(int i = 0; i < 4; ++i) {
        cont_counts[i]                = 0;
        continue_with_correct_base[i] = false;
        if(counts[i] <= (uint64_t)_ec.min_count())
          continue;
        check_code = i;
        // Check that it continues at least one more base with that quality
        _tmp_mer     = mer.kmer();
        dir_mer nmer = _tmp_mer;
        nmer.replace(0, check_code);
        // Does not matter what we shift, check all alternative anyway.
        nmer.shift(0);

        uint64_t   ncounts[4];
        int        nucode = 0;
        int        nlevel;
        const int ncount = _ec.mer_database()->get_best_alternatives(nmer, ncounts, nucode, nlevel);
        if(ncount > 0 && nlevel >= level) {
          continue_with_correct_base[i] = read_nbase_code >= 0 && ncounts[read_nbase_code] > 0;
          success                       = true;
          cont_counts[i]                = counts[i];
        }
      }

      if(success) {
        // We found at least one alternative base that continues now
        // we look for all alernatives that have a count closest to
        // the previous count first we determine the count that is the
        // closest to the current count but in the special case of
        // prev_count == 1 we simply pick the largest count
        check_code           = -1;
        uint32_t _prev_count = prev_count<=(uint64_t)_ec.min_count() ? std::numeric_limits<uint32_t>::max() : prev_count;
        int      min_diff    = std::numeric_limits<int>::max();
        for(int  i = 0; i < 4; ++i) {
          candidate_continuations[i] = false;
          if(cont_counts[i] > 0)
            min_diff = std::min(min_diff, (int)std::abs((long)cont_counts[i] - (long)_prev_count));
        }

        //we now know the count that is the closest, now we determine how many alternatives have this count
        for(uint32_t  i = 0; i < 4; i++) {
          if(std::abs((long)cont_counts[i] - (long)_prev_count) == min_diff){
            candidate_continuations[i] = true;
            ++ncandidate_continuations;
            check_code=i;
          }
        }

        //do we have more than one good candidate? if we do then check which one continues with the correct base
        if(ncandidate_continuations>1 && read_nbase_code >= 0)
          for(int  i = 0; i < 4; ++i){
            if(candidate_continuations[i]){
              if(!continue_with_correct_base[i])
                --ncandidate_continuations;
              else
                check_code = i;
            }
          }

        //fail if we still have more than one candidate
        if(ncandidate_continuations != 1)
          check_code = -1;

        if(check_code >= 0) {
          switch(log_substitution(mer, out, log, cpos, ori_code, check_code, error)) {
          case OK: break;
          case TRUNCATE: goto done;
          case ERROR: return 0;
          }
        }
      }
      if(ori_code < 0 && check_code < 0) {// if invalid base and no good sub found
        log.truncation(cpos);
        goto done;
      }
      *out++ = mer.base(0);
    }

  done:
    return out.ptr();
  }

  char* homo_trim(const char* start, char* out_start, char* out_end,
		  forward_log& fwd_log, backward_log& bwd_log, const char** error) {
    int   max_homo_score = std::numeric_limits<int>::min();
    char* max_pos        = 0;
    int   homo_score     = 0;
    char* ptr            = out_end - 1;
    char  pbase          = mer_dna::code(*ptr);

    for(--ptr; ptr >= out_start; --ptr) {
      char cbase = mer_dna::code(*ptr);
      homo_score += ((pbase == cbase) << 1) - 1; // Add 1 if same as last, -1 if not
      pbase       = cbase;
      if(homo_score > max_homo_score) {
	max_homo_score = homo_score;
	max_pos        = ptr;
      }
    }

    if(max_homo_score < _ec.homo_trim())
      return out_end; // Not a high score -> return without truncation
    // assert(max_pos >= out_start);
    // assert(max_pos >= start);
    if(max_pos < out_start) {
      *error = error_homopolymer;
      return 0;
    }
    fwd_log.force_truncate(forward_counter(max_pos - start));
    bwd_log.force_truncate(backward_counter(max_pos - start));
    fwd_log.truncation(forward_counter(max_pos - start));
    return max_pos;
  }

  void insure_length_buffer(size_t len) {
    if(len <= _buff_size)
      return;

    _buff_size = len > 2 * _buff_size ? len + 100 : 2 * _buff_size;
    _buffer    = (char *)realloc(_buffer, _buff_size);
    if(!_buffer)
      throw std::runtime_error(err::msg() << "Buffer allocation failed, size " << _buffer << err::no);
  }

  bool find_starting_mer(kmer_t &mer, const char * &input, const char *end, char * &out,
			 const char** error) {
    while(input < end) {
      for(int i = 0; input < end && i < (int)mer_dna::k(); ++i) {
	char base = *input++;
	*out++ = base;
	if(!mer.shift_left(base))
	  i = -1;        // If an N, skip to next k-mer
      }
      int found = 0;
      while(input < end) {
	bool contaminated = _ec.contaminant()->is_contaminant(mer.canonical(), _tmp_mer_dna);
	if(contaminated && !_ec.trim_contaminant()) {
	  *error = error_contaminant;
	  return false;
	}

	if(!contaminated) {
	  hval_t val = _ec.mer_database()->get_val(mer.canonical());

	  found = (int)val >= _ec.anchor() ? found + 1 : 0;
	  if(found >= _ec.good())
	    return true;
	}

	char base = *input++;
	*out++ = base;
	if(!mer.shift_left(base))
	  break;
      }
    }

    *error = error_no_starting_mer;
    return false;
  }
};

const char* error_correct_instance::error_contaminant     = "Contaminated read";
const char* error_correct_instance::error_no_starting_mer = "No high quality mer";
const char* error_correct_instance::error_homopolymer     = "Entire read is an homopolymer";

unsigned int compute_poisson_cutoff__(const val_array_raw& counts, double collision_prob, double poisson_threshold) {
  auto     it_end   = counts.end();
  uint64_t distinct = 0;
  uint64_t total    = 0;
  for(auto it = counts.begin(); it != it_end; ++it) {
    if((*it & 0x1) && (*it >= 2)) {
      distinct += 1;
      total    += *it >> 1;
    }
  }
  const double coverage = (double)total / (double)distinct;
  vlog << "distinct mers:" << distinct << " total mers:" << total << " estimated coverage:" << coverage;
  const double lambda = coverage * collision_prob;
  vlog << "lambda:" << lambda << " collision_prob:" << collision_prob << " poisson_threshold:" << poisson_threshold;
  for(unsigned int x = 2; x < 1000; ++x)
    if(poisson_term(lambda, x) < poisson_threshold)
      return x + 1;
  return 0;
}

unsigned int compute_poisson_cutoff(const val_array_raw& counts, double collision_prob, double poisson_threshold) {
  vlog << "Computing Poisson cutoff";
  unsigned int res = compute_poisson_cutoff__(counts, collision_prob, poisson_threshold);
  return res;
}

int main(int argc, char *argv[])
{
  args.parse(argc, argv);

  if(args.qual_cutoff_char_given && args.qual_cutoff_char_arg.size() != 1)
    args_t::error("The qual-cutoff-char must be one ASCII character.");
  if(args.qual_cutoff_value_given && args.qual_cutoff_value_arg > (uint32_t)std::numeric_limits<char>::max())
    args_t::error("The qual-cutoff-value must be in the range 0-127.");
  const char qual_cutoff = args.qual_cutoff_char_given ? args.qual_cutoff_char_arg[0] :
    (args.qual_cutoff_value_given ? (char)args.qual_cutoff_value_arg : std::numeric_limits<char>::max());

  verbose_log::verbose = args.verbose_flag;
  vlog << "Loading mer database";
  database_query mer_database(args.db_arg, args.no_mmap_flag);
  mer_dna::k(mer_database.header().key_len() / 2);

  // Open contaminant database.
  std::unique_ptr<contaminant_check> contaminant;
  contaminant.reset(new contaminant_no_database());
  if(args.contaminant_given) {
    vlog << "Loading contaminant sequences";
    std::ifstream contaminant_file(args.contaminant_arg);
    if(!contaminant_file.good())
      err::die(err::msg() << "Failed to open contaminant file '" << args.contaminant_arg << "'");
    jellyfish::file_header header(contaminant_file);
    if(header.format() != binary_dumper::format)
      err::die(err::msg() << "Contaminant format expected '" << binary_dumper::format << "'");
    if(mer_dna::k() * 2 != header.key_len())
      err::die(err::msg() << "Contaminant mer length (" << (header.key_len() / 2)
          << ") different than correction mer length (" << mer_dna::k() << ")");
    binary_reader reader(contaminant_file, &header);
    contaminant.reset(new contaminant_database(reader, header.size()));
  }

  stream_manager streams(args.sequence_arg.cbegin(), args.sequence_arg.cend(), 1);

  const unsigned int cutoff =   args.cutoff_given ?
    args.cutoff_arg :
    compute_poisson_cutoff(mer_database.vals(), args.apriori_error_rate_arg / 3,
                           args.poisson_threshold_arg / args.apriori_error_rate_arg);
  vlog << "Using cutoff of " << cutoff;
  if(cutoff == 0 && !args.cutoff_given)
    err::die("Cutoff computation failed. Pass it explicitly with -p switch.");

  error_correct_instance::ec_t correct(args.thread_arg, streams);
  correct.skip(args.skip_arg).good(args.good_arg)
    .anchor(args.anchor_count_arg)
    .prefix(args.output_given ? (std::string)args.output_arg : "")
    .min_count(args.min_count_arg)
    .cutoff(cutoff)
    .qual_cutoff(qual_cutoff)
    .window(args.window_arg)
    .error(args.error_arg)
    .gzip(args.gzip_flag)
    .mer_database(&mer_database)
    .contaminant(contaminant.get())
    .trim_contaminant(args.trim_contaminant_flag)
    .homo_trim(args.homo_trim_given ? args.homo_trim_arg : std::numeric_limits<int>::min())
    .collision_prob(args.apriori_error_rate_arg / 3)
    .poisson_threshold(args.poisson_threshold_arg)
    .no_discard(args.no_discard_flag);
  vlog << "Correcting reads";
  correct.do_it(args.thread_arg);
  vlog << "Done";

  return 0;
 }

