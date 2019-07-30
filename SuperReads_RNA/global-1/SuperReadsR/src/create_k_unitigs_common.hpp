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

#ifndef _CREATE_K_UNITIGS_COMMON_H_
#define _CREATE_K_UNITIGS_COMMON_H_

#include <jellyfish/thread_exec.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>

#include <jflib/multiplexed_io.hpp>
#include <jellyfish/atomic_field.hpp>

using jellyfish::mer_dna;
using jellyfish::thread_exec;
using jellyfish::mer_dna_bloom_counter;
using jellyfish::mer_dna_bloom_filter;

// Wrapper around mer_iterator to satisfy common interface
typedef std::vector<const char*> file_vector;
typedef jellyfish::stream_manager<file_vector::const_iterator> stream_manager;
typedef jellyfish::mer_overlap_sequence_parser<stream_manager> sequence_parser;
typedef jellyfish::mer_iterator<sequence_parser, mer_dna> mer_iterator;
class read_mers {
  mer_iterator stream_;

public:
  read_mers(sequence_parser& parser, int id) : stream_(parser, true) { }
  operator bool() const { return (void*)stream_ != 0; }
  const mer_dna* operator->() const { return stream_.operator->(); }
  const mer_dna& operator*() const { return stream_.operator*(); }
  read_mers& operator++() {
    ++stream_;
    return *this;
  }
};

// Wrapper around bloom filter class to have set compatible insert
// method.
class mer_bloom {
  mer_dna_bloom_filter bf_;

public:
  mer_bloom(double fp, size_t size) : bf_(fp, size) { }
  std::pair<unsigned int, bool> insert(const mer_dna& m) {
    unsigned int r = bf_.insert(m);
    return std::make_pair(r, r == 0);
  }
};


// Insert a mer in a set and return true if the k-mer is new in the
// set.
template<typename set_type>
bool insert_canonical(set_type& set, const mer_dna& mer) {
  return set.insert(mer.get_canonical()).second;
}

/* Read k-mers and store them in a map. The map_type must have the
   operator[]. The content returned must have the prefix ++. All this
   needs to be multi-thread safe.
 */
template<typename map_type, typename parser_type, typename stream_type>
class populate_mer_set : public thread_exec {
  int          mer_len_;
  parser_type& parser_;
  map_type&    set_;

public:
  populate_mer_set(int mer_len, map_type& set, parser_type& parser) :
    mer_len_(mer_len), parser_(parser), set_(set)
  { }

  void start(int thid) {
    for(stream_type stream(parser_, thid); stream; ++stream)
      ++set_[*stream];
    set_.done();
  }
};

/* - mer_counts_type maps k-mer to counts. Has operator[] returning
     the count.

   - used_type is a set type with operator insert. (set compatible)
   - end_points_type is a set type with operator insert. (set compatible)
   - parser_type is the type to generate iterator of k-mers
   - stream_type is the k-mer iteratorp type
   - args_type is the type of switches on command line.
 */
template<typename mer_counts_type, typename used_type, typename end_points_type,
         typename parser_type, typename stream_type,
         typename args_type>
class create_k_unitig : public jellyfish::thread_exec {
  const mer_counts_type&        counts_; // Counts for k-mers
  used_type&                    used_mers_; // Mark all k-mers whether they have been visited already
  end_points_type&              end_points_; // End points of k-unitigs, to ouput only once
  int                           threads_;
  parser_type&                  parser_;
  jflib::o_multiplexer          output_multiplexer_;
  jflib::atomic_field<uint64_t> unitig_id_;
  const args_type&              args_;

  enum direction { forward = 1, backward = -1 };
  static direction rev_direction(direction dir) { return (direction)-dir; }

public:
  create_k_unitig(const mer_counts_type& counts, used_type& used, end_points_type& ends,
                  int threads, parser_type& parser, std::ostream& output,
                  const args_type& args) :
    counts_(counts),
    used_mers_(used),
    end_points_(ends),
    threads_(threads),
    parser_(parser),
    output_multiplexer_(&output, 3 * threads, 4096),
    unitig_id_(0),
    args_(args)
  { }

  virtual void start(int thid) {
    //    mer_stream<mer_dna, read_parser> stream(mer_dna::k(), parser_);
    stream_type     stream(parser_, thid);
    jflib::omstream output(output_multiplexer_);
    mer_dna         current;
    mer_dna         continuation;
    mer_dna         tmp;

    for( ; stream; ++stream) {
      auto is_new = used_mers_.insert(*stream);
      if(!is_new.second)
        continue;
      // Never start a unitig on low count
      if(counts_[*stream] < args_.quality_threshold_arg)
        continue;
      current = *stream;
      //printf("Current mer %s\n",(current.to_str()).c_str());
      // Grow unitig if a starting (branching) mer
      if(starting_mer(forward, current)) {
        //printf("Starting forward\n");
        grow_unitig(backward, current, output);
      } else if(starting_mer(backward, current)) {
        //printf("Starting reverse\n");
        grow_unitig(forward, current, output);
      }
      // Unique continuation on both sides -> middle of k-unitig: do nothing
    }
  }

private:
  // Check all k-mers extending in one direction. If unique
  // continuation, store it in cont and return true. A unique
  // continuation may have a count less than the min quality
  // threshold, as will be reported in the count output argument.
  //
  // On the other hand, continuation with count less than the min
  // quality threshold do not create a branch compared to a high
  // quality continuation. I.e., if the threshold is 2, and the counts
  // are as follow:
  //
  // 0, 1, 0, 0 -> unique continuation, count of 1
  // 0, 1, 0, 1 -> no unique continuation
  // 2, 0, 0, 0 -> unique continuation, count of 2
  // 2, 1, 0, 0 -> unique continuation, count of 2
  // 2, 0, 3, 0 -> no unique continuation
  //
  // Otherwise return false and the value of cont is undetermined. If
  // true is returned and count is not NULL, the count of the unique
  // continuation mer is stored in the pointer.
  bool next_mer(const direction dir, const mer_dna& start, mer_dna& cont,
                unsigned int* count = 0) {
    int     index;
    mer_dna cont_comp(start);
    cont_comp.reverse_complement();
    cont = start;

    if(dir == forward) {
      cont.shift_left(0);
      cont_comp.shift_right(0);
      index = 0;
    } else {
      cont.shift_right(0);
      cont_comp.shift_left(0);
      index = cont.k() - 1;
    }
    auto base = cont.base(index); // Point to first or last base. Correct base to change
    auto base_comp = cont_comp.base(cont.k() - 1 - index);

    int          nb_cont    = 0, nb_low_cont = 0;
    int          code       = 0, low_code = 0;
    unsigned int cont_count = 0, low_cont_count = 0;
    for(int i = 0; i < 4; ++i) {
      base      = i;
      base_comp = mer_dna::complement(i);

      unsigned int cont_count_ = counts_[cont < cont_comp ? cont : cont_comp];
      if(cont_count_ >= args_.quality_threshold_arg) {
        if(++nb_cont > 1)
          return false;
        code       = i;
        cont_count = cont_count_;
      } else if(cont_count_ > 0) {
        ++nb_low_cont;
        low_code       = i;
        low_cont_count = cont_count_;
      }
    }

    if(nb_cont == 1) {
      base = code;
      if(count)
        *count = cont_count;
      return true;
    } else if(nb_cont == 0 && nb_low_cont == 1) {
        base = low_code;
        if(count)
          *count = low_cont_count;
        return true;
    }
    return false;
  }

  // Return true if m is a starting mer in the given dir: a mer is a
  // starting mer if it is branching forward, or backward, or dries
  // out, maybe after some number of low count mer to skip.
  bool starting_mer(direction dir, mer_dna m) {
    int     low_cont = args_.cont_on_low_arg;
    mer_dna tmp1, tmp2;

    while(true) {
      unsigned int count = 0;
      if(!next_mer(dir, m, tmp1, &count))
        return true;
      if(count >= args_.quality_threshold_arg) {
        if(!next_mer(rev_direction(dir), tmp1, tmp2, &count))
          return true;
        break;
      }
      if(--low_cont < 0)
        return true;
      m = tmp1;
    }
    return false;
  }

  void grow_unitig(const direction dir, const mer_dna& start, jflib::omstream& output) {
    bool start_new = insert_canonical(end_points_, start);
    if(!start_new)
      return;
    //printf("Starting unitig from %s\n",(start.to_str()).c_str());
    mer_dna            mer1(start);
    mer_dna            mer2;
    mer_dna            mer3;
    mer_dna           *current = &mer1;
    mer_dna           *cont    = &mer2;
    unsigned int       count   = 0;
    unsigned int       low_run = 0;
    unsigned int       index   = dir == forward ? 0 : start.k() - 1;
    std::string        seq;
    std::set<mer_dna>  set; // Current set of used mers to avoid endless loop

    insert_canonical(set, *current);

    while(true) {
      insert_canonical(used_mers_, *current);
      //if(!insert_canonical(set, *current))
        //break; // loop. Don't output anything  //AZ modified to avoid loss of k-unitigs from loops
      if(!next_mer(dir, *current, *cont, &count))
        break;
      //printf("Count next fwd %d mer %s\n",count,(current->to_str()).c_str());
      if(!next_mer(rev_direction(dir), *cont, mer3))
        break;
      //printf("Count next rev %d mer %s\n",count,(cont->to_str()).c_str());
      //printf("Mer3 %s\n",(mer3.to_str()).c_str());
      // This can happen (only) with continuation on low. It does not
      // create a branch as far as next_mer is concerned if one low
      // count and one high count, but it still a branch in this case:
      // there are two way to go through that region
      if(mer3 != *current)
        break;

      if(!insert_canonical(set, *cont)) //AZ we hit a loop -- do not extend
          break;
          
      seq += (char)cont->base(index);
      //printf("Seq cont is %s low_run is %d\n",seq.c_str(),low_run);
      if(count < args_.quality_threshold_arg) {
        if(++low_run > args_.cont_on_low_arg)
          break;
      } else
        low_run = 0;

      std::swap(current, cont);
    }
    //printf("Done %s\n",seq.c_str());
    // Erase trailing low quality bases if any and reset current to be
    // the actual last k-mer (with only good quality bases). Needed
    // for the test of already written k-unitigs to be accurate.
    if(low_run > 0) {
      seq.erase(seq.size() - std::min((unsigned int)seq.size(), low_run));
      if(seq.size() >= current->k()) {
        if(dir == forward) {
          *current = seq.substr(seq.size() - current->k());
        } else {
          std::string end = seq.substr(seq.size() - current->k());
          std::string rev(end.rbegin(), end.rend());
          *current = rev;
        }
      } else {
        // Sequence does not contain a full k-mer. Need to recreate it
        // from the starting k-mer and seq by shifting.
        *current = start;
        for(auto it = seq.begin(); it != seq.end(); ++it)
          if(dir == forward) {
            current->shift_left(*it);
          } else {
            current->shift_right(*it);
          }
      }
    }

    // If the last k-mer has been used in a k-unitigs already, this
    // means two threads are working on the same unitigs, starting
    // from opposite ends. Output only if the current thread has the
    // "largest" end k-mer.
    bool end_new = insert_canonical(end_points_, *current);
    if(!end_new) {
      if(start.get_canonical() < current->get_canonical())
        return;
    }
    // Output results
    if(start.k() + seq.length() < args_.min_len_arg)
      return;
    uint64_t id = (unitig_id_ += 1) - 1;
    output << ">" << id << " length:" << (start.k() + seq.size()) << "\n";
    if(dir == backward) {
      std::string reversed(seq.rbegin(), seq.rend());
      output << reversed << start << "\n";
    } else {
      output << start << seq << "\n";
    }
    output << jflib::endr;
  }
};

#endif /* _CREATE_K_UNITIGS_COMMON_H_ */
