/***** This code was generated by Yaggo. Do not edit ******/

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



#ifndef __CMDLINE_PARSE_HPP__
#define __CMDLINE_PARSE_HPP__

#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdexcept>
#include <string>
#include <limits>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>

class cmdline_parse {
 // Boiler plate stuff. Conversion from string to other formats
  static bool adjust_double_si_suffix(double &res, const char *suffix) {
    if(*suffix == '\0')
      return true;
    if(*(suffix + 1) != '\0')
      return false;

    switch(*suffix) {
    case 'a': res *= 1e-18; break;
    case 'f': res *= 1e-15; break;
    case 'p': res *= 1e-12; break;
    case 'n': res *= 1e-9;  break;
    case 'u': res *= 1e-6;  break;
    case 'm': res *= 1e-3;  break;
    case 'k': res *= 1e3;   break;
    case 'M': res *= 1e6;   break;
    case 'G': res *= 1e9;   break;
    case 'T': res *= 1e12;  break;
    case 'P': res *= 1e15;  break;
    case 'E': res *= 1e18;  break;
    default: return false;
    }
    return true;
  }

  static double conv_double(const char *str, ::std::string &err, bool si_suffix) {
    char *endptr = 0;
    errno = 0;
    double res = strtod(str, &endptr);
    if(errno) {
      err.assign(strerror(errno));
      return (double)0.0;
    }
    bool invalid =
      si_suffix ? !adjust_double_si_suffix(res, endptr) : *endptr != '\0';
    if(invalid) {
      err.assign("Invalid character");
      return (double)0.0;
    }
    return res;
  }

  static int conv_enum(const char* str, ::std::string& err, const char* const strs[]) {
    int res = 0;
    for(const char* const* cstr = strs; *cstr; ++cstr, ++res)
      if(!strcmp(*cstr, str))
        return res;
    err += "Invalid constant '";
    err += str;
    err += "'. Expected one of { ";
    for(const char* const* cstr = strs; *cstr; ++cstr) {
      if(cstr != strs)
        err += ", ";
      err += *cstr;
    }
    err += " }";
    return -1;
  }

  template<typename T>
  static bool adjust_int_si_suffix(T &res, const char *suffix) {
    if(*suffix == '\0')
      return true;
    if(*(suffix + 1) != '\0')
      return false;

    switch(*suffix) {
    case 'k': res *= (T)1000; break;
    case 'M': res *= (T)1000000; break;
    case 'G': res *= (T)1000000000; break;
    case 'T': res *= (T)1000000000000; break;
    case 'P': res *= (T)1000000000000000; break;
    case 'E': res *= (T)1000000000000000000; break;
    default: return false;
    }
    return true;
  }

  template<typename T>
  static T conv_int(const char *str, ::std::string &err, bool si_suffix) {
    char *endptr = 0;
    errno = 0;
    long long int res = strtoll(str, &endptr, 0);
    if(errno) {
      err.assign(strerror(errno));
      return (T)0;
    }
    bool invalid =
      si_suffix ? !adjust_int_si_suffix(res, endptr) : *endptr != '\0';
    if(invalid) {
      err.assign("Invalid character");
      return (T)0;
    }
    if(res > ::std::numeric_limits<T>::max() ||
       res < ::std::numeric_limits<T>::min()) {
      err.assign("Value out of range");
      return (T)0;
    }
    return (T)res;
  }

  template<typename T>
  static T conv_uint(const char *str, ::std::string &err, bool si_suffix) {
    char *endptr = 0;
    errno = 0;
    while(isspace(*str)) { ++str; }
    if(*str == '-') {
      err.assign("Negative value");
      return (T)0;
    }
    unsigned long long int res = strtoull(str, &endptr, 0);
    if(errno) {
      err.assign(strerror(errno));
      return (T)0;
    }
    bool invalid =
      si_suffix ? !adjust_int_si_suffix(res, endptr) : *endptr != '\0';
    if(invalid) {
      err.assign("Invalid character");
      return (T)0;
    }
    if(res > ::std::numeric_limits<T>::max()) {
      err.assign("Value out of range");
      return (T)0;
    }
    return (T)res;
  }

  template<typename T>
  static ::std::string vec_str(const std::vector<T> &vec) {
    ::std::ostringstream os;
    for(typename ::std::vector<T>::const_iterator it = vec.begin();
        it != vec.end(); ++it) {
      if(it != vec.begin())
        os << ",";
      os << *it;
    }
    return os.str();
  }

  class string : public ::std::string {
  public:
    string() : ::std::string() {}
    explicit string(const ::std::string &s) : std::string(s) {}
    explicit string(const char *s) : ::std::string(s) {}
    int as_enum(const char* const strs[]) {
      ::std::string err;
      int res = conv_enum((const char*)this->c_str(), err, strs);
      if(!err.empty())
        throw ::std::runtime_error(err);
      return res;
    }


    uint32_t as_uint32_suffix() const { return as_uint32(true); }
    uint32_t as_uint32(bool si_suffix = false) const {
      ::std::string err;
      uint32_t res = conv_uint<uint32_t>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to uint32_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    uint64_t as_uint64_suffix() const { return as_uint64(true); }
    uint64_t as_uint64(bool si_suffix = false) const {
      ::std::string err;
      uint64_t res = conv_uint<uint64_t>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to uint64_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int32_t as_int32_suffix() const { return as_int32(true); }
    int32_t as_int32(bool si_suffix = false) const {
      ::std::string err;
      int32_t res = conv_int<int32_t>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int32_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int64_t as_int64_suffix() const { return as_int64(true); }
    int64_t as_int64(bool si_suffix = false) const {
      ::std::string err;
      int64_t res = conv_int<int64_t>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int64_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int as_int_suffix() const { return as_int(true); }
    int as_int(bool si_suffix = false) const {
      ::std::string err;
      int res = conv_int<int>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    long as_long_suffix() const { return as_long(true); }
    long as_long(bool si_suffix = false) const {
      ::std::string err;
      long res = conv_int<long>((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to long_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    double as_double_suffix() const { return as_double(true); }
    double as_double(bool si_suffix = false) const {
      ::std::string err;
      double res = conv_double((const char*)this->c_str(), err, si_suffix);
      if(!err.empty()) {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to double_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
  };

public:
  int                            mer_arg;
  bool                           mer_given;
  uint64_t                       nb_mers_arg;
  bool                           nb_mers_given;
  int                            threads_arg;
  bool                           threads_given;
  const char *                   output_arg;
  bool                           output_given;
  uint32_t                       quality_threshold_arg;
  bool                           quality_threshold_given;
  uint32_t                       cont_on_low_arg;
  bool                           cont_on_low_given;
  uint64_t                       min_len_arg;
  bool                           min_len_given;
  bool                           gzip_flag;
  double                         false_positive_arg;
  bool                           false_positive_given;
  const char *                   load_arg;
  bool                           load_given;
  ::std::vector<const char *>    input_arg;
  typedef ::std::vector<const char *>::iterator input_arg_it;
  typedef ::std::vector<const char *>::const_iterator input_arg_const_it;

  enum {
    START_OPT = 1000,
    LOAD_OPT
  };

  cmdline_parse() :
    mer_arg(0), mer_given(false),
    nb_mers_arg(0), nb_mers_given(false),
    threads_arg((int)1), threads_given(false),
    output_arg(""), output_given(false),
    quality_threshold_arg((uint32_t)2), quality_threshold_given(false),
    cont_on_low_arg((uint32_t)0), cont_on_low_given(false),
    min_len_arg(0), min_len_given(false),
    gzip_flag(false),
    false_positive_arg((double)0.01), false_positive_given(false),
    load_arg(""), load_given(false),
    input_arg()
  { }

  cmdline_parse(int argc, char* argv[]) :
    mer_arg(0), mer_given(false),
    nb_mers_arg(0), nb_mers_given(false),
    threads_arg((int)1), threads_given(false),
    output_arg(""), output_given(false),
    quality_threshold_arg((uint32_t)2), quality_threshold_given(false),
    cont_on_low_arg((uint32_t)0), cont_on_low_given(false),
    min_len_arg(0), min_len_given(false),
    gzip_flag(false),
    false_positive_arg((double)0.01), false_positive_given(false),
    load_arg(""), load_given(false),
    input_arg()
  { parse(argc, argv); }

  void parse(int argc, char* argv[]) {
    static struct option long_options[] = {
      {"mer", 1, 0, 'm'},
      {"nb-mers", 1, 0, 'n'},
      {"threads", 1, 0, 't'},
      {"output", 1, 0, 'o'},
      {"quality-threshold", 1, 0, 'q'},
      {"cont-on-low", 1, 0, 'c'},
      {"min-len", 1, 0, 'l'},
      {"gzip", 0, 0, 'g'},
      {"false-positive", 1, 0, 'f'},
      {"load", 1, 0, LOAD_OPT},
      {"help", 0, 0, 'h'},
      {"usage", 0, 0, 'U'},
      {"version", 0, 0, 'V'},
      {0, 0, 0, 0}
    };
    static const char *short_options = "hVUm:n:t:o:q:c:l:gf:";

    ::std::string err;
#define CHECK_ERR(type,val,which) if(!err.empty()) { ::std::cerr << "Invalid " #type " '" << val << "' for [" which "]: " << err << "\n"; exit(1); }
    while(true) {
      int index = -1;
      int c = getopt_long(argc, argv, short_options, long_options, &index);
      if(c == -1) break;
      switch(c) {
      case ':':
        ::std::cerr << "Missing required argument for "
                  << (index == -1 ? ::std::string(1, (char)optopt) : std::string(long_options[index].name))
                  << ::std::endl;
        exit(1);
      case 'h':
        ::std::cout << usage() << "\n\n" << help() << std::endl;
        exit(0);
      case 'U':
        ::std::cout << usage() << "\nUse --help for more information." << std::endl;
        exit(0);
      case 'V':
        print_version();
        exit(0);
      case '?':
        ::std::cerr << "Use --usage or --help for some help\n";
        exit(1);
      case 'm':
        mer_given = true;
        mer_arg = conv_int<int>((const char*)optarg, err, false);
        CHECK_ERR(int_t, optarg, "-m, --mer=int")
        break;
      case 'n':
        nb_mers_given = true;
        nb_mers_arg = conv_uint<uint64_t>((const char*)optarg, err, true);
        CHECK_ERR(uint64_t, optarg, "-n, --nb-mers=uint64")
        break;
      case 't':
        threads_given = true;
        threads_arg = conv_int<int>((const char*)optarg, err, false);
        CHECK_ERR(int_t, optarg, "-t, --threads=int")
        break;
      case 'o':
        output_given = true;
        output_arg = optarg;
        break;
      case 'q':
        quality_threshold_given = true;
        quality_threshold_arg = conv_uint<uint32_t>((const char*)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-q, --quality-threshold=uint32")
        break;
      case 'c':
        cont_on_low_given = true;
        cont_on_low_arg = conv_uint<uint32_t>((const char*)optarg, err, false);
        CHECK_ERR(uint32_t, optarg, "-c, --cont-on-low=uint32")
        break;
      case 'l':
        min_len_given = true;
        min_len_arg = conv_uint<uint64_t>((const char*)optarg, err, false);
        CHECK_ERR(uint64_t, optarg, "-l, --min-len=uint64")
        break;
      case 'g':
        gzip_flag = true;
        break;
      case 'f':
        false_positive_given = true;
        false_positive_arg = conv_double((const char*)optarg, err, false);
        CHECK_ERR(double_t, optarg, "-f, --false-positive=double")
        break;
      case LOAD_OPT:
        load_given = true;
        load_arg = optarg;
        break;
      }
    }

    // Check that required switches are present
    if(!mer_given)
      error("[-m, --mer=int] required switch");
    if(!nb_mers_given)
      error("[-n, --nb-mers=uint64] required switch");

    // Parse arguments
    if(argc - optind < 0)
      error("Requires at least 0 argument.");
    for( ; optind < argc; ++optind) {
      input_arg.push_back(argv[optind]);
    }
  }
  static const char * usage() { return "Usage: cmdline_parse [options] input:path+"; }
  class error {
    int code_;
    std::ostringstream msg_;

    // Select the correct version (GNU or XSI) version of
    // strerror_r. strerror_ behaves like the GNU version of strerror_r,
    // regardless of which version is provided by the system.
    static const char* strerror__(char* buf, int res) {
      return res != -1 ? buf : "Invalid error";
    }
    static const char* strerror__(char* buf, char* res) {
      return res;
    }
    static const char* strerror_(int err, char* buf, size_t buflen) {
      return strerror__(buf, strerror_r(err, buf, buflen));
    }
    struct no_t { };

  public:
    static no_t no;
    error(int code = EXIT_FAILURE) : code_(code) { }
    explicit error(const char* msg, int code = EXIT_FAILURE) : code_(code)
      { msg_ << msg; }
    error(const std::string& msg, int code = EXIT_FAILURE) : code_(code)
      { msg_ << msg; }
    error& operator<<(no_t) {
      char buf[1024];
      msg_ << ": " << strerror_(errno, buf, sizeof(buf));
      return *this;
    }
    template<typename T>
    error& operator<<(const T& x) { msg_ << x; return (*this); }
    ~error() {
      ::std::cerr << "Error: " << msg_.str() << "\n"
                  << usage() << "\n"
                  << "Use --help for more information"
                  << ::std::endl;
      exit(code_);
    }
  };
  static const char * help() { return
    "Create k-unitigs with support for large k (k>31).\n\n\n\n"
    "Options (default value in (), *required):\n"
    " -m, --mer=int                           *k-mer size\n"
    " -n, --nb-mers=uint64                    *Estimated number of distinct k-mers\n"
    " -t, --threads=int                        Number of threads (1)\n"
    " -o, --output=path                        Ouput file (stdout)\n"
    " -q, --quality-threshold=uint32           Count threshold for high-quality mers (2)\n"
    " -c, --cont-on-low=uint32                 Max length of low quality mer run (0)\n"
    " -l, --min-len=uint64                     Minimum length of k-unitig to output (k+1)\n"
    " -g, --gzip                               Gzip output file. Ignored if -o not given. (false)\n"
    " -f, --false-positive=double              False positive rate in bloom filter (0.01)\n"
    "     --load=path                          Load jellyfish bloom counter\n"
    " -U, --usage                              Usage\n"
    " -h, --help                               This message\n"
    " -V, --version                            Version";
  }
  static const char* hidden() { return ""; }
  void print_version(::std::ostream &os = std::cout) const {
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.0"
#endif
    os << PACKAGE_VERSION << "\n";
  }
  void dump(::std::ostream &os = std::cout) {
    os << "mer_given:" << mer_given << " mer_arg:" << mer_arg << "\n";
    os << "nb_mers_given:" << nb_mers_given << " nb_mers_arg:" << nb_mers_arg << "\n";
    os << "threads_given:" << threads_given << " threads_arg:" << threads_arg << "\n";
    os << "output_given:" << output_given << " output_arg:" << output_arg << "\n";
    os << "quality_threshold_given:" << quality_threshold_given << " quality_threshold_arg:" << quality_threshold_arg << "\n";
    os << "cont_on_low_given:" << cont_on_low_given << " cont_on_low_arg:" << cont_on_low_arg << "\n";
    os << "min_len_given:" << min_len_given << " min_len_arg:" << min_len_arg << "\n";
    os << "gzip_flag:" << gzip_flag << "\n";
    os << "false_positive_given:" << false_positive_given << " false_positive_arg:" << false_positive_arg << "\n";
    os << "load_given:" << load_given << " load_arg:" << load_arg << "\n";
    os << "input_arg:" << vec_str(input_arg) << "\n";
  }
};
#endif // __CMDLINE_PARSE_HPP__"
