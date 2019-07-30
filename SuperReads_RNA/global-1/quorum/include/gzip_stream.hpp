#ifndef __GZIP_STREAM__
#define __GZIP_STREAM__

#include <iostream>
#include <stdio.h>
#include <ext/stdio_filebuf.h>
#include <stdexcept>

template<typename _CharT, typename _Traits = std::char_traits<_CharT> >
class basic_gzipstream : public std::ostream {
  typedef __gnu_cxx::stdio_filebuf<_CharT> stdbuf;
public:
  basic_gzipstream(const char *filename) : std::ostream(open_gzip(filename)), closed(false) { }
  virtual ~basic_gzipstream() {
    close();
    delete rdbuf();
  }
  void close() {
    if(!closed) {
      flush();
      pclose(((stdbuf*)rdbuf())->file());
      closed = true;
    }
  }

private:
  static stdbuf *open_gzip(const char *filename) {
    std::string command("gzip -1 > '");
    command += filename;
    command += "'";
    FILE *f = popen(command.c_str(), "w");
    if(!f)
      throw std::runtime_error("popen failed");
    return new stdbuf(f, std::ios::out);
  }
  bool closed;
};

typedef basic_gzipstream<char> gzipstream;

#endif
