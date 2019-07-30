#ifndef _TEST_MISC_H_
#define _TEST_MISC_H_

#include <unistd.h>

struct file_unlink {
  std::string path;
  bool do_unlink;
  explicit file_unlink(const char* s, bool d = true) : path(s), do_unlink(d) { }
  explicit file_unlink(const std::string& s, bool d = true) : path(s), do_unlink(d) { }

  ~file_unlink() {
    if(do_unlink)
      unlink(path.c_str());
  }
};


#endif /* _TEST_MISC_H_ */
