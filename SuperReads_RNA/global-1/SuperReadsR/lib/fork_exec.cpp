#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <errno.h>
#include <stdarg.h>
#include <string.h>

#include <vector>

int child_exec(int fd, const char *file, char *const argv[]) {
  while(close(fd) == -1)
    if(errno != EINTR)
      return -1;

  execvp(file, argv);
  return -1;
}

int parent_read(int fd_to_close, int fd) {
  while(close(fd_to_close) == -1)
    if(errno != EINTR)
      return -1;

  ssize_t bytes_read = 0;
  while((bytes_read = read(fd, &errno, sizeof(errno))) == -1)
    if(errno != EINTR)
        return -1;

  return bytes_read == 0 ? 0 : -1;
}

pid_t fork_execvp(const char* file, char* const argv[]) {
  int     err[2];
  pid_t   pid;
  ssize_t werr = 0;

  if(pipe(err) == -1)
    goto error_none_close;
  if(fcntl(err[1], F_SETFD, FD_CLOEXEC) == -1)
    goto error_both_close;

  switch((pid = fork())) {
  case -1: goto error_both_close;
  case 0:
    child_exec(err[0], file, argv);
    /* Get here only if failure. If write failes, then we are SOL: no
       way to return error to parent! */
    do {
      werr = write(err[1], &errno, sizeof(errno));
      if(werr == -1 && errno == EINTR)
        continue;
    } while(0);
    close(err[1]);
    exit(0);

  default:
    if(parent_read(err[1], err[0]) == -1)
      goto error_one_close;
  }

  close(err[0]);
  return pid;

 error_both_close:
  close(err[1]);
 error_one_close:
  close(err[0]);
 error_none_close:
  return -1;
}

pid_t fork_execlp(const char* file, const char* arg, ...) {
  va_list ap;
  char* x;
  std::vector<char*> argv;

  va_start(ap, arg);
  do {
    x = va_arg(ap, char*);
    argv.push_back(x);
  } while(x != 0);
  va_end(ap);

  return fork_execvp(file, argv.data());
}

pid_t fork_exec(char* cmd) {
  char* saveptr = 0;
  char* current = strtok_r(cmd, " ", &saveptr);
  std::vector<char*> argv;

  for( ; current; current = strtok_r(0, " ", &saveptr))
    argv.push_back(current);
  argv.push_back((char*)0);

  return fork_execvp(argv[0], argv.data());
}

struct dup_str{
  char* data;
  dup_str(const char* str) : data(strdup(str)) { }
  ~dup_str() { free(data); }
};

pid_t fork_exec(const char* cmd) {
  dup_str cmd_duped(cmd);
  if(!cmd_duped.data)
    return -1;
  return fork_exec(cmd_duped.data);
}
