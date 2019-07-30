#include <sys/mman.h>
#include <fcntl.h>           /* For O_* constants */
#include <sys/stat.h>        /* For mode constants */
#include <semaphore.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdlib.h>

#include <algorithm>
#include <functional>
#include <memory>

#include <src/semaphore_cmdline.hpp>

class semaphore_type {
protected:
  sem_t*      _sem;

public:
  semaphore_type(sem_t* sem) :
    _sem(sem)
  { }

  virtual ~semaphore_type() {
    sem_close(_sem);
  }

  bool good() { return _sem != SEM_FAILED; }
  void wait() {
    while(true) {
      int res = sem_wait(_sem);
      if(res == 0 || errno != EINTR) return;
    }
  }
  void post() {
    sem_post(_sem);
  }

  int value() {
    int sval = 0;
    sem_getvalue(_sem, &sval);
    return sval;
  }

  virtual void unlink() = 0;
};

class shared_semaphore : public semaphore_type {
  typedef semaphore_type super;
  const std::string _path;

public:
  shared_semaphore(const std::string& path, unsigned int value) :
    super(sem_open(path.c_str(), O_CREAT, S_IRUSR|S_IWUSR, value)),
    _path(path)
  { }

  virtual ~shared_semaphore() { }

  virtual void unlink() {
    sem_unlink(_path.c_str());
  }
};

template<typename T>
class auto_unmap {
  T* _val;
public:
  auto_unmap(int fd) :
    _val((T*)mmap(0, sizeof(T), PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0))
  { }
  ~auto_unmap() {
    if(_val != MAP_FAILED)
      ::munmap(_val, sizeof(T));
  }

  T* get() { return _val; }
  T* release() {
    T* res = _val;
    _val = (T*)MAP_FAILED;
    return res;
  }
};

class file_semaphore : public semaphore_type {
  typedef semaphore_type super;
  const std::string _path;

  static std::pair<int, sem_t*> sem_file_creat(const char* path, unsigned int value) {
    const char* const prefix     = "tmp_file_semaphore_XXXXXX";
    const char*       last_slash = strrchr(path, '/');

    std::unique_ptr<char[]> tmp_sem;
    if(!last_slash) {
      // No slash -> in current directory
      tmp_sem.reset(new char[strlen(prefix) + 1]);
      strcpy(tmp_sem.get(), prefix);
    } else {
      ssize_t len = last_slash - path;
      tmp_sem.reset(new char[len + 1 + strlen(prefix) + 1]);
      strncpy(&tmp_sem[0], path, len);
      tmp_sem[len] = '/';
      strcpy(&tmp_sem[len + 1], prefix);
    }

    int fd = mkstemp(tmp_sem.get());
    if(fd == -1)
      return std::make_pair(-1, SEM_FAILED);
    // Make sure the temporary file is unlinked in all cases
    std::unique_ptr<const char, int(*)(const char*)> unlink_tmp(tmp_sem.get(), ::unlink);

    {
      sem_t sem;
      memset(&sem, '\0', sizeof(sem));
      if(write(fd, &sem, sizeof(sem)) != sizeof(sem))
        return std::make_pair(-1, SEM_FAILED);
    }

    // Memory map and make sure it is unmapped, unless we succeed in creating the semaphore
    auto_unmap<sem_t> ptr(fd);
    close(fd);
    if(ptr.get() == MAP_FAILED)
      return std::make_pair(-1, SEM_FAILED);
    if(sem_init(ptr.get(), 1, value) == -1)
      return std::make_pair(-1, SEM_FAILED);

    if(link(tmp_sem.get(), path) == -1) {
      if(errno == EEXIST)
        return std::make_pair(0, SEM_FAILED);
      else
        return std::make_pair(-1, SEM_FAILED);
    }

    return std::make_pair(0, ptr.release());
  }

  static std::pair<int, sem_t*> sem_file_try(const char* path) {
    int fd = open(path, O_RDWR);
    if(fd == -1) {
      if(errno == ENOENT)
        return std::make_pair(0, SEM_FAILED);
      else
        return std::make_pair(-1, SEM_FAILED);
    }

    sem_t* sem = (sem_t*)mmap(0, sizeof(sem_t), PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);
    if(sem == MAP_FAILED)
      return std::make_pair(-1, SEM_FAILED);

    return std::make_pair(0, sem);
  }

  static sem_t* sem_file_open(const char* path, unsigned int value) {
    std::pair<int, sem_t*> res;
    while(true) {
      res = sem_file_try(path);
      if(res.first == -1 || res.second != SEM_FAILED) return res.second;
      res = sem_file_creat(path, value);
      if(res.first == -1 || res.second != SEM_FAILED) return res.second;
    }
  }

public:
  file_semaphore(const std::string& path, unsigned int value) :
    super(sem_file_open(path.c_str(), value)),
    _path(path)
  { }

  virtual ~file_semaphore() {
    ::munmap(super::_sem, sizeof(sem_t));
  }

  virtual void unlink() {
    ::unlink(_path.c_str());
  }
};

int run_command(std::vector<const char*>& cmd) {
  cmd.push_back((const char*)0);
  pid_t child = fork();
  if(child == -1) return -1;
  if(child == 0) {
    execvp(cmd[0], (char* const*)cmd.data());
    exit(EXIT_FAILURE);
  }

  int status;
  wait(&status);
  return status;
}

int main(int argc, char* argv[]) {
  args_t args(argc, argv);

  if(args.wait_flag || args.print_flag) {
    if(args.cmd_arg.size() > 0 || args.command_given)
      args_t::error() << "No command is allowed with the --wait and --print flag";
  } else {
    if(args.cmd_arg.size() == 0 && !args.command_given)
      args_t::error("Expected a command to run");
    if(args.cmd_arg.size() > 0 && args.command_given)
      args_t::error() << "Pass a command with the -c switch or through the command line";
  }

  std::unique_ptr<semaphore_type> sem;

  if(args.file_given) {
    sem.reset(new file_semaphore(args.file_arg, args.jobs_arg));
  } else {
    std::string sem_path("/");
    if(args.id_given) {
      sem_path += args.id_arg;
    } else {
      const char* tty = ttyname(0);
      if(!tty)
        args_t::error() << "Not connected to a terminal. Must pass an id.";
      sem_path += tty;
    }
    // Replace / by _
    std::replace_if(sem_path.begin() + 1, sem_path.end(), std::bind(std::equal_to<char&>(), '/', std::placeholders::_1), '_');
    sem.reset(new shared_semaphore(sem_path, args.jobs_arg));
  }

  if(!sem->good())
    args_t::error() << "Failed to create or access semaphore" << args_t::error::no;

  if(args.wait_flag) {
    // Wait for semaphore to be free and release it
    for(unsigned int i = 0; i < args.jobs_arg; ++i)
      sem->wait();
    sem->unlink();
  } else if(args.print_flag) {
    std::cout << sem->value() << std::endl;
  } else {
    sem->wait();
    // Go into the background and run the command
    pid_t child = fork();
    if(child == -1)
      args_t::error() << "Failed to go into background" << args_t::error::no;
    if(child == 0) {
      int status;
      if(args.command_given)
        status = system(args.command_arg);
      else
        status = run_command(args.cmd_arg);
      sem->post();
      if(status == -1 || !WIFEXITED(status))
        return EXIT_FAILURE;
      return WEXITSTATUS(status);
    }
  }

  return EXIT_SUCCESS;
}
