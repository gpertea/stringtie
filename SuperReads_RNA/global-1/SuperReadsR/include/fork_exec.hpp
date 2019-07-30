/** @file fork_exec.hpp
 * Defines a family of functions: fork_exec*. They behave like
 * calling fork followed by exec. The return value and errno reflect
 * the success of both the fork and exec system call (even though
 * technically the exec system call was performed by a child
 * process). In other word, if fork_exec* returns a valid PID
 * (i.e. not -1), then it is guaranteed that the child process is
 * executing. The status returned by `wait` or `waitpid` is whatever
 * is returned by the program executed.
 */

/**
 * Fork and exec, like execvp. Return the pid of the child process in
 * case of success or -1 in case of failure (of either fork, execvp or
 * some other syscall: close, pipe). Errno is set properly. See the
 * man page for `exec`.
 *
 * @param file Path to executable.
 * @param argv Null terminated array of arguments (C strings)
 * @return -1 on error. Pid of child process on success.
 */
pid_t fork_execvp(const char* file, char *const argv[]);

/**
 * Fork and exec, like execlp. 
 *
 * @param file Path to executable
 * @param arg 0th argument (executable name), followed by a 0 terminated list of arguments.
 * @return -1 on error. Pid of child process on success.
 * @see fork_execvp.
 */
pid_t fork_execlp(const char* file, const char* arg, ...);

/**
 * Split command on white spaces, then fork and exec. Warning: this
 * command modifies its input argument.
 *
 * @param cmd Command to execute.
 * @return -1 on error. Pid of child process on success.
 * @see fork_execvp
 */
pid_t fork_exec(char* cmd);

/**
 * Split command on white spaces, then fork and exec. @see fork_execvp.
 */
pid_t fork_exec(const char* cmd);
