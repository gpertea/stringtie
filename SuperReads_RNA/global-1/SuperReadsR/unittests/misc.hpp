#ifndef __MISC_HPP__
#define __MISC_HPP__

#include <cstring>
#include <stdexcept>
#include <pthread.h>

void pdo(unsigned int n, void* (*f)(void*), void *data);

#endif /* __MISC_HPP__ */
