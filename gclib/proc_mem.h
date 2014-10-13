#ifndef PROC_MEM_H
#define PROC_MEM_H
#include <stdio.h>
// a Linux-specific way to report the memory usage of the current process
void get_mem_usage(double& vm_usage, double& resident_set);

void print_mem_usage(FILE* fout=stderr);

#endif
