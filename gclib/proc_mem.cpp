#include "proc_mem.h"
#ifdef __APPLE__
#include<mach/mach.h>
void get_mem_usage(double& vm_usage, double& resident_set) {
  vm_usage=0;
  resident_set=0;
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS == task_info(mach_task_self(),
         TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count))

    {
    vm_usage=double(t_info.virtual_size)/1024;
    resident_set=double(t_info.resident_size)/1024;
    }
// resident size is in t_info.resident_size;
// virtual size is in t_info.virtual_size;
}

#elif defined(_WIN32) || defined(_WIN64)
#include "windows.h"
#include "psapi.h"
void get_mem_usage(double& vm_usage, double& resident_set) {
	PROCESS_MEMORY_COUNTERS_EX pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), (PPROCESS_MEMORY_COUNTERS)&pmc, sizeof(pmc));
	//SIZE_T virtualMemUsedByMe = pmc.PrivateUsage;
    //SIZE_T physMemUsedByMe = pmc.WorkingSetSize;
    vm_usage=(double)pmc.PrivateUsage;
    resident_set=(double)pmc.WorkingSetSize;
}

#else //assume Linux
#include <unistd.h>
#include <string>
#include <ios>
#include <fstream>
void get_mem_usage(double& vm_usage, double& resident_set) {
   using std::ios_base;
   using std::ifstream;
   using std::string;
   vm_usage     = 0.0;
   resident_set = 0.0;
   // 'file' stat seems to give the most reliable results
   ifstream stat_stream("/proc/self/stat",ios_base::in);
   // dummy vars for leading entries in stat that we don't care about
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}
#endif

//////////////////////////////////////////////////////////////////////////////
// get_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0


void print_mem_usage(FILE* fout) {
  double vs, rs;
  get_mem_usage(vs,rs);
  vs/=1024;
  rs/=1024;
  fprintf(fout, "VMSize: %6.1fMB\tRSize: %6.1fMB\n", vs, rs);
  }
