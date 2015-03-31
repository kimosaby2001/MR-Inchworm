
#include "mpi.h"

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "iostream"

#include <sys/sysinfo.h>

using namespace std;

void printMemory(int iter, long long *Mem1, long long *Mem2, int n);
void getMemory(int iter, int rank, int nprocs, int VmSm);
long long getTotalVirtualMemory();
long long getFreeVirtualMemory();
long long getTotalSystemMemory();
long long getFreeSystemMemory();

