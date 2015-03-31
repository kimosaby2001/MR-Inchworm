
#include "profiling.h"

using namespace std;

void printMemory(int iter, long long *Mem1, long long *Mem2, int n)
{
  cerr << iter;
  for(int i=0; i<n; i++) cerr << "\t" << Mem1[i] -  Mem2[i];
  cerr << endl;
}

void getMemory(int iter, int rank, int nprocs, int VmSm)
{
   long long totalMem, freeMem;

   if(VmSm == 0) {
     totalMem = getTotalVirtualMemory();
     freeMem  = getFreeVirtualMemory();
   } else {
     totalMem = getTotalSystemMemory();
     freeMem  = getFreeSystemMemory();
   }

   long long *tBytes_1 =  (long long *) malloc(sizeof(long long) * nprocs);
   long long *tBytes_2 =  (long long *) malloc(sizeof(long long) * nprocs);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Gather(&totalMem, 1, MPI_LONG_LONG, tBytes_1, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
   MPI_Gather(&freeMem,  1, MPI_LONG_LONG, tBytes_2, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);
   if(rank == 0) printMemory( iter, tBytes_1, tBytes_2, nprocs );

   free(tBytes_1);
   free(tBytes_2);
}

long long getTotalVirtualMemory()
{
    struct sysinfo memInfo;
    
    sysinfo (&memInfo);
    long long totalVirtualMem = memInfo.totalram;
    totalVirtualMem += memInfo.totalswap;
    totalVirtualMem *= memInfo.mem_unit;
    return totalVirtualMem;
}

long long getFreeVirtualMemory()
{
    struct sysinfo memInfo;

    sysinfo (&memInfo);
    long long totalVirtualMem = memInfo.freeram;
    totalVirtualMem += memInfo.freeswap;
    totalVirtualMem *= memInfo.mem_unit;
    return totalVirtualMem;
}


long long getTotalSystemMemory()
{
  struct sysinfo memInfo;

  sysinfo (&memInfo);
  long long totalVirtualMem = memInfo.totalram;
  totalVirtualMem *= memInfo.mem_unit;
  return totalVirtualMem;
}


long long getFreeSystemMemory()
{
  struct sysinfo memInfo;

  sysinfo (&memInfo);
  long long totalVirtualMem = memInfo.freeram;
  totalVirtualMem *= memInfo.mem_unit;
  return totalVirtualMem;

}



