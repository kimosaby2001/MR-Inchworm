

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/stat.h"

#define __STDC_LIMIT_MACROS
#include "mapreduce.h"
#include "keyvalue.h"
#include "blockmacros.h"
#include "typedefs.h"

#include "fstream"
#include "iostream"
#include "sstream"

#include "vector"
#include "map"
#include "set"
#include "algorithm"
#include "time.h"
#include "math.h"

#include "mrSubroutine.h"

#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"

#include "IRKE.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"


using namespace MAPREDUCE_NS;
using namespace std;


void Execute(const char * command) {
    int ret = system(command);
    if (ret != 0) {
        cout << "COMMAND: " << command << endl;
        cout << "Died with exit code " << ret << endl;
        cout << "Exiting." << endl;
        exit(-1);
    }
}


// ###############################################
//   Main routine
//


int main(int narg, char **args)
{
  MPI_Init(&narg,&args);

  Data data;

  MPI_Comm_rank(MPI_COMM_WORLD,&data.me);
  MPI_Comm_size(MPI_COMM_WORLD,&data.nprocs);

  data.prune_error_kmers = true;
  data.min_ratio_non_error = 0.05f;
  data.min_kmer_count = 1;
  data.min_edge_count = 1;
  data.min_any_entropy = 0.0;
  data.kmer_length = 25;
  data.DS = false;

  data.MAX_RECURSION = 1;
  data.MIN_SEED_ENTROPY = 1.5;
  data.MIN_SEED_COVERAGE = 2;
  
  data.PACMAN = false;
  data.CRAWL = false;
  data.crawl_length = 1; 

  data.MIN_CONNECTIVITY_RATIO = 0.0;
  data.MIN_ASSEMBLY_LENGTH = data.kmer_length;
  data.MIN_ASSEMBLY_COVERAGE = 2;
  data.WRITE_COVERAGE = false;

  data.num_threads = 16;

  data.seed = 123456789;
  srand48(data.seed+data.me);

  data.ncluster = 1;

  int pbits = 0;
  while ((1 << pbits) < data.nprocs) pbits++;
  data.pshift = 63 - pbits;
  int hbits = pbits + 1;
  data.lmask = ALLBITS >> hbits;

  data.nthresh = 1000;


  MapReduce *mrKmers = new MapReduce(MPI_COMM_WORLD);
  mrKmers->memsize = 2048;
  mrKmers->verbosity = 1;
  mrKmers->timer = 1;

  MapReduce *mrE = new MapReduce(MPI_COMM_WORLD);
  mrE->memsize = 2048;
  mrE->verbosity = 1;
  mrE->timer = 1;

  MapReduce *mrV = new MapReduce(MPI_COMM_WORLD);
  mrV->memsize = 2048;
  mrV->verbosity = 1;
  mrV->timer = 1; 

  MapReduce *mrZ = new MapReduce(MPI_COMM_WORLD);
  mrZ->memsize = 2048;
  mrZ->verbosity = 1;
  mrZ->timer = 1;


  MPI_Barrier(MPI_COMM_WORLD);
  double tstart = MPI_Wtime();  

  int nkmers = mrKmers->map(narg-1,&args[1],0,1,0,fileread_RNAseq,&data);
//  int nkmers = mrKmers->map(narg-1,&args[1],0,1,0,fileread,&data);
  int nfiles = mrKmers->mapfilecount;

  mrKmers->collate(NULL);

  data.flag = 0;
  mrKmers->reduce(reduce_kmers_RNAseq,&data);
//  mrKmers->reduce(reduce_merge_kmers,&data);

  MPI_Barrier(MPI_COMM_WORLD);
  double tstop = MPI_Wtime();

  int flagall;
  MPI_Allreduce(&data.flag,&flagall,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if(data.me == 0) cerr <<  "Number of kmers =  " << flagall << " Time took for counting kmers = " << tstop - tstart << endl << endl;

  tstart = MPI_Wtime();
//  mrE->map(mrKmers, map_kmers_edge_vertex, &data);
//  mrE->collate(NULL);
//  mrE->reduce(reduce_filter_nonExisting_edge,NULL);

  mrE->map(narg-1,&args[1],0,1,0,fileread_RNAseq_map_Edge,&data);
  mrE->collate(NULL);
  mrE->reduce(reduce_Edge_from_RNAseq,&data);

  tstop = MPI_Wtime();
  if(data.me == 0) {
	cerr << "Time took for all possible connections of kmers = " << tstop - tstart << endl;
	cerr << "-----------------------------------------------------------" << endl;
  }

  tstart = MPI_Wtime();
  data.forward_direction = true;
  mrE->map(mrE, map_kmer_edge, &data);
  mrE->add(mrKmers);
  mrE->collate(NULL);
  mrE->reduce(reduce_count_to_edge, &data);

  data.forward_direction = false;
  mrE->map(mrE, map_kmer_edge, &data);
  mrE->add(mrKmers);
  mrE->collate(NULL);
  mrE->reduce(reduce_count_to_edge, &data);

//  data.forward_direction = true;
//  mrE->map(mrE, map_kmer_edge, &data);
//  mrE->collate(NULL);
//  mrE->reduce(reduce_filter_edge, &data);

//  data.forward_direction = false;
//  mrE->map(mrE, map_kmer_edge, &data);
//  mrE->collate(NULL);
//  mrE->reduce(reduce_filter_edge, &data);

  mrE->map(mrE, map_filter_edge_by_count, &data);

  tstop = MPI_Wtime();
  if(data.me == 0) {
        cerr << "Time took for filtering edges = " << tstop - tstart << endl;
        cerr << "---------------------------------------------------------" << endl;
  }

  tstart = MPI_Wtime();
  mrV->map(mrE,edge_to_vertices,NULL);
  mrV->collate(NULL);
  mrV->reduce(reduce_self_zone,NULL);
  
  int niterates = 0;

  while(1) {

    niterates++;

    mrZ->map(mrE,map_edge_vert,NULL);
    mrZ->add(mrV);
    mrZ->collate(NULL);
    mrZ->reduce(reduce_edge_zone,NULL);

    mrZ->collate(NULL);
    data.flag = 0;
    mrZ->reduce(reduce_zone_winner,&data);
    flagall;
    MPI_Allreduce(&data.flag,&flagall,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    if (flagall == 0) break;

    mrV->map(mrV, map_invert_multi, &data);
    mrV->map(mrZ, map_zone_multi, &data, 1);
    mrV->collate(NULL);
    mrV->reduce(reduce_zone_reassign,&data);

   if(data.me == 0) 
	cerr <<  niterates << " th iteration swithed the number of " << flagall << " zones" <<endl << endl;

  } 

  if(data.me == 0)  cerr << "-----------------------------------------------------------" << endl;

  mrZ->map(mrV,map_strip,NULL);

  mrZ->add(mrKmers);
  mrZ->collate(NULL);

  data.flag = 0;
  mrZ->reduce(reduce_zone_kmer_count,&data);

  unsigned long long flagall_ull = 0;
  MPI_Allreduce(&data.flag,&flagall_ull,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);

  mrZ->collate(NULL);

  tstop = MPI_Wtime();
  if(data.me == 0) {
	cerr << "Total number of kmers with zoneID after clustering = " << flagall_ull << endl;
	cerr << "Time took for clustering of kmers using connected component finding algorithms = " << tstop - tstart << endl;
	cerr << "-----------------------------------------------------------" << endl;
  }

  //data.flag = 0;
  //flagall_ull = 0;

  //stringstream out_filename;
  //out_filename << "kmer_" << data.me << ".fa";
  //data.outFile.open(out_filename.str().c_str());
  //mrZ->reduce(reduce_print_clustered_kmers, &data);
  //data.outFile.close();

  //MPI_Allreduce(&data.flag,&flagall_ull,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  //if(data.me == 0) cerr << "Total number of kmers after clustering = " << flagall_ull << endl;


  tstart = MPI_Wtime();
  stringstream out_filename;
  out_filename << "iworm_" << data.me << ".fa";
  data.outFile.open(out_filename.str().c_str());
  
  data.flag = 0;
  mrZ->reduce(reduce_run_inchworm, &data);
  data.outFile.close();
  tstop = MPI_Wtime();
  if(data.me == 0) cerr << "Time took for parallel iworm contigs construction = " << tstop - tstart << endl;

  flagall_ull = 0;
  MPI_Allreduce(&data.flag,&flagall_ull,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if(data.me == 0) cerr << "number of inchworm contigs = " << flagall_ull << endl;

  MPI_Barrier(MPI_COMM_WORLD);

  delete mrKmers;
  delete mrE;
  delete mrV;
  delete mrZ;

  MPI_Finalize();

}

// ########################################################


