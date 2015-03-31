

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
#include "profiling.h"

#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"

#include "IRKE.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"
#include "argProcessor.hpp"

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
  data.min_any_entropy = 0.4;
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

  int page_size = 1024;

  int num_args = 1;
  try {
       ArgProcessor in_args(narg, args);

       if (in_args.isArgSet("-K")) {
            data.kmer_length = in_args.getIntVal("-K");
	    num_args += 2;
            if(data.me==0) cerr << "Kmer length set to: " << data.kmer_length << endl;
       }

       if (in_args.isArgSet("--minKmerCount")) {
            data.min_kmer_count = in_args.getIntVal("--minKmerCount");
	    num_args += 2;
            if(data.me==0) cerr << "Min kmer coverage set to: " << data.min_kmer_count << endl;
       }

       if (in_args.isArgSet("--minEdgeCount")) {
            data.min_edge_count = in_args.getIntVal("--minEdgeCount");
	    num_args += 2;
            if(data.me==0) cerr << "Min edge coverage set to: " << data.min_edge_count << endl;
       }

       if (in_args.isArgSet("-L")) {
            data.MIN_ASSEMBLY_LENGTH = in_args.getIntVal("-L");
	    num_args += 2;
            if(data.me==0) cerr << "Min assembly length set to: " << data.MIN_ASSEMBLY_LENGTH << endl;
       }

       if (in_args.isArgSet("--min_assembly_coverage")) {
            data.MIN_ASSEMBLY_COVERAGE = in_args.getIntVal("--min_assembly_coverage");
	    num_args += 2;
            if(data.me==0) cerr << "Min assembly coverage set to: " << data.MIN_ASSEMBLY_COVERAGE << endl;
       }
      
       if (in_args.isArgSet("--min_con_ratio")) {
            data.MIN_CONNECTIVITY_RATIO = in_args.getFloatVal("--min_con_ratio");
	    num_args += 2;
       }

       if (in_args.isArgSet("--DS")) {
            data.DS = true;
	    num_args++;
            if(data.me==0) cerr << "double stranded mode set" << endl;
       }

       if (in_args.isArgSet("--min_seed_entropy")) {
            data.MIN_SEED_ENTROPY = in_args.getFloatVal("--min_seed_entropy");
	    num_args += 2;
            if(data.me==0) cerr << "Min seed entropy set to: " << data.MIN_SEED_ENTROPY << endl;
       }

       if (in_args.isArgSet("--min_seed_coverage")) {
            data.MIN_SEED_COVERAGE = in_args.getIntVal("--min_seed_coverage");
	    num_args += 2;
            if(data.me==0) cerr << "min seed coverage set to: " << data.MIN_SEED_COVERAGE << endl;
       }

       if (in_args.isArgSet("--min_any_entropy")) {
            data.min_any_entropy = in_args.getFloatVal("--min_any_entropy");
	    num_args += 2;
            if(data.me==0) cerr << "min entropy set to: " << data.min_any_entropy << endl;
       }

       if (in_args.isArgSet("--no_prune_error_kmers")) {
            data.prune_error_kmers = false;
	    num_args++;
       }

       if (data.prune_error_kmers && in_args.isArgSet("--min_ratio_non_error")) {
            data.min_ratio_non_error = in_args.getFloatVal("--min_ratio_non_error");
	    num_args += 2;
            if(data.me==0) cerr << "Set to prune kmers below min ratio non-erro: " << data.min_ratio_non_error << endl;
       }

       if (in_args.isArgSet("--coverage_outfile")) {
            data.WRITE_COVERAGE = true;
            data.COVERAGE_OUTPUT_FILENAME = in_args.getStringVal("--coverage_outfile");
	    num_args += 2;
       }

       if(in_args.isArgSet("--PageSize")) {
            page_size = in_args.getIntVal("--PageSize");
            num_args += 2;
            if(data.me==0) cerr << "Page size for map reduce object set to: " << page_size << endl << endl;
       }
 
  }

  catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
  }


  data.seed = 123456789;
  srand48(data.seed+data.me);

  int pbits = 0;
  while ((1 << pbits) < data.nprocs) pbits++;
  data.pshift = 63 - pbits;
  int hbits = pbits + 1;
  data.lmask = ALLBITS >> hbits;

  data.nthresh = 1000;

  MapReduce *mrKmers = new MapReduce(MPI_COMM_WORLD);
  mrKmers->memsize = page_size;
  mrKmers->verbosity = 1;
  mrKmers->timer = 0;

  MapReduce *mrE = new MapReduce(MPI_COMM_WORLD);
  mrE->memsize = page_size;
  mrE->verbosity = 1;
  mrE->timer = 0;

  MapReduce *mrV = new MapReduce(MPI_COMM_WORLD);
  mrV->memsize = page_size;
  mrV->verbosity = 1;
  mrV->timer = 0; 

  MapReduce *mrZ = new MapReduce(MPI_COMM_WORLD);
  mrZ->memsize = page_size;
  mrZ->verbosity = 1;
  mrZ->timer = 0;

// ###########################################################

  int VmSm = 1;
  int iter = 0;

  MPI_Barrier(MPI_COMM_WORLD);

  double tstart = MPI_Wtime();  

  int nkmers = mrKmers->map(narg-num_args,&args[num_args],0,1,0,fileread_RNAseq,&data);

  iter++;
  getMemory(iter, data.me, data.nprocs, VmSm); 

  int nfiles = mrKmers->mapfilecount;
  mrKmers->collate(NULL);

  iter++;
  getMemory(iter, data.me, data.nprocs, VmSm);

  data.flag = 0;
  mrKmers->reduce(reduce_kmers_RNAseq,&data);

  iter++;
  getMemory(iter, data.me, data.nprocs, VmSm);

  double tstop = MPI_Wtime();

  unsigned long long flagall = 0;
  MPI_Allreduce(&data.flag,&flagall,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  if(data.me == 0) cerr <<  "Number of kmers =  " << flagall << " Time took for counting kmers = " << tstop - tstart << endl << endl;

  tstart = MPI_Wtime();

  mrE->map(narg-num_args,&args[num_args],0,1,0,fileread_RNAseq_map_Edge,&data);
  nfiles = mrE->mapfilecount;

  iter++;
  getMemory(iter, data.me, data.nprocs, VmSm);

  mrE->collate(NULL);

  iter++;
  getMemory(iter, data.me, data.nprocs, VmSm);

  data.flag = 0;
  mrE->reduce(reduce_Edge_from_RNAseq,&data);

  iter++;
  getMemory(iter, data.me, data.nprocs, VmSm);

  flagall = 0;
  MPI_Allreduce(&data.flag,&flagall,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);

  tstop = MPI_Wtime();
  if(data.me == 0) cerr << "number of edges = " << flagall << " Time took for all possible connections of kmers = " << tstop - tstart << endl << endl;

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

  data.forward_direction = true;
  mrE->map(mrE, map_kmer_edge, &data);
  mrE->collate(NULL);
  mrE->reduce(reduce_keep_dominantEdge, NULL);

  data.forward_direction = false;
  mrE->map(mrE, map_kmer_edge, &data);
  mrE->collate(NULL);
  mrE->reduce(reduce_keep_dominantEdge, NULL);


  tstop = MPI_Wtime();
  if(data.me == 0) cerr << "Time took for filtering edges = " << tstop - tstart << endl << endl;

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

    flagall = 0;
    MPI_Allreduce(&data.flag,&flagall,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    if (flagall == 0) break;

    mrV->map(mrV, map_invert_multi, &data);
    mrV->map(mrZ, map_zone_multi, &data, 1);

    mrV->collate(NULL);
    mrV->reduce(reduce_zone_reassign,&data);

   if(data.me == 0) cerr <<  niterates << " th iteration swithed the number of " << flagall << " zones" <<endl << endl;

  } 

  mrZ->map(mrV,map_strip,NULL);
  mrZ->add(mrKmers);

  mrZ->collate(NULL);
  data.flag = 0;
  mrZ->reduce(reduce_zone_kmer_count,&data);

  flagall = 0;
  MPI_Allreduce(&data.flag,&flagall,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);

  mrZ->collate(NULL);

  tstop = MPI_Wtime();
  if(data.me == 0) {
	cerr << "Total number of kmers with zoneID after clustering = " << flagall << endl;
	cerr << "Time took for clustering of kmers using connected component finding algorithms = " << tstop - tstart << endl;
	cerr << "-----------------------------------------------------------" << endl;
  }

  tstart = MPI_Wtime();
  stringstream out_filename;
  out_filename << "iworm_" << data.me << ".fa";
  data.outFile.open(out_filename.str().c_str());
  
  data.flag = 0;
  mrZ->reduce(reduce_run_inchworm, &data);
  data.outFile.close();
  tstop = MPI_Wtime();
  if(data.me == 0) cerr << "Time took for parallel iworm contigs construction = " << tstop - tstart << endl;

  flagall = 0;
  MPI_Allreduce(&data.flag,&flagall,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if(data.me == 0) cerr << "number of inchworm contigs = " << flagall << endl;

  MPI_Barrier(MPI_COMM_WORLD);

  delete mrKmers;
  delete mrE;
  delete mrV;
  delete mrZ;

  MPI_Finalize();

}

// ########################################################


