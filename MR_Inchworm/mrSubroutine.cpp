
#include "IRKE.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"


#include "typedefs.h"
#include "mrSubroutine.h"
#include "profiling.h"

using namespace MAPREDUCE_NS;
using namespace std;


// ###################################################################
void map_strip(uint64_t itask, char *key, int keybytes,
               char *value, int valuebytes, KeyValue *kv, void *ptr)
{
  uint64_t zone = *(uint64_t *) value;
  zone &= INT64MAX;
  kv->add(key,keybytes,(char *) &zone,sizeof(uint64_t));
}

// #####################################################################
void reduce_zone_kmer_count(char *key, int keybytes,
                            char *multivalue, int nvalues,
                            int *valuebytes, KeyValue *kv, void *ptr)
{
    int check = 0;
    Data *data = (Data *) ptr;
    int i;

    char *value;
    uint64_t zone;
    KmerCount kmer_count;
    kmer_count.kmer = *(kmer_int_type_t *) key;
    kmer_count.count=0;

    uint64_t nvalues_total;
    CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
    BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

     value = multivalue;
     for(i=0; i<nvalues; i++) {
       if(valuebytes[i] == sizeof(uint64_t))
       {
	 check++;
	 zone = *(uint64_t *) value;
	 break;
       }
       value += valuebytes[i];
     }

     value = multivalue;
     for(i=0; i<nvalues; i++) {
       if(valuebytes[i] == sizeof(unsigned int))
       {
	 check++;
	 kmer_count.count = *(unsigned int *) value;
	 break;
       }
       value += valuebytes[i];
     }

    END_BLOCK_LOOP

    if((check == 2) && (kmer_count.count > 0)) {
      data->flag++;
      kv->add((char *) &zone, sizeof(uint64_t), (char *) &kmer_count, sizeof(KmerCount));
    }

}

// ##############################################################################

void reduce_run_inchworm(char *key, int keybytes,
                         char *multivalue, int nvalues,
                         int *valuebytes, KeyValue *kv, void *ptr)
{
 if(nvalues >= 0) {
    Data *data = (Data *) ptr;
    IRKE irke(data->kmer_length, data->MAX_RECURSION, data->MIN_SEED_ENTROPY, data->MIN_SEED_COVERAGE,
              data->min_any_entropy, data->PACMAN, data->CRAWL, data->crawl_length, data->DS);

    KmerCount kmer_count;

    uint64_t zone = *(uint64_t *) key;
    char *value;
    uint64_t nvalues_total;
    CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
    BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    value = multivalue;
    for(int i=0; i<nvalues; i++) {
        kmer_count = *(KmerCount *) value;
        irke.add_kmer(kmer_count.kmer, kmer_count.count);
        value += valuebytes[i];

    }

    END_BLOCK_LOOP

//if(data->me == 5) cerr << "kcounter size = " << irke.get_size() << endl;

  if(irke.get_size() > 0) {

    if (data->min_kmer_count > 1 || data->min_any_entropy > 0 || data->prune_error_kmers)
      irke.prune_some_kmers(data->min_kmer_count, data->min_any_entropy, data->prune_error_kmers, data->min_ratio_non_error);

      irke.populate_sorted_kmers_list();

     data->flag += irke.compute_sequence_assemblies( (unsigned long long) zone, data->outFile, 
  						    data->MIN_CONNECTIVITY_RATIO, data->MIN_ASSEMBLY_LENGTH,
                                             	    data->MIN_ASSEMBLY_COVERAGE, data->WRITE_COVERAGE, 
						    data->COVERAGE_OUTPUT_FILENAME );

  }

 }

}

// ##############################################################################
void reduce_print_clustered_kmers(char *key, int keybytes,
                                  char *multivalue, int nvalues,
                                  int *valuebytes, KeyValue *kv, void *ptr)
{
//  if(nvalues >= 0) {

    char *value;
    Data *data = (Data *) ptr;
    KmerCount kmer_count;

    uint64_t zone = *(uint64_t *) key;
    uint64_t nvalues_total;
    CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
    BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    value = multivalue;
    for(int i=0; i<nvalues; i++) {
        kmer_count = *(KmerCount *) value;
        string seq = decode_kmer_from_intval(kmer_count.kmer, data->kmer_length);

        //data->outFile << ">" << zone << "_" << i << "_" << kmer_count.count << endl << seq << endl;
        value += valuebytes[i];
	data->flag++;

	kv->add(key, keybytes, (char *) &kmer_count.kmer, sizeof(kmer_int_type_t));
    }

    END_BLOCK_LOOP
//  }
}


// ##################################################################################
void map_invert_multi(uint64_t itask, char *key, int keybytes,
                      char *value, int valuebytes,
                      KeyValue *kv, void *ptr)
{
  uint64_t z = *(uint64_t *) value;

  if (z >> 63) {
    Data *data = (Data *) ptr;
    uint64_t iproc = static_cast<uint64_t> (data->nprocs * drand48());
    uint64_t znew = z | (iproc << data->pshift);
    kv->add((char *) &znew,sizeof(uint64_t),key,keybytes);
  } else kv->add(value,valuebytes,key,keybytes);
}


// ##################################################################################
void map_zone_multi(uint64_t itask, char *key, int keybytes,
                    char *value, int valuebytes,
                    KeyValue *kv, void *ptr)
{
  uint64_t z = *(uint64_t *) key;

  if (z >> 63) {
    Data *data = (Data *) ptr;
    uint64_t zstrip = z & INT64MAX;
    kv->add((char *) &zstrip,sizeof(uint64_t),value,valuebytes);
    int nprocs = data->nprocs;
    int pshift = data->pshift;
    uint64_t znew;
    for (uint64_t iproc = 0; iproc < nprocs; iproc++) {
      znew = zstrip | (iproc << pshift);
      znew |= HIBIT;
      kv->add((char *) &znew,sizeof(uint64_t),value,valuebytes);
    }
  } else kv->add(key,keybytes,value,valuebytes);
}


// ###################################################################################
void reduce_zone_reassign(char *key, int keybytes,
                                  char *multivalue, int nvalues,
                                  int *valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  int nthresh = data->nthresh;
  uint64_t lmask = data->lmask;

  int i;
  int hnew;
  char *value;
  uint64_t znew;

  uint64_t zcount = 0;
  uint64_t zone = *(uint64_t *) key;
  int hkey = zone >> 63;
  zone &= lmask;
  int hwinner = 0;

  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  value = multivalue;
  for (i = 0; i < nvalues; i++) {
    if (valuebytes[i] != sizeof(uint64_t)) {
      znew = *(uint64_t *) value;
      hnew = znew >> 63;
      znew &= INT64MAX;
      if (znew < zone) {
        zone = znew;
        hwinner = hnew;
      }
      zcount++;
    }
    value += valuebytes[i];
  }

  END_BLOCK_LOOP

  if (hkey || hwinner) zone |= HIBIT;
  else if (nvalues_total-zcount > nthresh) zone |= HIBIT;

  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  value = multivalue;
  for (i = 0; i < nvalues; i++) {
    if (valuebytes[i] == sizeof(uint64_t))
      kv->add(value,valuebytes[i],(char *) &zone,sizeof(uint64_t));
    value += valuebytes[i];
  }

  END_BLOCK_LOOP
}


// #########################################################################################
void reduce_zone_winner(char *key, int keybytes,
                        char *multivalue, int nvalues,
                        int *valuebytes, KeyValue *kv, void *ptr)
{
  uint64_t *z = (uint64_t *) multivalue;
  uint64_t z0 = z[0] & INT64MAX;
  uint64_t z1 = z[1] & INT64MAX;

  if (z0 == z1) return;

  Data *data = (Data *) ptr;
  data->flag++;

  PAD *pad = &(data->pad);

  if (z0 > z1) {
    pad->zone = z[1];
    kv->add((char *) &z[0],sizeof(uint64_t),(char *) pad,sizeof(PAD));
  } else {
    pad->zone = z[0];
    kv->add((char *) &z[1],sizeof(uint64_t),(char *) pad,sizeof(PAD));
  }

}


// ###########################################################################################
void map_edge_vert(uint64_t itask, char *key, int keybytes,
                   char *value, int valuebytes,
                   KeyValue *kv, void *ptr)
{
  EDGE *edge = (EDGE *) key;
  kv->add((char *) &edge->vi,sizeof(VERTEX),key,sizeof(EDGE));
  kv->add((char *) &edge->vj,sizeof(VERTEX),key,sizeof(EDGE));
}


// ###########################################################################################
void reduce_edge_zone(char *key, int keybytes,
                      char *multivalue, int nvalues,
                      int *valuebytes, KeyValue *kv, void *ptr)
{
  int i;
  char *value;

  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  value = multivalue;
  for (i = 0; i < nvalues; i++) {
    if (valuebytes[i] == sizeof(uint64_t)) break;
    value += valuebytes[i];
  }
  if (i < nvalues) break;

  END_BLOCK_LOOP

  uint64_t zone = *(uint64_t *) value;

  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  value = multivalue;
  for (i = 0; i < nvalues; i++) {
    if (valuebytes[i] != sizeof(uint64_t))
      kv->add(value,valuebytes[i],(char *) &zone,sizeof(uint64_t));
    value += valuebytes[i];
  }

  END_BLOCK_LOOP
}


// #############################################################################################
void edge_to_vertices(uint64_t itask, char *key, int keybytes, char *value,
                      int valuebytes, KeyValue *kv, void *ptr)
{
  EDGE *edge = (EDGE *) key;
  kv->add((char *) &edge->vi,sizeof(VERTEX),NULL,NULL);
  kv->add((char *) &edge->vj,sizeof(VERTEX),NULL,NULL);
}


// ##############################################################################################
void reduce_self_zone(char *key, int keybytes,
                      char *multivalue, int nvalues,
                      int *valuebytes, KeyValue *kv, void *ptr)
{
  kv->add(key,keybytes,key,keybytes);
}


// ###############################################################################################
void output_edge(uint64_t itask, char *key, int keybytes, char *value,
                 int valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;

  EDGE edge = *(EDGE *) key;
  data->outFile << edge.vi << "\t" << edge.vj << "\t" << edge.ci << "\t" << edge.cj << endl;
  kv->add(key,keybytes, value, valuebytes);
}


// #################################################################################################

void reduce_filter_edge(char *key, int keybytes,
                        char *multivalue, int nvalues,
                        int *valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  char *value;
  int i;
  uint64_t nvalues_total;
  kmer_int_type_t kmer = *(kmer_int_type_t *) key;
 
  float ratio_count; 
  unsigned int dominant_count = data->min_kmer_count;

  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  value = multivalue;
  for (i = 0; i < nvalues; i++) {
	EDGE edge = *(EDGE *) value;
	if(data->forward_direction) {
	    if(edge.cj > dominant_count) dominant_count = edge.cj;
	} else {
	    if(edge.ci > dominant_count) dominant_count = edge.ci;
	}
	value += valuebytes[i];	
  } 

  value = multivalue;
  for (i = 0; i < nvalues; i++) {

	EDGE edge = *(EDGE *) value;
//	if(data->forward_direction)
//	   ratio_count = edge.cj / dominant_count; 
//	else
//	   ratio_count = edge.ci / dominant_count;

//	if(ratio_count >= data->min_ratio_non_error) 
//	  kv->add((char *) &edge, sizeof(EDGE), NULL, NULL);	

	if( data->forward_direction ) {
		if(edge.cj == dominant_count) 
		  kv->add((char *) &edge, sizeof(EDGE), NULL, NULL);
	} else {
		if(edge.ci == dominant_count)
                  kv->add((char *) &edge, sizeof(EDGE), NULL, NULL);
	}	

	value += valuebytes[i];
  }

  END_BLOCK_LOOP

}

// ##################################################################################################
void map_filter_edge_by_count(uint64_t itask, char *key, int keybytes,
                              char *value, int valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  EDGE edge = *(EDGE *) key;
   if(edge.ce >= data->min_edge_count)
	kv->add(key, keybytes, NULL, NULL);
}

// #################################################################################################

void map_kmer_edge (uint64_t itask, char *key, int keybytes,
                    char *value, int valuebytes,
                    KeyValue *kv, void *ptr)
{
   Data *data = (Data *) ptr;
   kmer_int_type_t kmer;
   EDGE edge = *(EDGE *) key;

   if(data->forward_direction) kmer = edge.vi;
   else                        kmer = edge.vj;

   kv->add((char *) &kmer, sizeof(kmer_int_type_t), key, keybytes);
}

// ########################################################################################

void reduce_filter_edge_by_count(char *key, int keybytes,
                                 char *multivalue, int nvalues,
                                 int *valuebytes, KeyValue *kv, void *ptr)
{
   char *value;
   Data *data = (Data *) ptr;
   EDGE edge;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   value = multivalue;
   if(nvalues > 1) {
     for (int i = 0; i < nvalues; i++) {
       edge = *(EDGE *) value;
       if(edge.ce >= data->min_edge_count) kv->add((char *) &edge, sizeof(EDGE), NULL, NULL);
       value += valuebytes[i]; 
     }
   } else {
     edge = *(EDGE *) value;
     kv->add((char *) &edge, sizeof(EDGE), NULL, NULL);
   }
 
   END_BLOCK_LOOP

}

// ########################################################################################

void reduce_keep_dominantEdge(char *key, int keybytes,
                                 char *multivalue, int nvalues,
                                 int *valuebytes, KeyValue *kv, void *ptr)
{
   char *value;
   EDGE edge;
   int i;

   uint64_t nvalues_total;

   if(nvalues > 1) {

    unsigned int dominant_count = 0;

    CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
    BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    value = multivalue;
    for (i = 0; i < nvalues; i++) {
	edge = *(EDGE *) value;
	if(edge.ce > dominant_count) dominant_count = edge.ce;
	value += valuebytes[i];
    }

    value = multivalue;
    for (i = 0; i < nvalues; i++) {
        edge = *(EDGE *) value;
	if(edge.ce == dominant_count) {
	   kv->add((char *) &edge, sizeof(EDGE), NULL, NULL);
	}
	value += valuebytes[i];
    }

    END_BLOCK_LOOP

   } else {

    CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
    BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

    value = multivalue;
    edge = *(EDGE *) value;
    kv->add((char *) &edge, sizeof(EDGE), NULL, NULL);

    END_BLOCK_LOOP

   }

}

// #########################################################################################

void reduce_count_to_edge(char *key, int keybytes,
                          char *multivalue, int nvalues,
                          int *valuebytes, KeyValue *kv, void *ptr)
{
   char *value;
   Data *data = (Data *) ptr;
   int i;

   uint64_t nvalues_total;
   CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
   BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

   value = multivalue;
   for (i=0; i<nvalues; i++) {
        if(valuebytes[i] == sizeof(unsigned int)) break;
        value += valuebytes[i];
   }
   unsigned int count = *(unsigned int *) value;

   value = multivalue;
   for (i = 0; i < nvalues; i++) {
        if(valuebytes[i] == sizeof(EDGE)) {
	  EDGE edge = *(EDGE *) value;
   	  if(data->forward_direction) edge.ci = count;
   	  else                        edge.cj = count;
   	  kv->add((char *) &edge, sizeof(EDGE), NULL, NULL);
	} 
        value += valuebytes[i];
   }

   END_BLOCK_LOOP

}


// ##############################################################################################
//
void map_kmers_edge_vertex(uint64_t itask, char *key, int keybytes,
                           char *value, int valuebytes,
                           KeyValue *kv, void *ptr)
{
  EDGE edge;
  vector<kmer_int_type_t> candidates;

  Data *data = (Data *) ptr;
  kmer_int_type_t kmer = *(kmer_int_type_t *) key;
  unsigned int kmer_length = data->kmer_length;

  edge.ci = 0;
  edge.cj = 0;

  if(data->DS) {

    candidates = get_Fkmer_candidates(kmer, kmer_length);
    for(int i=0; i<candidates.size(); i++) {
        edge.vi = kmer;
        edge.vj = candidates[i];
        kv->add((char *) &edge, sizeof(EDGE), (char *) &kmer, sizeof(kmer_int_type_t));

        edge.vj = kmer;
        edge.vi = candidates[i];
        kv->add((char *) &edge, sizeof(EDGE), (char *) &kmer, sizeof(kmer_int_type_t));
    }

    candidates = get_Rkmer_candidates(kmer, kmer_length);
    for(int i=0; i<candidates.size(); i++) {
        edge.vi = kmer;
        edge.vj = candidates[i];
        kv->add((char *) &edge, sizeof(EDGE), (char *) &kmer, sizeof(kmer_int_type_t));

        edge.vj = kmer;
        edge.vi = candidates[i];
        kv->add((char *) &edge, sizeof(EDGE), (char *) &kmer, sizeof(kmer_int_type_t));
    }

    kmer_int_type_t rev_kmer = revcomp_val(kmer, kmer_length);

    candidates = get_Fkmer_candidates(rev_kmer, kmer_length);
    for(int i=0; i<candidates.size(); i++) {
        edge.vi = kmer;
        edge.vj = candidates[i];
        kv->add((char *) &edge, sizeof(EDGE), (char *) &kmer, sizeof(kmer_int_type_t));

        edge.vj = kmer;
        edge.vi = candidates[i];
        kv->add((char *) &edge, sizeof(EDGE), (char *) &kmer, sizeof(kmer_int_type_t));
    }

    candidates = get_Rkmer_candidates(rev_kmer, kmer_length);
    for(int i=0; i<candidates.size(); i++) {
        edge.vi = kmer;
        edge.vj = candidates[i];
        kv->add((char *) &edge, sizeof(EDGE), (char *) &kmer, sizeof(kmer_int_type_t));

        edge.vj = kmer;
        edge.vi = candidates[i];
        kv->add((char *) &edge, sizeof(EDGE), (char *) &kmer, sizeof(kmer_int_type_t));
    }

  } else {

    edge.vi = kmer;
    candidates = get_Fkmer_candidates(kmer, kmer_length);
    for(int i=0; i<candidates.size(); i++) {
        edge.vj = candidates[i];
        kv->add((char *) &edge, sizeof(EDGE), (char *) &kmer, sizeof(kmer_int_type_t));
    }

    edge.vj = kmer;
    candidates = get_Rkmer_candidates(kmer, kmer_length);
    for(int i=0; i<candidates.size(); i++) {
        edge.vi = candidates[i];
        kv->add((char *) &edge, sizeof(EDGE), (char *) &kmer, sizeof(kmer_int_type_t));
    }

 }

}


// ###########################################################################################
//
void reduce_filter_nonExisting_edge(char *key, int keybytes,
                                    char *multivalue, int nvalues,
                                    int *valuebytes, KeyValue *kv, void *ptr)
{
  if(nvalues > 1) kv->add(key, keybytes, NULL, NULL);
}


// #############################################################################################
//
void reduce_merge_kmers(char *key, int keybytes,
                        char *multivalue, int nvalues,
                        int *valuebytes, KeyValue *kv, void *ptr)
{
  int i;
  char *value;
  unsigned int count;

  Data *data = (Data *) ptr;

  uint64_t nvalues_total;
  CHECK_FOR_BLOCKS(multivalue,valuebytes,nvalues,nvalues_total)
  BEGIN_BLOCK_LOOP(multivalue,valuebytes,nvalues)

  count = 0;
  value = multivalue;
  for (i = 0; i < nvalues; i++) {
     if(valuebytes[i] == sizeof(unsigned int))
        count += *(unsigned int *) value;
     value += valuebytes[i];
  }
  END_BLOCK_LOOP

//  kmer_int_type_t kmer = *(kmer_int_type_t *) key;
//  float entropy = compute_entropy(kmer, data->kmer_length);

// if( (count >= data->min_kmer_count) && (entropy >= data->min_any_entropy) )
    kv->add(key, keybytes, (char *) &count, sizeof(unsigned int));

}


// ##################################################################################################
//
void fileread(int itask, char *fname, KeyValue *kv, void *ptr)
{

  Data *data = (Data *) ptr;

  struct stat stbuf;
  int flag = stat(fname,&stbuf);
  if (flag < 0) {
    printf("ERROR: Could not query file size\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  string filename(fname);
  Fasta_reader fasta_reader(filename);

  while (true) {
        Fasta_entry fe = fasta_reader.getNext();
        string seq = fe.get_sequence();
        if (seq == "") break;

        unsigned int count = atoi(fe.get_header().c_str());
        kmer_int_type_t kmer_val = kmer_to_intval(seq);

        if(data->DS)
                kmer_val = get_DS_kmer_val(kmer_val, data->kmer_length);

        kv->add((char *) &kmer_val,sizeof(kmer_int_type_t), (char *) &count, sizeof(unsigned int));

  }
}

// #################################################################################################

void fileread_RNAseq(int itask, char *fname, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  struct stat stbuf;
  int flag = stat(fname,&stbuf);
  if (flag < 0) {
    printf("ERROR: Could not query file size\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  kmer_int_type_t max_kmer = get_maximum_kmer_intval(data->kmer_length);

  string filename(fname);
  Fasta_reader fasta_reader(filename);

  while (true) {
        Fasta_entry fe = fasta_reader.getNext();
        string seq = fe.get_sequence();
        if (seq == "") break;

        int nKmer = seq.length() - data->kmer_length + 1;
        for(int i=0; i<nKmer; i++) {
           string seq_kmer = seq.substr(i,data->kmer_length);
           if( contains_non_gatc(seq) ) {  
	      kmer_int_type_t kmer_val = kmer_to_intval(seq_kmer);
              if(data->DS) kmer_val = get_DS_kmer_val(kmer_val, data->kmer_length);

	      float entropy = compute_entropy(kmer_val, data->kmer_length);
              if( entropy > data->min_any_entropy ) {
		// int dummy = 0;
                // kv->add((char *) &kmer_val,sizeof(kmer_int_type_t),(char *) &dummy,sizeof(int));
                 kv->add((char *) &kmer_val,sizeof(kmer_int_type_t),(char *) &kmer_val,sizeof(kmer_int_type_t));
	      }
           }

        }
  }

}

// ####################################################################

void fileread_RNAseq_map_Edge(int itask, char *fname, KeyValue *kv, void *ptr)
{

  Data *data = (Data *) ptr;
  EDGE edge;

  struct stat stbuf;
  int flag = stat(fname,&stbuf);
  if (flag < 0) {
    printf("ERROR: Could not query file size\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  kmer_int_type_t max_kmer = get_maximum_kmer_intval(data->kmer_length);

  string filename(fname);
  Fasta_reader fasta_reader(filename);

  edge.ci = 0;
  edge.cj = 0;
  edge.ce = 0;

  while (true) {
	Fasta_entry fe = fasta_reader.getNext();
        string seq = fe.get_sequence();
        if (seq == "") break;

	int nKmer = seq.length() - data->kmer_length + 1;
        for(int i=0; i<(nKmer-1); i++) {
	   string seq_kmer_i = seq.substr(i,data->kmer_length);
	   string seq_kmer_j = seq.substr(i+1,data->kmer_length);

	   if( contains_non_gatc(seq_kmer_i) && contains_non_gatc(seq_kmer_j) ) {
	     kmer_int_type_t kmer_val_i = kmer_to_intval(seq_kmer_i);
	     kmer_int_type_t kmer_val_j = kmer_to_intval(seq_kmer_j);	

		if(data->DS) {
		   kmer_val_i =  get_DS_kmer_val(kmer_val_i, data->kmer_length); 
		   kmer_val_j =  get_DS_kmer_val(kmer_val_j, data->kmer_length);
		}		

		float entropy_i = compute_entropy(kmer_val_i, data->kmer_length);
                float entropy_j = compute_entropy(kmer_val_j, data->kmer_length);

		if( (entropy_i > data->min_any_entropy) && (entropy_j > data->min_any_entropy) ) {	
		  edge.vi = kmer_val_i;
		  edge.vj = kmer_val_j;
		  //int dummy = 0;
	   	  //kv->add((char *) &edge,sizeof(EDGE),(char *) &dummy,sizeof(int));

		 kv->add((char *) &edge,sizeof(EDGE),(char *) &edge,sizeof(EDGE));
		}
	   }
	}

  }

}

// #####################################################################

void map_output_kmers(uint64_t itask, char *key, int keybytes, char *value,
                  int valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;

  kmer_int_type_t kmer = *(kmer_int_type_t *) key;
  unsigned int count = *(unsigned int *) value;

  data->outFile << ">" << kmer << endl
                << decode_kmer_from_intval(kmer , data->kmer_length) << endl;
}

// ########################################################################

void reduce_kmers_RNAseq(char *key, int keybytes,
                        char *multivalue, int nvalues,
                        int *valuebytes, KeyValue *kv, void *ptr)
{
    Data *data = (Data *) ptr;
    unsigned int count = (unsigned int)nvalues;
    kv->add(key, keybytes, (char *) &count, sizeof(unsigned int));
    data->flag++;
}

// #########################################################################

void reduce_Edge_from_RNAseq(char *key, int keybytes,
                            char *multivalue, int nvalues,
                            int *valuebytes, KeyValue *kv, void *ptr)
{
  Data *data = (Data *) ptr;
  EDGE edge = *(EDGE *) key;
  edge.ce = (unsigned int)nvalues;
  kv->add((char *) &edge, sizeof(EDGE), NULL, NULL);
  data->flag++;
}

