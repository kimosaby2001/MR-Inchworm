
#define __STDC_LIMIT_MACROS
#include "stdint.h"
#include "keyvalue.h"
#include "sequenceUtil.hpp"

using MAPREDUCE_NS::KeyValue;

#define ALLBITS UINT64_MAX
#define INT64MAX INT64_MAX
#define HIBIT UINT64_MAX-INT64_MAX


typedef struct {
  kmer_int_type_t kmer;
  unsigned int count;
} KmerCount;

typedef kmer_int_type_t VERTEX;

typedef struct {
  VERTEX vi,vj;
  unsigned int ci, cj;
  unsigned int ce;
  uint64_t zi, zj; 
} EDGE;

typedef struct {
    uint64_t zone,empty;
  } PAD;

struct Data {
        int n;
	uint64_t flag;
        unsigned int kmer_length;
        ofstream outFile;
        int pshift;
        uint64_t lmask;
        int seed,nthresh;
        int me,nprocs;
        PAD pad;

        bool forward_direction;  //  forward = true, reverse = false
        float min_ratio_non_error;
        unsigned int min_kmer_count;
        unsigned int min_edge_count;

        bool DS;
	bool PACMAN;
	bool CRAWL;
	bool prune_error_kmers;

	unsigned int MAX_RECURSION;
	float MIN_SEED_ENTROPY;
	unsigned int MIN_SEED_COVERAGE;
	float min_any_entropy;
	unsigned int crawl_length;
	
	float MIN_CONNECTIVITY_RATIO;
	unsigned int MIN_ASSEMBLY_LENGTH;
	unsigned int MIN_ASSEMBLY_COVERAGE;
	bool WRITE_COVERAGE;
	string COVERAGE_OUTPUT_FILENAME;	
};



