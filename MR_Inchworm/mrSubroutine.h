#include "mpi.h"

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/stat.h"

#define __STDC_LIMIT_MACROS
#include "mapreduce.h"
#include "keyvalue.h"
#include "blockmacros.h"

#include "fstream"
#include "iostream"
#include "sstream"

#include "vector"
#include "map"
#include "set"
#include "algorithm"
#include "time.h"
#include "math.h"

#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"


using namespace MAPREDUCE_NS;
using namespace std;



void fileread(int, char *, KeyValue *, void *);

void fileread_RNAseq(int itask, char *fname, KeyValue *kv, void *ptr);

void fileread_RNAseq_map_Edge(int itask, char *fname, KeyValue *kv, void *ptr);

void reduce_kmers_RNAseq(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_Edge_from_RNAseq(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void map_output_kmers(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);

void reduce_merge_kmers(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void map_kmers_edge_vertex(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);

void reduce_filter_nonExisting_edge(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void map_kmer_edge (uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);

void reduce_filter_edge_by_count(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_keep_dominantEdge(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_count_to_edge(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_filter_edge(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void map_filter_edge_by_count(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);

void output_edge(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);


void edge_to_vertices(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);

void reduce_self_zone(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);


void map_edge_vert(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);

void reduce_edge_zone(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_zone_winner(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);


void map_invert_multi(uint64_t itask, char *key, int keybytes ,char *value, int valuebytes, KeyValue *kv, void *ptr);

void map_zone_multi(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);

void reduce_zone_reassign(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void map_strip(uint64_t itask, char *key, int keybytes, char *value, int valuebytes, KeyValue *kv, void *ptr);

void reduce_zone_kmer_count(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_print_clustered_kmers(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

void reduce_run_inchworm(char *key, int keybytes, char *multivalue, int nvalues, int *valuebytes, KeyValue *kv, void *ptr);

