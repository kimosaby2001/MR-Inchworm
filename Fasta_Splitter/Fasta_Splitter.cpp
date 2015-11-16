#include <mpi.h>

#include <stdio.h>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <map>
#include <fstream>

#include "CommandLineParser.h"
#include "FileParser.h"
#include "mutil.h"

#ifndef BIG_FILE_DEFINES_H
#define BIG_FILE_DEFINES_H

#ifdef __linux
     #ifndef _LARGEFILE64_SOURCE
     #define _LARGEFILE64_SOURCE
     #endif
     #define _FILE_OFFSET_BITS 64
     #define STAT_NAME stat64
#else
     #define STAT_NAME stat
#endif

#endif //BIG_FILE_DEFINES_H

void Execute_command(const char * command) {
    int ret = system(command);
    if (ret != 0) {
        cout << "COMMAND: " << command << endl;
        cout << "Died with exit code " << ret << endl;
        cout << "Exiting." << endl;
        exit(-1);
    }
}

string exec(const char * cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    string result = "";
    while(!feof(pipe)) {
    	if(fgets(buffer, 128, pipe) != NULL)
    		result += buffer;
    }
    pclose(pipe);
    return result;
}

static unsigned long long MAX_CLUSTER_SIZE = 201;

using namespace std;

int main (int argc, char* argv[]) {

   int rank,numranks;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numranks);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   stringstream cmdstr, cmd_out;

   commandArg<string> fastaCmmd("-r"    ,"read fasta");
   commandArg<string> headerCmmd("-i"   ,"header file");
   commandArg<string> outputCmmd("-o"   ,"Dir for splited fasta");

   commandLineParser P(argc,argv);
   P.registerArg(fastaCmmd);
   P.registerArg(headerCmmd);
   P.registerArg(outputCmmd);
 
   P.parse();

   string readFile = P.GetStringValueFor(fastaCmmd);
   string headerFile = P.GetStringValueFor(headerCmmd);
   string sKmerDir = P.GetStringValueFor(outputCmmd);

   cmdstr.str("");
   cmdstr << "wc -l " << headerFile;
//   cerr << cmdstr.str() << endl;
   string CMD_out = exec(cmdstr.str().c_str());
 
   cmd_out.str(""); 
   cmd_out << CMD_out;
   string dummy_str;
   cmd_out >> dummy_str;
   unsigned long long numLine = strtoull(dummy_str.c_str(), NULL, 10);

   unsigned long long nlines = numLine / numranks;
   unsigned long long remain = numLine % (unsigned long long)numranks;
 
   int START[MAX_CLUSTER_SIZE], END[MAX_CLUSTER_SIZE];

   for(unsigned long long i=0; i<(unsigned long long)numranks; i++) {
	unsigned long long interval = nlines;
	if(i<remain) { interval++; }	

	if(i == 0) START[i] = 1;
        else       START[i] = END[i-1] + 1;
        END[i] = START[i] + interval - 1;
   }

   cmdstr.str("");
   cmdstr << "sed '" << START[rank] << "q;d' " << headerFile;
   CMD_out = exec(cmdstr.str().c_str());
   stringstream cmd_start(CMD_out);
   string dummy_start;
   cmd_start >> dummy_start;
   unsigned long long startLine = strtoull(dummy_start.c_str(), NULL, 10);  

   unsigned long long endLine;
   if(rank < numranks-1) {
	unsigned long long dummy = END[rank]+1;
	cmdstr.str("");
   	cmdstr << "sed '" << dummy << "q;d' " << headerFile;
	
	CMD_out = exec(cmdstr.str().c_str());
        stringstream cmd_end(CMD_out);
        string dummy_end;
        cmd_end >> dummy_end;
        endLine = strtoull(dummy_end.c_str(), NULL, 10);
	endLine = endLine - 1;

  } else {
        cmdstr.str("");
        cmdstr << "wc -l " << readFile;

        CMD_out = exec(cmdstr.str().c_str());
        stringstream cmd_end(CMD_out);
        string dummy_end;
        cmd_end >> dummy_end;
        endLine = strtoull(dummy_end.c_str(), NULL, 10);
  }

   cmdstr.str("");
   cmdstr << "sed '" << startLine << "," << endLine << "p;d' " << readFile << " > " << sKmerDir << "/sKmer_" << rank << ".fa";
   Execute_command(cmdstr.str().c_str()); 
   
   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Finalize();
   return(0);

}


