###################################################################
#
# The default compiler is GNU mpic++.
# Run
#   make 
# to build MR-Inchworm and Fasta_Splitter 
# 

TARGETS=mr_inchworm fasta_splitter fasta_tool

all: ${TARGETS}
	sh install_tests.sh

mr_inchworm:
	cd MR_Inchworm && $(MAKE) -f Makefile.mpicc

fasta_splitter:
	cd Fasta_Splitter && $(MAKE)

fasta_tool:
	cd fastool && $(MAKE)

clean:
	cd MR_Inchworm && $(MAKE) -f Makefile.mpicc clean
	cd Fasta_Splitter && $(MAKE) clean 
	cd fastool && $(MAKE) clean
