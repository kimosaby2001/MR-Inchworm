###################################################################
#
# The default compiler is GNU mpig++.
# Run
#   make 
# to build MR-Inchworm and Fasta_Splitter 
# 

TARGETS=mr_inchworm fasta_splitter

all: ${TARGETS}

mr_inchworm:
	cd MR_Inchworm && $(MAKE) -f Makefile.mpicc

fasta_splitter:
	cd Fasta_Splitter && $(MAKE)

clean:
	cd MR_Inchworm && $(MAKE) -f Makefile.mpicc clean
	cd Fasta_Splitter && $(MAKE) clean 
