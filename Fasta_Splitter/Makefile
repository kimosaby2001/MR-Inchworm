
OBJS = Fasta_Splitter.o FileParser.o mutil.o StringUtil.o ErrorHandling.o
EXECUTABLE = Fasta_Splitter

CXX=mpic++
CXXFLAG=-Wall -O3 
COMPILE=$(CXX) $(CXXFLAG) -c

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS) 
	$(CXX) -o $(EXECUTABLE) $(OBJS)

%.o: %.cpp
	$(COMPILE) -o $@ $<

clean:
	-rm $(OBJS) $(EXECUTABLE)

