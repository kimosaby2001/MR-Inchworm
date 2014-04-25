CC     = gcc
CFLAGS = -O2 -std=c99 -Wno-deprecated

all:kseq.h fastool.c
	$(CC) $(CFLAGS) fastool.c -o fastool -lz 

clean:
	rm -f *.o fastool
