CC=gcc
CFLAGS=-g -Wall -Wextra
INC=-I./include
LIB=-L./lib -lhts -lz 

all: vcf2msa htslib

.PHONY: all

vcf2msa: vcf2msa.c htslib
	$(CC) $(CFLAGS) -o $@ $< $(INC) $(LIB)

htslib: 
	cd htslib && $(MAKE) &&  $(MAKE) prefix=.. install

clean:
	rm -f vcf2msa
	rm -rf include bin lib share
