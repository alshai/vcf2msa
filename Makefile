CC=gcc
CFLAGS=-Wall -Wextra
INC=-I./include
LIB=-L./lib -l:libhts.a -llzma -lz -lpthread -lbz2

all: htslib vcf2msa 

.PHONY: htslib

vcf2msa: vcf2msa.cpp
	$(CC) $(CFLAGS) -o $@ $^ $(INC) $(LIB)

htslib: 
	cd htslib && $(MAKE) &&  $(MAKE) prefix=.. install

clean:
	rm -f vcf2msa
	rm -rf include bin lib share
