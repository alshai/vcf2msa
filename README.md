# vcf2msa

Conversion of a multi-sample VCF file to a multiple sequence alignment in FASTA
format using [htslib](https://github.com/samtools/htslib)

## Compilation:

```
make
```

## Usage

```
./vcf2msa <reference fasta> <vcf>
```

output: within the current directory, saves a fasta file per sample contained in the VCF file.

ie. if the VCF file contains samples `S1`, `S2`, and `S2`, `S1.fa`, `S2.fa` and `S3.fa` will be produced.
The concatenation of these fasta files is the multiple sequence alignment. 

## TODO: 

- check name of input sequence against the sequence name associated with each VCF record
- output the reference sequence as part of the MSA
- check that the input VCF is sorted
- output one fasta file instead of multiple fasta files?
- support nested variants
- support subsetting of samples (currently does the MSA for ALL samples)
- support or ignore complex variants, ie. `<CNV>`, `<DUP>`, etc.
