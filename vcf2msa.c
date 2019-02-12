#include <cstdio>
#include <zlib.h>
#include <vector>
#include <sys/stat.h>
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/khash.h"
KSEQ_INIT(gzFile, gzread);
KHASH_MAP_INIT_INT(seqdict, kseq_t*);
typedef khash_t(seqdict) seqdict_t;

/* 
 * assuming that vcf file is SORTED (do bcftools sort)
 */

int vcf2msa(kseq_t* seq, htsFile* vcfp, bcf_hdr_t* hdr, int nfps, FILE** fps);
void open_vcf(const char* fname, htsFile** vcfp, bcf_hdr_t** hdr);
void align_alleles(bcf1_t* rec, char* ref, char* alt);

void open_vcf(const char* fname, htsFile** vcfp, bcf_hdr_t** hdr) {
    *vcfp = hts_open(fname, "r");
    *hdr = bcf_hdr_read(*vcfp);
    if (vcfp == NULL || hdr == NULL) {
        fprintf(stderr, "failure to read vcf!\n");
        exit(1);
    }
}

void align_alleles(bcf1_t* rec, char* ref, char* alt) {
    char *longer, *shorter;
    strcpy(ref, rec->d.allele[0]);
    strcpy(alt, rec->d.allele[1]);
    longer = ref; 
    shorter = alt;
    if (strlen(ref) < strlen(alt)) {
        longer = alt;
        shorter = ref;
    }
    for (size_t i = 0; i < strlen(longer) - strlen(shorter); ++i) {
        strcat(shorter, "-");
    }
}


int vcf2msa(kseq_t* seq, htsFile* vcfp, bcf_hdr_t* hdr, int nfps, FILE** fps) {
    bcf1_t* rec = bcf_init();
    int npos = 0, ppos = 0;
    int nsmpl = bcf_hdr_nsamples(hdr);
    int32_t* gts = NULL;
    int ngt;
    int ngt_arr = 0;
    char ref[512];
    char alt[512];
    for (int s = 0; s < nsmpl; ++s) {
        fprintf(fps[s*2],     ">%s_1.%s\n", hdr->samples[s], seq->name.s);
        fprintf(fps[(s*2)+1], ">%s_2.%s\n", hdr->samples[s], seq->name.s);
    }
    int nvars = 0, nskipvars = 0;
    while (ppos < seq->seq.l && !bcf_read(vcfp, hdr, rec)) {
        bcf_unpack(rec, BCF_UN_ALL);
        if (rec->n_allele > 2) continue;
        align_alleles(rec, ref, alt);
        npos = rec->pos;
        ++nvars;
        if (npos < ppos) { // overlapping variant
            ++nskipvars;
            // fprintf(stderr, "skipping overlapping variant at %d->%d: ", ppos, npos);
            // fprintf(stderr, "pvar = %s, nvar = %s\n", prec->d.allele[0], rec->d.allele[0]); // TODO: remove this
            continue;
        }
        ngt = bcf_get_genotypes(hdr, rec, &gts, &ngt_arr);
        if (ngt != nsmpl * 2) {
            fprintf(stderr, "error - genotype missing for a sample! expected %d, got %d\n", nsmpl*2, ngt);
            exit(1);
        }
        // TODO: add newline every 60bp
        for (int s = 0; s < nsmpl*2; ++s) { 
            char* allele = bcf_gt_allele(gts[s]) ? alt : ref;
            fwrite(seq->seq.s + ppos, sizeof(char), npos-ppos, fps[s]);
            fwrite(allele, sizeof(char), strlen(allele), fps[s]);
        }
        ppos = npos + strlen(rec->d.allele[0]);
    }
    for (int s = 0; s < nsmpl*2; ++s) {
        fwrite(seq->seq.s + ppos, sizeof(char), seq->seq.l-ppos, fps[s]);
        fwrite("\n", sizeof(char), 1, fps[s]);
    }
    fprintf(stderr, "nvars: %d, skipped %d (%.3f%%)\n", nvars, nskipvars, (double) nskipvars / (double) nvars);
    return 0;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: ./vcf2msa <fasta> <vcf> <output prefix>\n");
        exit(1);
    }
    char* fa_fname = argv[1];
    char* vcf_fname = argv[2];

    gzFile fafp; 
    kseq_t *seq;
    fafp = gzopen(fa_fname, "r");
    seq = kseq_init(fafp);

    // peek at header to prepare file pointers
    htsFile* vcfp;
    bcf_hdr_t* hdr;
    open_vcf(vcf_fname, &vcfp, &hdr);
    int nsmpl = bcf_hdr_nsamples(hdr);
    // std::vector<FILE*> outfps(nsmpl*2);
    FILE** outfps = (FILE**) malloc(sizeof(FILE*) * nsmpl * 2);
    for (int i = 0; i < nsmpl; ++i) {
        char* fname1 = strdup(hdr->samples[i]);
        fname1 = strcat(fname1, "_1.fa");
        char* fname2 = strdup(hdr->samples[i]);
        fname2 = strcat(fname2, "_2.fa");
        outfps[i*2] = fopen(fname1, "w");
        outfps[(i*2)+1] = fopen(fname2, "w");
        free(fname1);
        free(fname2);
    }
    hts_close(vcfp);

    // go through each sequence and append haplotype to the files
    // TODO: don't reopen the vcf every time
    while (kseq_read(seq) >= 0) {
        open_vcf(vcf_fname, &vcfp, &hdr);
        vcf2msa(seq, vcfp, hdr, nsmpl*2, outfps);
    }
    kseq_destroy(seq);
    for (int i = 0; i < nsmpl*2; ++i) {
        fclose(outfps[i]);
    }
    free(outfps);
    return 0;
}
