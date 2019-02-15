#include <stdio.h>
#include <zlib.h>
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/khash.h"
KSEQ_INIT(gzFile, gzread);
KHASH_MAP_INIT_INT(seqdict, kseq_t*);
typedef khash_t(seqdict) seqdict_t;

/* 
 * assuming that vcf file is SORTED (do bcftools sort)
 */

void vcf2msa(kseq_t* seq, htsFile* vcfp, bcf_hdr_t* hdr, FILE** fps);
void open_vcf(const char* fname, htsFile** vcfp, bcf_hdr_t** hdr);
void align_alleles(bcf1_t* rec, char* ref, char* alt);
void open_smpl_fps(bcf_hdr_t* hdr, int* nsmpls, FILE*** fps);

void open_vcf(const char* fname, htsFile** vcfp, bcf_hdr_t** hdr) {
    *vcfp = hts_open(fname, "r");
    *hdr = bcf_hdr_read(*vcfp);
    if (vcfp == NULL || hdr == NULL) {
        fprintf(stderr, "failure to read vcf!\n");
        exit(1);
    }
}

void close_vcf(htsFile* vcfp, bcf_hdr_t* hdr) {
    bcf_hdr_destroy(hdr);
    hts_close(vcfp);
}

void open_smpl_fps(bcf_hdr_t* hdr, int* nsmpl, FILE*** fps) {
    int n = bcf_hdr_nsamples(hdr);
    *fps = (FILE**) malloc(sizeof(FILE*) * n * 2);
    // TODO: maybe we want the filenames to be dynamically allocated?
    char fname1[512];
    char fname2[512];
    for (int i = 0; i < n; ++i) {
        strcpy(fname1, hdr->samples[i]);
        strcpy(fname2, hdr->samples[i]);
        strcat(fname1, "_1.fa");
        strcat(fname2, "_2.fa");
        (*fps)[i*2] = fopen(fname1, "w");
        (*fps)[(i*2)+1] = fopen(fname2, "w");
    }
    *nsmpl = n;
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
    size_t diff = strlen(longer) - strlen(shorter);
    for (size_t i = 0; i < diff; ++i) {
        strcat(shorter, "-");
    }
}

/* INPUT
 * seq : reference sequence
 * vcfp : pointer to VCF file
 * hdr : pointer to header of VCF
 * fps: pointers to output files - must be size of at least 2*nsmpls
 *      where nsmpls is number of samples in the vcf file
 *
 * OUTPUT/EFFECTS:
 * sequences will be written to each file in fps
 */
void vcf2msa(kseq_t* seq, htsFile* vcfp, bcf_hdr_t* hdr, FILE** fps) {
    bcf1_t* rec = bcf_init();
    size_t npos = 0, ppos = 0;
    int nsmpl = bcf_hdr_nsamples(hdr);
    int ngt;
    int ngt_arr = 0;
    int nvars = 0, nskipvars = 0;
    int seqid = bcf_hdr_name2id(hdr, seq->name.s);
    int32_t* gts = NULL;
    /* TODO: dynamically allocate ref/alt to support larger indels */
    char ref[512];
    char alt[512];
    for (int s = 0; s < nsmpl; ++s) {
        fprintf(fps[s*2],     ">%s_1.%s\n", hdr->samples[s], seq->name.s);
        fprintf(fps[(s*2)+1], ">%s_2.%s\n", hdr->samples[s], seq->name.s);
    }
    while (ppos < seq->seq.l && !bcf_read(vcfp, hdr, rec) && !bcf_unpack(rec, BCF_UN_ALL)) {
        if (seqid != rec->rid || rec->n_allele > 2) {
            continue;
        }
        npos = rec->pos;
        ++nvars;
        if (npos < ppos) { // overlapping variant
            ++nskipvars;
            continue;
        }
        align_alleles(rec, ref, alt);
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
    free(gts);
    bcf_destroy(rec);
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
    FILE** outfps;
    int nsmpl;
    open_smpl_fps(hdr, &nsmpl, &outfps);
    close_vcf(vcfp, hdr);

    // TODO: don't reopen the vcf every time?
    while (kseq_read(seq) >= 0) {
        open_vcf(vcf_fname, &vcfp, &hdr);
        vcf2msa(seq, vcfp, hdr, outfps);
        close_vcf(vcfp, hdr);
    }
    kseq_destroy(seq);
    for (int i = 0; i < nsmpl*2; ++i) {
        fclose(outfps[i]);
    }
    free(outfps);
    gzclose(fafp);
    return 0;
}
