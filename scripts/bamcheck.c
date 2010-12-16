/* 
    Author: petr.danecek@sanger
    gcc -Wall -g -O2 -I ~/cvs/samtools bamcheck.c -o bamcheck -lm -lz -L ~/cvs/samtools -lbam
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "sam.h"

void error(const char *format, ...)
{
    if ( !format )
    {
        printf("Usage: bamcheck [OPTIONS] file.bam\n\n");
    }
    else
    {
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }
    exit(-1);
}

int main(int argc, char *argv[])
{
    char *bam_fname = NULL;
    samfile_t *sam = NULL;
    char in_mode[5];
    int i;

    strcpy(in_mode, "rb");

    // Parse command line arguments
    for (i=1; i<argc; i++)
    {
        if (!bam_fname) { bam_fname=argv[i]; continue; }
        error(NULL);
    }
    if ( !bam_fname )
        error(NULL);


    // Init structures
    if ((sam = samopen(argv[optind], in_mode, NULL)) == 0) 
        error("Failed to open: %s\n", bam_fname);
    size_t nquals = 95;
    size_t nbases = 200;
    size_t *quals_1st = calloc(nquals*nbases,sizeof(size_t));
    size_t *quals_2nd = calloc(nquals*nbases,sizeof(size_t));
    size_t *quals;
    bam1_t *bam_line = bam_init1();
    int max_len   = 30;
    int max_qual  = 40;
    size_t total_len = 0;
    size_t nreads_1st = 0;
    size_t nreads_2nd = 0;
    size_t nreads_unmapped = 0;
    double sum_qual = 0;

    // Collect statistics
    int ret,seq_len;
    while ((ret = samread(sam,bam_line)) >= 0) 
    {
        seq_len = bam_line->core.l_qseq;
        if ( !seq_len ) continue;
        if ( seq_len > nbases )
            error("TODO: read length too long %d>%d\n",seq_len,nbases);
        if ( max_len<seq_len )
            max_len = seq_len;

        uint8_t *bam_quals = bam1_qual(bam_line);

        int offset = 0;
        if ( bam_line->core.flag&BAM_FREAD2 )
        {
            // The 2nd pair is stored differently. The sequences may have different length
            //  and their indexes must be counted from the end. The last base of the read
            //  is in fact the first, highest quality base.
            offset = nbases-seq_len;
            quals  = quals_2nd;
            nreads_2nd++;
        }
        else
        {
            quals = quals_1st;
            nreads_1st++;
        }

        if ( bam_line->core.flag&BAM_FUNMAP )
            nreads_unmapped++;

        for (i=0; i<seq_len; i++)
        {
            uint8_t qual = bam_quals[i];
            if ( qual>=nquals )
                error("TODO: quality too high %d>=%d\n", quals[i],nquals);
            if ( qual>max_qual )
                max_qual = qual;

            quals[ (i+offset)*nquals+qual ]++;
            sum_qual += qual;
        }
        total_len += seq_len;
    }

    // Output statistics
    printf("# This file was produced by bamcheck. Columns correspond to qualities (0,1,2,etc.) and rows to cycles (bases).\n");
    printf("# First column is the cycle (base) index. Negative indexes correspond to the second pair of the paired-end sequencing.\n");
    printf("#\n");
    printf("# Summary stats\n");
    printf("#\tsequences: %ld\n", nreads_1st+nreads_2nd);
    printf("#\t1st fragments: %ld\n", nreads_1st);
    printf("#\tlast fragments: %ld\n", nreads_2nd);
    printf("#\tunmapped: %ld\n", nreads_unmapped);
    printf("#\ttotal length: %ld\n", total_len);
    printf("#\taverage length: %ld\n", (nreads_1st+nreads_2nd)?total_len/(nreads_1st+nreads_2nd):0);
    printf("#\tmaximum length: %d\n", max_len);
    printf("#\taverage quality: %.1f\n", total_len?sum_qual/total_len:0);
    printf("#\n");
    int ibase,iqual;
    if ( max_len<nbases ) max_len++;
    if ( max_qual+1<nquals ) max_qual++;
    for (ibase=max_len-1; ibase>=0; ibase--)
    {
        printf("%d",-(ibase+1));
        for (iqual=0; iqual<=max_qual; iqual++)
        {
            printf(" %ld", quals_2nd[(nbases-1-ibase)*nquals+iqual]);
        }
        printf("\n");
    }
    for (ibase=0; ibase<max_len; ibase++)
    {
        printf("%d",ibase+1);
        for (iqual=0; iqual<=max_qual; iqual++)
        {
            printf(" %ld", quals_1st[ibase*nquals+iqual]);
        }
        printf("\n");
    }

	return 0;
}



#if 0
    // Conversion from uint8_t coding to ACGT
    char *bam_nt16_2nd_table = "=ACMGRSVTWYHKDBN";
    uint8_t *seq  = bam1_seq(bam_line);
    for (i=0; i<core->l_qseq; i++) printf("%c",bam_nt16_2nd_table[bam1_seqi(seq,i)]);
#endif

