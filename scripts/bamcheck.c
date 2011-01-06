/* 
    Author: petr.danecek@sanger
    gcc -Wall -g -O2 -I ~/cvs/samtools bamcheck.c -o bamcheck -lm -lz -L ~/cvs/samtools -lbam

    Assumptions and approximations:
        - GC content % calculation assumes that all reads have the same length (this can be fixed quite easily)
        - GC-depth does not split reads, the starting position determines which bin is incremented
*/

#define _ISOC99_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "sam.h"

#define GC_DEPTH 1

#define BWA_MIN_RDLEN 35
#define IS_PAIRED(bam) ((bam)->core.flag&BAM_FPAIRED && !((bam)->core.flag&BAM_FUNMAP) && !((bam)->core.flag&BAM_FMUNMAP))
#define IS_UNMAPPED(bam) ((bam)->core.flag&BAM_FUNMAP)
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE)

typedef struct
{
    float gc;
    uint32_t depth;
}
gc_depth_t;

typedef struct
{
    // Parameters
    int trim_qual;      // bwa trim quality

    // Dimensions of the quality histogram holder (quals_1st,quals_2nd), GC content holder (gc_1st,gc_2nd),
    //  insert size histogram holder
    int nquals;         // The number of quality bins 
    int nbases;         // The maximum sequence length the allocated array can hold
    int nisize;         // The maximum insert size that the allocated array can hold

    // Arrays for the histogram data
    uint64_t *quals_1st, *quals_2nd;
    uint64_t *gc_1st, *gc_2nd;
    uint64_t *isize;

    // The extremes encountered
    int max_len;            // Maximum read length
    int max_qual;           // Maximum quality
    float isize_main_bulk;  // There are always some unrealistically big insert sizes, report only the main part

    // Summary numbers
    uint64_t total_len;
    uint64_t nreads_1st;
    uint64_t nreads_2nd;
    uint64_t nreads_unmapped;
    uint64_t nreads_unpaired;
    uint64_t nreads_paired;
    uint64_t nbases_mapped;
    uint64_t nbases_mapped_cigar;
    uint64_t nbases_trimmed;  // bwa trimmed bases
    uint64_t nmismatches;

    // GC-depth related data
    uint32_t ngcd, igcd;        // The maximum number of GC depth bins and index of the current bin
    gc_depth_t *gcd;            // The GC-depth bins holder
    int gcd_bin_size;           // The size of GC-depth bin
    int32_t tid, pos;           // Position of the current bin

    // Auxiliary data
    double sum_qual;            // For calculation average quality value 
}
stats_t;

void error(const char *format, ...);

// Calculate the number of bases in the read trimmed by BWA
int bwa_trim_read(int trim_qual, uint8_t *quals, int len, int reverse) 
{
    if ( len<=BWA_MIN_RDLEN ) return 0;

    int max_trimmed = len-BWA_MIN_RDLEN;
    int l, sum=0, max_sum=0, max_l=0;

    for (l=0; l<max_trimmed; l++)
    {
        sum += trim_qual - quals[ reverse ? len-1-l : l];
        if ( sum<0 ) break;
        if ( sum>max_sum )
        {
            max_sum = sum;
            max_l   = l+1;
        }
    }
    return max_l;
}

void collect_stats(bam1_t *bam_line, stats_t *stats)
{
    int seq_len = bam_line->core.l_qseq;
    if ( !seq_len ) return;
    if ( seq_len >= stats->nbases )
        error("TODO: read length too long %d>=%d\n",seq_len,stats->nbases);
    if ( stats->max_len<seq_len )
        stats->max_len = seq_len;

    // Count GC
    uint8_t *seq  = bam1_seq(bam_line);
    int gc_count  = 0;
    int i;
    for (i=0; i<seq_len; i++)
    {
        // Conversion from uint8_t coding to ACGT
        //      -12-4---8-------
        //      =ACMGRSVTWYHKDBN
        if ( bam1_seqi(seq,i)==2 || bam1_seqi(seq,i)==4 ) gc_count++;
    }

    // Determine which array (1st or 2nd read) will these stats go to,
    //  trim low quality bases from end the same way BWA does, 
    //  fill GC histogram
    uint64_t *quals;
    uint8_t *bam_quals = bam1_qual(bam_line);
    if ( bam_line->core.flag&BAM_FREAD2 )
    {
        quals  = stats->quals_2nd;
        stats->nreads_2nd++;
        stats->gc_2nd[gc_count]++;
    }
    else
    {
        quals = stats->quals_1st;
        stats->nreads_1st++;
        stats->gc_1st[gc_count]++;
    }
    int reverse = IS_REVERSE(bam_line);
    if ( stats->trim_qual ) 
        stats->nbases_trimmed += bwa_trim_read(stats->trim_qual, bam_quals, seq_len, reverse);

    // Quality histogram and average quality
    for (i=0; i<seq_len; i++)
    {
        uint8_t qual = bam_quals[ reverse ? seq_len-i-1 : i];
        if ( qual>=stats->nquals )
            error("TODO: quality too high %d>=%d\n", quals[i],stats->nquals);
        if ( qual>stats->max_qual )
            stats->max_qual = qual;

        quals[ i*stats->nquals+qual ]++;
        stats->sum_qual += qual;
    }

    // Look at the flags and increment appropriate counters (mapped, paired, etc)
    if ( IS_UNMAPPED(bam_line) )
        stats->nreads_unmapped++;
    else
    {
        // The insert size is tricky, because for long inserts the libraries are
        // prepared differently and the pairs point in other direction. BWA does
        // not set the paired flag for them.  Similar thing is true also for 454
        // reads. Therefore, do the insert size stats for all unmapped reads.
        int32_t isize = bam_line->core.isize;
        if ( isize<0 ) isize = -isize;
        if ( IS_PAIRED(bam_line) && isize!=0 )
        {
            stats->nreads_paired++;
            if ( isize >= stats->nisize ) { isize=stats->nisize-1; }
            stats->isize[isize]++;
        }
        else
            stats->nreads_unpaired++;

        // Number of mismatches 
        uint8_t *nm = bam_aux_get(bam_line,"NM");
        if (nm) 
            stats->nmismatches += bam_aux2i(nm);

        // Number of mapped bases from cigar 
        if ( bam_line->core.n_cigar == 0) 
            error("FIXME: mapped read with no cigar?\n");
        for (i=0; i<bam_line->core.n_cigar; i++) 
        {
            // Conversion from uint32_t to MIDNSHP
            //  01-----
            //  MIDNSHP
            if ( (bam1_cigar(bam_line)[i]&BAM_CIGAR_MASK)==0 || (bam1_cigar(bam_line)[i]&BAM_CIGAR_MASK)==1 )
                stats->nbases_mapped_cigar += bam1_cigar(bam_line)[i]>>BAM_CIGAR_SHIFT;
        }

        stats->nbases_mapped += seq_len;

#if GC_DEPTH
        // GC-depth graph
        if ( stats->tid==-1 || stats->tid != bam_line->core.tid || bam_line->core.pos - stats->pos > stats->gcd_bin_size )
        {
            // First pass or new chromosome
            stats->tid = bam_line->core.tid;
            stats->pos = bam_line->core.pos;
            stats->igcd++;
        }
        if ( stats->igcd >= stats->ngcd )
            error("The genome too long?? [%ud]\n", stats->igcd);
        stats->gcd[ stats->igcd ].gc += (float) gc_count / seq_len;
        stats->gcd[ stats->igcd ].depth++;
#endif
    }

    stats->total_len += seq_len;
}

// Sort by GC and depth
#define GCD_t(x) ((gc_depth_t *)x)
static int gcd_cmp(const void *a, const void *b)
{
    if ( GCD_t(a)->gc < GCD_t(b)->gc ) return -1;
    if ( GCD_t(a)->gc > GCD_t(b)->gc ) return 1;
    if ( GCD_t(a)->depth < GCD_t(b)->depth ) return -1;
    if ( GCD_t(a)->depth > GCD_t(b)->depth ) return 1;
    return 0;
}

float gcd_percentile(gc_depth_t *gcd, int N, int p)
{
    float n,d;
    int k;

    n = p*(N+1)/100;
    k = n;
    if ( k<=0 ) 
        return gcd[0].depth;
    if ( k>=N ) 
        return gcd[N-1].depth;

    d = n - k;
    return gcd[k-1].depth + d*(gcd[k].depth - gcd[k-1].depth);
}

void output_stats(stats_t *stats)
{
    // Calculate average insert size and standard deviation. Save the number of reads with nonzero insert sizes
    //  to determine the main bulk
    int isize;
    uint64_t nisize = 0;
    double avg_isize=0, sd_isize=0;
    for (isize=1; isize<stats->nisize; isize++)
    {
        nisize += stats->isize[isize];
        avg_isize += isize * stats->isize[isize];
    }
    avg_isize /= nisize;
    for (isize=1; isize<stats->nisize; isize++)
        sd_isize += stats->isize[isize] * (isize-avg_isize)*(isize-avg_isize) / nisize;
    sd_isize = sqrt(sd_isize);

    printf("# This file was generated by bamcheck.\n");
    printf("# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.\n");
    printf("SN\tsequences:\t%ld\n", stats->nreads_1st+stats->nreads_2nd);
    printf("SN\t1st fragments:\t%ld\n", stats->nreads_1st);
    printf("SN\tlast fragments:\t%ld\n", stats->nreads_2nd);
    printf("SN\treads mapped:\t%ld\n", stats->nreads_paired+stats->nreads_unpaired);
    printf("SN\treads unmapped:\t%ld\n", stats->nreads_unmapped);
    printf("SN\treads unpaired:\t%ld\n", stats->nreads_unpaired);
    printf("SN\treads paired:\t%ld\n", stats->nreads_paired);
    printf("SN\ttotal length:\t%ld\n", stats->total_len);
    printf("SN\tbases mapped:\t%ld\n", stats->nbases_mapped);
    printf("SN\tbases mapped (cigar):\t%ld\n", stats->nbases_mapped_cigar);
    printf("SN\tbases trimmed:\t%ld\n", stats->nbases_trimmed);
    printf("SN\tmismatches:\t%ld\n", stats->nmismatches);
    printf("SN\terror rate:\t%e\n", (float)stats->nmismatches/stats->nbases_mapped_cigar);
    printf("SN\taverage length:\t%ld\n", (stats->nreads_1st+stats->nreads_2nd)?stats->total_len/(stats->nreads_1st+stats->nreads_2nd):0);
    printf("SN\tmaximum length:\t%d\n", stats->max_len);
    printf("SN\taverage quality:\t%.1f\n", stats->total_len?stats->sum_qual/stats->total_len:0);
    printf("SN\tinsert size average:\t%.1f\n", avg_isize);
    printf("SN\tinsert size standard deviation:\t%.1f\n", sd_isize);

    int ibase,iqual;
    if ( stats->max_len<stats->nbases ) stats->max_len++;
    if ( stats->max_qual+1<stats->nquals ) stats->max_qual++;
    printf("# First Fragment Qualitites. Use `grep ^FFQ | cut -f 2-` to extract this part.\n");
    printf("# Columns correspond to qualities and rows to cycles. First column is the cycle number.\n");
    for (ibase=0; ibase<stats->max_len; ibase++)
    {
        printf("FFQ\t%d",ibase+1);
        for (iqual=0; iqual<=stats->max_qual; iqual++)
        {
            printf("\t%ld", stats->quals_1st[ibase*stats->nquals+iqual]);
        }
        printf("\n");
    }
    printf("# Last Fragment Qualitites. Use `grep ^LFQ | cut -f 2-` to extract this part.\n");
    printf("# Columns correspond to qualities and rows to cycles. First column is the cycle number.\n");
    for (ibase=0; ibase<stats->max_len; ibase++)
    {
        printf("LFQ\t%d",ibase+1);
        for (iqual=0; iqual<=stats->max_qual; iqual++)
        {
            printf("\t%ld", stats->quals_2nd[ibase*stats->nquals+iqual]);
        }
        printf("\n");
    }
    printf("# GC Content of first fragments. Use `grep ^GCF | cut -f 2-` to extract this part.\n");
    for (ibase=0; ibase<stats->max_len; ibase++)
    {
        printf("GCF\t%.2f\t%ld\n", (float)ibase*100./stats->max_len,stats->gc_1st[ibase]);
    }
    printf("# GC Content of last fragments. Use `grep ^GCL | cut -f 2-` to extract this part.\n");
    for (ibase=0; ibase<stats->max_len; ibase++)
    {
        printf("GCL\t%.2f\t%ld\n", (float)ibase*100./stats->max_len,stats->gc_2nd[ibase]);
    }
    printf("# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part.\n");
    double bulk = 0;
    for (isize=1; isize<stats->nisize; isize++)
    {
        printf("IS\t%d\t%ld\n", isize,stats->isize[isize]);
        bulk += stats->isize[isize];
        if ( bulk/nisize > stats->isize_main_bulk ) 
            break;
    }

#if GC_DEPTH
    // Calculate average GC content, then sort by GC and depth
    printf("# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. The columns are: GC%%, unique sequence percentiles, 10th, 25th, 50th, 75th and 90th depth percentile\n");
    uint32_t igcd;
    for (igcd=0; igcd<stats->igcd; igcd++)
    {
        if ( stats->gcd[igcd].depth ) 
            stats->gcd[igcd].gc = round(100. * stats->gcd[igcd].gc / stats->gcd[igcd].depth);
    }
    qsort(stats->gcd, stats->igcd+1, sizeof(gc_depth_t), gcd_cmp);
    // for (igcd=0; igcd<stats->igcd; igcd++)
    //     printf("raw_GCD\t%d\t%f\n",stats->gcd[igcd].depth,stats->gcd[igcd].gc);
    igcd = 0;
    while ( igcd < stats->igcd )
    {
        // Calculate percentiles (10,25,50,75,90th) for the current GC content and print
        uint32_t nbins=0, itmp=igcd;
        float gc = stats->gcd[igcd].gc;
        while ( itmp<stats->igcd && fabs(stats->gcd[itmp].gc-gc)<0.1 )
        {
            nbins++;
            itmp++;
        }
        printf("GCD\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n", gc, (igcd+nbins+1)*100./(stats->igcd+1),
                gcd_percentile(&(stats->gcd[igcd]),nbins,10), 
                gcd_percentile(&(stats->gcd[igcd]),nbins,25), 
                gcd_percentile(&(stats->gcd[igcd]),nbins,50), 
                gcd_percentile(&(stats->gcd[igcd]),nbins,75), 
                gcd_percentile(&(stats->gcd[igcd]),nbins,90) 
              );
        igcd += nbins;
    }
#endif
}

void error(const char *format, ...)
{
    if ( !format )
    {
        printf("Usage: bamcheck [OPTIONS] file.bam\n\n");
        printf("Options:\n");
        printf("    -i, --insert-size <int>         Maximum insert size [8000]\n");
        printf("    -q, --trim-quality <int>        The BWA trimming parameter [0]\n");
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

    stats_t stats;
    stats.nquals = 95;
    stats.nbases = 200;
    stats.nisize = 8000;
    stats.max_len   = 30;
    stats.max_qual  = 40;
    stats.total_len = 0;
    stats.nreads_1st = 0;
    stats.nreads_2nd = 0;
    stats.nreads_unmapped = 0;
    stats.nreads_unpaired = 0;
    stats.nreads_paired   = 0;
    stats.nbases_mapped   = 0;
    stats.nbases_mapped_cigar = 0;
    stats.nmismatches     = 0;
    stats.sum_qual = 0;
    stats.isize_main_bulk = 0.99;
    stats.trim_qual = 0;
    stats.nbases_trimmed = 0;
    stats.gcd_bin_size = 20000;
    stats.ngcd         = 3e5;     // 300k of 20k bins is enough to hold a genome 6Gbp big
    stats.tid = stats.pos = -1;
    stats.igcd = 0;

    strcpy(in_mode, "rb");

    // Parse command line arguments
    for (i=1; i<argc; i++)
    {
        if ( !strcmp(argv[i],"-i") || !strcmp(argv[i],"--insert-size") )
        {
            if ( ++i>=argc )
                error(NULL);
            if ( sscanf(argv[i],"%d",&(stats.nisize)) != 1 )
                error("Expected integer after -i, got [%s]\n", argv[i]);
            continue;
        }
        if ( !strcmp(argv[i],"-q") || !strcmp(argv[i],"--trim-quality") )
        {
            if ( ++i>=argc )
                error(NULL);
            if ( sscanf(argv[i],"%d",&(stats.trim_qual)) != 1 )
                error("Expected integer after -q, got [%s]\n", argv[i]);
            continue;
        }
        if (!bam_fname) { bam_fname=argv[i]; continue; }
        error(NULL);
    }
    if ( !bam_fname )
    {
        if ( isatty(fileno((FILE *)stdin)) )
            error(NULL);
        bam_fname = "-";
    }


    // Init structures
    if ((sam = samopen(bam_fname, in_mode, NULL)) == 0) 
        error("Failed to open: %s\n", bam_fname);
    bam1_t *bam_line = bam_init1();
    stats.quals_1st = calloc(stats.nquals*stats.nbases,sizeof(uint64_t));
    stats.quals_2nd = calloc(stats.nquals*stats.nbases,sizeof(uint64_t));
    stats.gc_1st    = calloc(stats.nbases,sizeof(uint64_t));
    stats.gc_2nd    = calloc(stats.nbases,sizeof(uint64_t));
    stats.isize     = calloc(stats.nisize,sizeof(uint64_t));
    stats.gcd       = calloc(stats.ngcd,sizeof(gc_depth_t));

    // Collect statistics
    while (samread(sam,bam_line) >= 0) 
    {
        collect_stats(bam_line,&stats);
    }

    output_stats(&stats);

	return 0;
}



