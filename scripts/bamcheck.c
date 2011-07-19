/* 
    Author: petr.danecek@sanger
    gcc -Wall -Winline -g -O2 -I ~/cvs/samtools bamcheck.c -o bamcheck -lm -lz -L ~/cvs/samtools -lbam

    Assumptions, approximations and other issues:
        - GC-depth does not split reads, the starting position determines which bin is incremented
        - coverage distribution ignores softclips and deletions
        - some stats require sorted BAMs
        - GC content graph can have saw-like pattern when BAM contains multiple read lengths. This is 
            unavoidable, consider for example uneven mixture of reads with lengths 3 and 4: 50% GC bin cannot 
            be accessed with the read length of 3 and 25% bin cannot be accessed with the read length of 4.
*/

#define _ISOC99_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "sam.h"
#include "faidx.h"
#include "khash.h"

#define BWA_MIN_RDLEN 35
#define IS_PAIRED(bam) ((bam)->core.flag&BAM_FPAIRED && !((bam)->core.flag&BAM_FUNMAP) && !((bam)->core.flag&BAM_FMUNMAP))
#define IS_UNMAPPED(bam) ((bam)->core.flag&BAM_FUNMAP)
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE)
#define IS_DUP(bam) ((bam)->core.flag&BAM_FDUP)

typedef struct 
{
    uint64_t len:32, line_len:16, line_blen:16;
    uint64_t offset;
} 
faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)

#ifndef _NO_RAZF
#include "razf.h"
#else
#ifdef _WIN32
#define ftello(fp) ftell(fp)
#define fseeko(fp, offset, whence) fseek(fp, offset, whence)
#else
extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);
#endif
#define RAZF FILE
#define razf_read(fp, buf, size) fread(buf, 1, size, fp)
#define razf_open(fn, mode) fopen(fn, mode)
#define razf_close(fp) fclose(fp)
#define razf_seek(fp, offset, whence) fseeko(fp, offset, whence)
#define razf_tell(fp) ftello(fp)
#endif

struct __faidx_t {
    RAZF *rz;
    int n, m;
    char **name;
    khash_t(s) *hash;
};


typedef struct
{
    float gc;
    uint32_t depth;
}
gc_depth_t;

typedef struct 
{
    int64_t pos;
    int size, start;
    int *buffer;
} 
round_buffer_t;

typedef struct
{
    // Parameters
    int trim_qual;      // bwa trim quality
    int rmdup;          // Exclude reads marked as duplicates from the stats

    // Dimensions of the quality histogram holder (quals_1st,quals_2nd), GC content holder (gc_1st,gc_2nd),
    //  insert size histogram holder
    int nquals;         // The number of quality bins 
    int nbases;         // The maximum sequence length the allocated array can hold
    int nisize;         // The maximum insert size that the allocated array can hold
    int ngc;            // The size of gc_1st and gc_2nd

    // Arrays for the histogram data
    uint64_t *quals_1st, *quals_2nd;
    uint64_t *gc_1st, *gc_2nd;
    uint64_t *isize;

    // The extremes encountered
    int max_len;            // Maximum read length
    int max_qual;           // Maximum quality
    float isize_main_bulk;  // There are always some unrealistically big insert sizes, report only the main part
    int is_sorted;

    // Summary numbers
    uint64_t total_len;
    uint64_t total_len_dup;
    uint64_t nreads_1st;
    uint64_t nreads_2nd;
    uint64_t nreads_dup;
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

    // Coverage distribution related data
    int ncov;                       // The number of coverage bins
    uint64_t *cov;                  // The coverage frequencies
    int cov_min,cov_max,cov_step;   // Minimum, maximum coverage and size of the coverage bins
    round_buffer_t cov_rbuf;        // Pileup round buffer

    // Auxiliary data
    double sum_qual;            // For calculation average quality value 
    samfile_t *sam;             // Unused
    faidx_t *fai;               // Reference sequence for GC-depth graph
}
stats_t;

void error(const char *format, ...);

// Coverage distribution methods
inline int coverage_idx(int min, int max, int n, int step, int depth)
{
    if ( depth < min )
        return 0;

    if ( depth > max )
        return n-1;

    return 1 + (depth - min) / step;
}

inline int round_buffer_lidx2ridx(int offset, int size, int64_t refpos, int64_t pos)
{
    return (offset + (pos-refpos) % size) % size;
}

void round_buffer_flush(stats_t *stats, int64_t pos)
{
    int ibuf,idp;

    if ( pos==stats->cov_rbuf.pos ) 
        return;

    if ( pos<0 || pos - stats->cov_rbuf.pos >= stats->cov_rbuf.size )
    {
        // flush the whole buffer
        for (ibuf=0; ibuf<stats->cov_rbuf.size; ibuf++)
        {
            if ( !stats->cov_rbuf.buffer[ibuf] ) 
                continue;

            idp = coverage_idx(stats->cov_min,stats->cov_max,stats->ncov,stats->cov_step,stats->cov_rbuf.buffer[ibuf]);
            stats->cov[idp]++;
            stats->cov_rbuf.buffer[ibuf] = 0;
        }
        stats->cov_rbuf.start = 0;
        stats->cov_rbuf.pos   = pos;
        return;
    }

    if ( pos < stats->cov_rbuf.pos ) 
        error("Expected coordinates in ascending order, got %ld after %ld\n", pos,stats->cov_rbuf.pos);

    int ifrom = stats->cov_rbuf.start;
    int ito = round_buffer_lidx2ridx(stats->cov_rbuf.start,stats->cov_rbuf.size,stats->cov_rbuf.pos,pos-1);
    if ( ifrom>ito )
    {
        for (ibuf=ifrom; ibuf<stats->cov_rbuf.size; ibuf++)
        {
            if ( !stats->cov_rbuf.buffer[ibuf] )
                continue;
            idp = coverage_idx(stats->cov_min,stats->cov_max,stats->ncov,stats->cov_step,stats->cov_rbuf.buffer[ibuf]);
            stats->cov[idp]++;
            stats->cov_rbuf.buffer[ibuf] = 0;
        }
        ifrom = 0;
    }
    for (ibuf=ifrom; ibuf<=ito; ibuf++)
    {
        if ( !stats->cov_rbuf.buffer[ibuf] )
            continue;
        idp = coverage_idx(stats->cov_min,stats->cov_max,stats->ncov,stats->cov_step,stats->cov_rbuf.buffer[ibuf]);
        stats->cov[idp]++;
        stats->cov_rbuf.buffer[ibuf] = 0;
    }
    stats->cov_rbuf.start = round_buffer_lidx2ridx(stats->cov_rbuf.start,stats->cov_rbuf.size,stats->cov_rbuf.pos,pos);
    stats->cov_rbuf.pos   = pos;
}

void round_buffer_insert_read(round_buffer_t *rbuf, int64_t from, int64_t to)
{
    if ( to-from >= rbuf->size )
        error("The read length too big (%d), please increase the buffer length (currently %d)\n", to-from+1,rbuf->size);
    if ( from < rbuf->pos )
        error("The reads are not sorted (%ld comes after %ld).\n", from,rbuf->pos);

    int ifrom,ito,ibuf;
    ifrom = round_buffer_lidx2ridx(rbuf->start,rbuf->size,rbuf->pos,from);
    ito   = round_buffer_lidx2ridx(rbuf->start,rbuf->size,rbuf->pos,to);
    if ( ifrom>ito )
    {
        for (ibuf=ifrom; ibuf<rbuf->size; ibuf++)
            rbuf->buffer[ibuf]++;
        ifrom = 0;
    }
    for (ibuf=ifrom; ibuf<=ito; ibuf++)
        rbuf->buffer[ibuf]++;
}

// Calculate the number of bases in the read trimmed by BWA
int bwa_trim_read(int trim_qual, uint8_t *quals, int len, int reverse) 
{
    if ( len<BWA_MIN_RDLEN ) return 0;

    // Although the name implies that the read cannot be trimmed to more than BWA_MIN_RDLEN,
    //  the calculation can in fact trim it to (BWA_MIN_RDLEN-1). (bwa_trim_read in bwa/bwaseqio.c).
    int max_trimmed = len - BWA_MIN_RDLEN + 1;
    int l, sum=0, max_sum=0, max_l=0;

    for (l=0; l<max_trimmed; l++)
    {
        sum += trim_qual - quals[ reverse ? l : len-1-l ];
        if ( sum<0 ) break;
        if ( sum>max_sum )
        {
            max_sum = sum;
            // This is the correct way, but bwa clips from some reason one base less
            // max_l   = l+1;
            max_l   = l;
        }
    }
    return max_l;
}

float fai_gc_content(faidx_t *fai,char *chr,int from, int to)
{
    khash_t(s) *h;
    khiter_t iter;
    faidx1_t val;
    uint32_t gc,count;
    char c;

    h = fai->hash;

    // ID of the sequence name
    iter = kh_get(s, h, chr);
    if (iter == kh_end(h)) 
        error("No such reference sequence [%s]?\n", chr);
    val = kh_value(h, iter);

    // Check the boundaries
    if (from > to ) 
        error("FIXME: %s:%d-%d\n", chr,from,to);
    if (from > 0) from--;
    if (from >= val.len)
        error("Was the bam file mapped with the reference sequence supplied? A read mapped beyond the chromosome (%s:%d, chromosome length %d).\n", chr,from+1,val.len);
    if (to >= val.len) to = val.len;

    // Position the razf reader
    razf_seek(fai->rz, val.offset + from / val.line_blen * val.line_len + from % val.line_blen, SEEK_SET);

    // Count GC content
    gc = count = 0;
    while (razf_read(fai->rz, &c, 1)==1 && count<to-from && !fai->rz->z_err)
    {
        if (isgraph(c))
        {
            if ( c=='G' || c=='g' || c=='C' || c=='C' )
            {
                gc++;
                count++;
            }
            else if ( c=='A' || c=='a' || c=='T' || c=='t' )
                count++;
        }
    }
    return (float)gc/count;
}


void realloc_buffers(stats_t *stats, int seq_len)
{
    int n = 2*(1 + seq_len - stats->nbases) + stats->nbases;

    stats->quals_1st = realloc(stats->quals_1st, n*stats->nquals*sizeof(uint64_t));
    if ( !stats->quals_1st )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,n*stats->nquals*sizeof(uint64_t));
    memset(stats->quals_1st + stats->nbases*stats->nquals, 0, (n-stats->nbases)*stats->nquals*sizeof(uint64_t));

    stats->quals_2nd = realloc(stats->quals_2nd, n*stats->nquals*sizeof(uint64_t));
    if ( !stats->quals_2nd )
        error("Could not realloc buffers, the sequence too long: %d (2x%ld)\n", seq_len,n*stats->nquals*sizeof(uint64_t));
    memset(stats->quals_1st + stats->nbases*stats->nquals, 0, (n-stats->nbases)*stats->nquals*sizeof(uint64_t));

    stats->nbases = n;

    // Realloc the coverage distribution buffer
    int *rbuffer = calloc(sizeof(int),seq_len*10);
    n = stats->cov_rbuf.size-stats->cov_rbuf.start;
    memcpy(rbuffer,stats->cov_rbuf.buffer+stats->cov_rbuf.start,n);
    if ( stats->cov_rbuf.start>1 )
        memcpy(rbuffer+n,stats->cov_rbuf.buffer,stats->cov_rbuf.start);
    stats->cov_rbuf.start = 0;
    free(stats->cov_rbuf.buffer);
    stats->cov_rbuf.buffer = rbuffer;
}

void collect_stats(bam1_t *bam_line, stats_t *stats)
{
    if ( stats->rmdup && IS_DUP(bam_line) )
        return;

    int seq_len = bam_line->core.l_qseq;
    if ( !seq_len ) return;
    if ( seq_len >= stats->nbases )
        realloc_buffers(stats,seq_len);
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
    int gc_idx = gc_count*(stats->ngc-1)/seq_len;

    // Determine which array (1st or 2nd read) will these stats go to,
    //  trim low quality bases from end the same way BWA does, 
    //  fill GC histogram
    uint64_t *quals;
    uint8_t *bam_quals = bam1_qual(bam_line);
    if ( bam_line->core.flag&BAM_FREAD2 )
    {
        quals  = stats->quals_2nd;
        stats->nreads_2nd++;
        stats->gc_2nd[gc_idx]++;
    }
    else
    {
        quals = stats->quals_1st;
        stats->nreads_1st++;
        stats->gc_1st[gc_idx]++;
    }
    int reverse = IS_REVERSE(bam_line);
    if ( stats->trim_qual>0 ) 
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
            //  01--4--
            //  MIDNSHP
            if ( (bam1_cigar(bam_line)[i]&BAM_CIGAR_MASK)==0 || (bam1_cigar(bam_line)[i]&BAM_CIGAR_MASK)==1 )
                stats->nbases_mapped_cigar += bam1_cigar(bam_line)[i]>>BAM_CIGAR_SHIFT;
        }

        stats->nbases_mapped += seq_len;

        if ( stats->tid==bam_line->core.tid && bam_line->core.pos<stats->pos )
            stats->is_sorted = 0;

        if ( stats->is_sorted )
        {
            // GC-depth graph
            if ( stats->tid==-1 || stats->tid != bam_line->core.tid || bam_line->core.pos - stats->pos > stats->gcd_bin_size )
            {
                // First pass or a new chromosome. Initialize the positions and get the reference GC content for this bin
                stats->tid = bam_line->core.tid;
                stats->pos = bam_line->core.pos;
                stats->igcd++;

                if ( stats->igcd >= stats->ngcd )
                    error("The genome too long?? [%ud]\n", stats->igcd);

                if ( stats->fai )
                    stats->gcd[ stats->igcd ].gc = fai_gc_content(stats->fai,stats->sam->header->target_name[stats->tid],stats->pos,stats->pos+stats->gcd_bin_size);

                round_buffer_flush(stats,-1);
            }
            stats->gcd[ stats->igcd ].depth++;
            // When no reference sequence is given, approximate the GC graph but determinig GC from each bin
            if ( !stats->fai )
                stats->gcd[ stats->igcd ].gc += (float) gc_count / seq_len;

            // Coverage distribution graph
            round_buffer_flush(stats,bam_line->core.pos);
            round_buffer_insert_read(&(stats->cov_rbuf),bam_line->core.pos,bam_line->core.pos+seq_len-1);
        }
    }

    stats->total_len += seq_len;
    if ( IS_DUP(bam_line) )
    {
        stats->total_len_dup += seq_len;
        stats->nreads_dup++;
    }
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
#undef GCD_t

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
    // Calculate average insert size and standard deviation (from the main bulk data only)
    int isize, ibulk=0;
    uint64_t nisize = 0;
    for (isize=1; isize<stats->nisize; isize++)
        nisize += stats->isize[isize];

    double bulk=0, avg_isize=0, sd_isize=0;
    for (isize=1; isize<stats->nisize; isize++)
    {
        bulk += stats->isize[isize];
        avg_isize += isize * stats->isize[isize];

        if ( bulk/nisize > stats->isize_main_bulk )
        {
            ibulk  = isize+1;
            nisize = bulk;
            break;
        }
    }
    avg_isize /= nisize ? nisize : 1;
    for (isize=1; isize<ibulk; isize++)
        sd_isize += stats->isize[isize] * (isize-avg_isize)*(isize-avg_isize) / nisize;
    sd_isize = sqrt(sd_isize);


    printf("# This file was generated by bamcheck.\n");
    printf("# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.\n");
    printf("SN\tsequences:\t%ld\n", stats->nreads_1st+stats->nreads_2nd);
    printf("SN\tis paired:\t%d\n", stats->nreads_1st&&stats->nreads_2nd ? 1 : 0);
    printf("SN\tis sorted:\t%d\n", stats->is_sorted ? 1 : 0);
    printf("SN\t1st fragments:\t%ld\n", stats->nreads_1st);
    printf("SN\tlast fragments:\t%ld\n", stats->nreads_2nd);
    printf("SN\treads mapped:\t%ld\n", stats->nreads_paired+stats->nreads_unpaired);
    printf("SN\treads unmapped:\t%ld\n", stats->nreads_unmapped);
    printf("SN\treads unpaired:\t%ld\n", stats->nreads_unpaired);
    printf("SN\treads paired:\t%ld\n", stats->nreads_paired);
    printf("SN\treads duplicated:\t%ld\n", stats->nreads_dup);
    printf("SN\ttotal length:\t%ld\n", stats->total_len);
    printf("SN\tbases mapped:\t%ld\n", stats->nbases_mapped);
    printf("SN\tbases mapped (cigar):\t%ld\n", stats->nbases_mapped_cigar);
    printf("SN\tbases trimmed:\t%ld\n", stats->nbases_trimmed);
    printf("SN\tbases duplicated:\t%ld\n", stats->total_len_dup);
    printf("SN\tmismatches:\t%ld\n", stats->nmismatches);
    printf("SN\terror rate:\t%e\n", (float)stats->nmismatches/stats->nbases_mapped_cigar);
    float avg_read_length = (stats->nreads_1st+stats->nreads_2nd)?stats->total_len/(stats->nreads_1st+stats->nreads_2nd):0;
    printf("SN\taverage length:\t%.0f\n", avg_read_length);
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
    for (ibase=0; ibase<stats->ngc; ibase++)
    {
        // Skip the zero-values. The discrete read length leaves unpleasent saw pattern in the graphs
        if ( !stats->gc_1st[ibase] ) continue;
        printf("GCF\t%.2f\t%ld\n", ibase*100./(stats->ngc-1),stats->gc_1st[ibase]);
    }
    printf("# GC Content of last fragments. Use `grep ^GCL | cut -f 2-` to extract this part.\n");
    for (ibase=0; ibase<stats->ngc; ibase++)
    {
        // Skip the zero-values. The discrete read length leaves unpleasent saw pattern in the graphs
        if ( !stats->gc_1st[ibase] ) continue;
        printf("GCL\t%.2f\t%ld\n", ibase*100./(stats->ngc-1),stats->gc_2nd[ibase]);
    }
    printf("# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part.\n");
    for (isize=1; isize<ibulk; isize++)
        printf("IS\t%d\t%ld\n", isize,stats->isize[isize]);


    printf("# Coverage distribution. Use `grep ^COV | cut -f 2-` to extract this part.\n");
    printf("COV\t[<%d]\t%d\t%ld\n",stats->cov_min,stats->cov_min-1,stats->cov[0]);
    int icov;
    for (icov=1; icov<stats->ncov-1; icov++)
        printf("COV\t[%d-%d]\t%d\t%ld\n",stats->cov_min + (icov-1)*stats->cov_step, stats->cov_min + icov*stats->cov_step-1,stats->cov_min + icov*stats->cov_step-1,stats->cov[icov]);
    printf("COV\t[%d<]\t%d\t%ld\n",stats->cov_min + (stats->ncov-2)*stats->cov_step-1,stats->cov_min + (stats->ncov-2)*stats->cov_step-1,stats->cov[stats->ncov-1]);


    // Calculate average GC content, then sort by GC and depth
    printf("# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. The columns are: GC%%, unique sequence percentiles, 10th, 25th, 50th, 75th and 90th depth percentile\n");
    uint32_t igcd;
    for (igcd=0; igcd<stats->igcd; igcd++)
    {
        if ( stats->fai )
            stats->gcd[igcd].gc = round(100. * stats->gcd[igcd].gc);
        else
            if ( stats->gcd[igcd].depth ) 
                stats->gcd[igcd].gc = round(100. * stats->gcd[igcd].gc / stats->gcd[igcd].depth);
    }
    qsort(stats->gcd, stats->igcd+1, sizeof(gc_depth_t), gcd_cmp);
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
        printf("GCD\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", gc, (igcd+nbins+1)*100./(stats->igcd+1),
                gcd_percentile(&(stats->gcd[igcd]),nbins,10) *avg_read_length/stats->gcd_bin_size,
                gcd_percentile(&(stats->gcd[igcd]),nbins,25) *avg_read_length/stats->gcd_bin_size, 
                gcd_percentile(&(stats->gcd[igcd]),nbins,50) *avg_read_length/stats->gcd_bin_size, 
                gcd_percentile(&(stats->gcd[igcd]),nbins,75) *avg_read_length/stats->gcd_bin_size, 
                gcd_percentile(&(stats->gcd[igcd]),nbins,90) *avg_read_length/stats->gcd_bin_size 
              );
        igcd += nbins;
    }
}

void error(const char *format, ...)
{
    if ( !format )
    {
        printf("Usage: bamcheck [OPTIONS] file.bam\n");
        printf("Options:\n");
        printf("    -c, --coverage <int>,<int>,<int>    Coverage distribution min,max,step [1,1000,1]\n");
        printf("    -d, --remove-dups                   Exlude from statistics reads marked as duplicates\n");
        printf("    -h, --help                          This help message\n");
        printf("    -i, --insert-size <int>             Maximum insert size [8000]\n");
        printf("    -m, --most-inserts <float>          Report only the main part of inserts [0.99]\n");
        printf("    -q, --trim-quality <int>            The BWA trimming parameter [0]\n");
        printf("    -r, --ref-seq <file>                Reference sequence (required for GC-depth calculation).\n");
        printf("    -s, --sam                           Input is SAM\n");
        printf("\n");
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
    stats.ngc    = 1000+1;
    stats.nquals = 95;
    stats.nbases = 300;
    stats.nisize = 8000;
    stats.max_len   = 30;
    stats.max_qual  = 40;
    stats.total_len = 0;
    stats.total_len_dup = 0;
    stats.nreads_1st = 0;
    stats.nreads_2nd = 0;
    stats.nreads_dup = 0;
    stats.nreads_unmapped = 0;
    stats.nreads_unpaired = 0;
    stats.nreads_paired   = 0;
    stats.nbases_mapped   = 0;
    stats.nbases_mapped_cigar = 0;
    stats.nmismatches     = 0;
    stats.sum_qual = 0;
    stats.isize_main_bulk = 0.99;   // There are always outliers at the far end
    stats.trim_qual = 0;
    stats.nbases_trimmed = 0;
    stats.gcd_bin_size = 20000;
    stats.ngcd         = 3e5;     // 300k of 20k bins is enough to hold a genome 6Gbp big
    stats.tid = stats.pos = -1;
    stats.igcd = 0;
    stats.fai  = NULL;
    stats.is_sorted = 1;
    stats.rmdup = 0;
    stats.cov_min  = 1;
    stats.cov_max  = 1000;
    stats.cov_step = 1;

    strcpy(in_mode, "rb");

    // Parse command line arguments
    for (i=1; i<argc; i++)
    {
        if ( !strcmp(argv[i],"-d") || !strcmp(argv[i],"--remove-dups") )
        {
            stats.rmdup = 1;
            continue;
        }
        if ( !strcmp(argv[i],"-s") || !strcmp(argv[i],"--sam") )
        {
            strcpy(in_mode, "r");
            continue;
        }
        if ( !strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") )
            error(NULL);
        if ( !strcmp(argv[i],"-r") || !strcmp(argv[i],"--ref-seq") )
        {
            if ( ++i>=argc )
                error(NULL);
            stats.fai = fai_load(argv[i]);
            if ( stats.fai==0 )
                error(NULL);
            continue;
        }
        if ( !strcmp(argv[i],"-c") || !strcmp(argv[i],"--coverage") )
        {
            if ( ++i>=argc )
                error(NULL);
            if ( sscanf(argv[i],"%d,%d,%d",&(stats.cov_min),&(stats.cov_max),&(stats.cov_step)) != 3 )
                error(NULL);
            continue;
        }
        if ( !strcmp(argv[i],"-i") || !strcmp(argv[i],"--insert-size") )
        {
            if ( ++i>=argc )
                error(NULL);
            if ( sscanf(argv[i],"%d",&(stats.nisize)) != 1 )
                error("Expected integer after -i, got [%s]\n", argv[i]);
            continue;
        }
        if ( !strcmp(argv[i],"-m") || !strcmp(argv[i],"--most-inserts") )
        {
            if ( ++i>=argc )
                error(NULL);
            if ( sscanf(argv[i],"%f",&(stats.isize_main_bulk)) != 1 )
                error("Expected float after -m, got [%s]\n", argv[i]);
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
    //  .. coverage bins and round buffer
    if ( stats.cov_step > stats.cov_max - stats.cov_min + 1 )
    {
        stats.cov_step = stats.cov_max-stats.cov_min;
        if ( stats.cov_step <= 0 )
            stats.cov_step = 1;
    }
    stats.ncov = 3 + (stats.cov_max-stats.cov_min) / stats.cov_step;
    stats.cov_max = stats.cov_min + ((stats.cov_max-stats.cov_min)/stats.cov_step +1)*stats.cov_step - 1;
    stats.cov = calloc(sizeof(uint64_t),stats.ncov);
    stats.cov_rbuf.size = stats.nbases;
    stats.cov_rbuf.buffer = calloc(sizeof(int32_t),stats.cov_rbuf.size);
    stats.cov_rbuf.pos = 0;
    stats.cov_rbuf.start = 0;
    // .. bam
    if ((sam = samopen(bam_fname, in_mode, NULL)) == 0) 
        error("Failed to open: %s\n", bam_fname);
    stats.sam = sam;
    bam1_t *bam_line = bam_init1();
    // .. arrays
    stats.quals_1st  = calloc(stats.nquals*stats.nbases,sizeof(uint64_t));
    stats.quals_2nd  = calloc(stats.nquals*stats.nbases,sizeof(uint64_t));
    stats.gc_1st     = calloc(stats.ngc,sizeof(uint64_t));
    stats.gc_2nd     = calloc(stats.ngc,sizeof(uint64_t));
    stats.isize      = calloc(stats.nisize,sizeof(uint64_t));
    stats.gcd        = calloc(stats.ngcd,sizeof(gc_depth_t));

    // Collect statistics
    while (samread(sam,bam_line) >= 0) 
    {
        collect_stats(bam_line,&stats);
    }
    round_buffer_flush(&stats,-1);

    output_stats(&stats);

	return 0;
}



