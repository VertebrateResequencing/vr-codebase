#!/usr/bin/env perl

use strict;
use warnings;

use Bio::Tools::GFF;

print STDERR "Parse GFF Test\n";
print STDERR "==============\n\n";



#my $bam_file = "/lustre/scratch101/pathogen/pathpipe-test/cp7/seq-pipelines/Citrobacter/rodentium/TRACKING/421/airLcr/SLX/airLcr_121145/4407_6#3/22.pe.raw.sorted.bam";

#my $gff_file = '/lustre/scratch103/pathogen/pathpipe/refs/Citrobacter/rodentium_ICC168/Citrobacter_rodentium_ICC168_v1.gff';

#my $gff_file = '../Citrobacter_rodentium_ICC168_v1.gff';
#my $gff_file = '/lustre/scratch103/pathogen/pathpipe/refs/Schistosoma/mansoni/Schistosoma_mansoni_3.1.gff';
#my $gff_file = '/lustre/scratch103/pathogen/pathpipe/refs/Schistosoma/mansoni/Schistosoma_mansoni_v5.gff';


my $gff_file = $ARGV[0] ? $ARGV[0]:'Unknown';
my $bam_file = $ARGV[1] ? $ARGV[1]:'Unknown';



print STDERR 'GFF File    : ',$gff_file,"\n";
print STDERR 'BAM File    : ',$bam_file,"\n";


my %GFF_Gene; 
my %GFF_Exon;
my %GFF_RPKM;


# Parse GFF File
print STDERR "Parsing GFF File...";

#my $count = 0; # debug
my $gff_parser = Bio::Tools::GFF->new(-gff_version => 3, -file => $gff_file);
while( my $feature = $gff_parser->next_feature())
{
    last unless defined($feature); # No more features
    next unless $feature->primary_tag eq 'CDS'; # Only CDS features.
    #$count++; # Debug.

    my ($gene_id, $seq_id,$gene_strand,$gene_start,$gene_end,$exon_length);

    # Get Gene ID and Sequence
    ($gene_id, my @junk) = $feature->get_tag_values('ID');
    $gene_id =~ s/^"|"$//g;

    # Discontinous features will have several lines with the same gene ID.
    unless(exists $GFF_Gene{$gene_id})
    {
        # Values from GFF
        $seq_id = $feature->seq_id();
        $gene_strand = $feature->strand;
        die("$gene_id: Feature on both strands.\n") unless $gene_strand;
        $gene_start  = $feature->start;
        $gene_end    = $feature->end;
        $exon_length = $feature->end - $feature->start + 1;
    }
    else
    {
        # Update Start, End and Exon_length of Gene
        ($seq_id,$gene_strand,$gene_start,$gene_end,$exon_length) = @{$GFF_Gene{$gene_id}};
        die("$gene_id: DNA sequence id doesn't match.\n") unless ($seq_id eq $feature->seq_id());
        die("$gene_id: DNA template strand doesn't match.\n") unless ($gene_strand == $feature->strand);
        $gene_start = $feature->start < $gene_start ? $feature->start : $gene_start;
        $gene_end   = $feature->end   > $gene_end   ? $feature->end   : $gene_end;
        $exon_length += $feature->end - $feature->start + 1;
    }


    #print "$gene_id,$seq_id,$gene_strand,$gene_start,$gene_end,$exon_length\n"; # Debug

    # Update Gene info
    $GFF_Gene{$gene_id} = [$seq_id,$gene_strand,$gene_start,$gene_end,$exon_length];

    # Add Exon Start and End.
    push @{$GFF_Exon{$gene_id}},[$feature->start,$feature->end];

}
print STDERR "done.\n";

# If no Bam File then Dump GFF and exit.
unless($ARGV[1])
{
    foreach my $gene_id (sort idsort keys %GFF_Gene)
    {
        print $gene_id,"\t",join("\t",@{$GFF_Gene{$gene_id}}),"\n";
    }
    exit;
}

# Get BAM reads here... 

print STDERR "Finding mapped reads from BAM File...";
my $flagstat = "samtools flagstat ".$bam_file." | grep mapped |";
die("Bam file doesn't exist") unless( -e $bam_file );
open(my $flagstat_fh, $flagstat) || die ("Cant run flagstat\n");
my $total_mapped_reads = <$flagstat_fh>;
$total_mapped_reads =~ /^(\d+)\s/;
$total_mapped_reads = $1;
close $flagstat_fh;
print STDERR "done.\n";

print STDERR 'Total Reads : ',$total_mapped_reads," (from flagstat)\n";

# Parse BAM.
print STDERR "Parsing BAM File...";
foreach my $gene_id (sort idsort keys %GFF_Gene)
{
    # Gene data.
    my ($seq_id,$gene_strand,$gene_start,$gene_end,$genelength) = @{$GFF_Gene{$gene_id}};

    #print $gene_id,"\t",join("\t",@{$GFF_Gene{$gene_id}}),"\n"; # Debug

    # Get Bam Slice ($offset bp on either side of gene)
    my $offset = 100;
    my $window_start = $gene_start - $offset;
    my $window_end   = $gene_end   + $offset;
    $window_start = $window_start < 1 ? 1:$window_start;
    #$window_end   = $window_end > $reference_size ? $reference_size:$window_end;

    my $bamfile_slice = "samtools view ".$bam_file." ".$seq_id.":".$window_start."-".$window_end." |";
    open( my $bam_fh, $bamfile_slice ) || die; 
    my $mapped_read = 0;
    my $mapped_read_antisense = 0;
    while(my $line = <$bam_fh>)
    {
        # Parse sam line
        my($qname, $flag, $rname, $read_position, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual) = split(/\t/,$line);
        my $read_strand = $flag & 16 ? -1:1; # $flag bit set to 0 for forward 1 for reverse.

        # Calculate read length from cigar
        my $read_length = 0;
        $cigar =~ s/(\d+)[MIS=X]/$read_length+=$1/eg;

        foreach my $exon (@{$GFF_Exon{$gene_id}})
        {
            # Match/Antisense Match
            my($exon_start,$exon_end) = @{$exon};
            if($read_position < $exon_end && ($read_position + $read_length - 1) >= $exon_start)
            {
                if($read_strand == $gene_strand)
                {
                    $mapped_read++;
                    last;
                }
                else
                {
                    $mapped_read_antisense++;
                    last;
                }
            }
        }
    }
    close $bam_fh;

    my $rpkm           = $mapped_read           / ( ($genelength/1000) * ($total_mapped_reads/1000000) );
    my $rpkm_antisense = $mapped_read_antisense / ( ($genelength/1000) * ($total_mapped_reads/1000000) );
    
    $rpkm           = sprintf("%.1f",$rpkm);
    $rpkm_antisense = sprintf("%.1f",$rpkm_antisense);

    $GFF_RPKM{$gene_id} = [$mapped_read, $mapped_read_antisense,$rpkm,$rpkm_antisense];

    #last if $count > 99; # Debug
} 
print STDERR "done.\n";


# Output RPKM.
print STDERR "Outputting RPKM.\n\n";
foreach my $gene_id (sort idsort keys %GFF_Gene)
{
    # Dump everything
    print STDOUT $gene_id,"\t",join("\t",@{$GFF_Gene{$gene_id}}),,"\t", join("\t",@{$GFF_RPKM{$gene_id}}), "\n";
}


exit;

# Sort on ID.
# ID is usually Sequence.Int
sub idsort
{
    my @A = split(/\./,$a,2);
    my @B = split(/\./,$b,2);

    $A[0] cmp $B[0] || $A[1] <=> $B[1];
}