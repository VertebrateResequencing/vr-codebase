#!/usr/bin/env perl

BEGIN{
    my $ROOT = '/nfs/turing/work3/swat/';
    unshift(@INC, "$ROOT/ensembl/53/ensembl/modules");
    unshift(@INC, "$ROOT/ensembl/53/ensembl-variation/modules");
    unshift(@INC, "$ROOT/ensembl/53/ensembl-functgenomics/modules");
    unshift(@INC, "$ROOT/bioperl/bioperl-1.2.3");
}

use strict;
use warnings;
no warnings 'uninitialized';
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::TranscriptAlleles qw(get_all_ConsequenceType);
use Bio::EnsEMBL::Variation::AlleleFeature;
use Data::Dumper;
use Getopt::Long;

$|++;

my %known_spp = ('human' => 'Homo_sapiens',
                 'mouse' => 'Mus_musculus'
                 );

my ($snpfile, $species, $help);

GetOptions(
    'snps=s'	    =>  \$snpfile,
    'spp|species=s' =>  \$species,
    'h|help'        =>  \$help,
    );


my $spp = $known_spp{lc($species)};

(-f $snpfile && $spp && !$help) or die <<USAGE;
    Usage: $0   
                --snps   <file of snps>
                --spp    <species - human or mouse>
                --help <this message>

Calls consequences of snps against ensembl 53.
Outputs one line for each consequence of each snp.

SNP file should be sorted by chr, start and be tab-delimited with fields:
chromosome
position
reference base
snp base (single allele)
[any additional fields]
 
The fields in the SNP file (including any additional fields) will be included at the start of each line of the output.

Consequence fields output at the end of the original line are: type, geneid, transid, transtrand, aa_start_in_bp, aa_end_in_bp, aa_number_in_pep, ref_codon, snp_codon

USAGE

# snp types copied from Bio::EnsEMBL::Variation::ConsequenceType
our %CONSEQUENCE_TYPES = ('ESSENTIAL_SPLICE_SITE' => 1,
			  'STOP_GAINED' => 2,
			  'STOP_LOST' => 4,
			  'COMPLEX_INDEL' => 8,
			  'FRAMESHIFT_CODING' => 16,
			  'NON_SYNONYMOUS_CODING' => 32,
			  'SPLICE_SITE' => 64,
			  'SYNONYMOUS_CODING' => 128,
			  'REGULATORY_REGION' => 256,
			  'WITHIN_MATURE_miRNA' =>512,
			  '5PRIME_UTR' => 1024,
			  '3PRIME_UTR' => 2048,
			  'UTR'        => 2094,
			  'INTRONIC' => 4096,
			  'WITHIN_NON_CODING_GENE' => 8192,
			  'UPSTREAM' => 16384,
			  'DOWNSTREAM' => 32768,
			  'INTERGENIC' => 65536,
			  '_'          => 65537,
			  );


# header
#print join "\t", ('chr', 'pos', 'ref','snp', 'freq', 'type','gene','transcript','strand','codon_start','codon_end','codon_number','ref_codon','snp_codon');
#print "\n";
my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db (	-host => 'ensdb-archive',
                                -user => 'ensro',
                                -port => 5304,
				-verbose => 0,
                            );

my $ga = $reg->get_adaptor($spp, "core", "Gene");
my $sa = $reg->get_adaptor($spp, "core", "Slice");
my $ta = $reg->get_adaptor($spp,"core","Transcript");
my ($highest_cs) = @{$sa->db->get_CoordSystemAdaptor->fetch_all()};
my $assembly = $highest_cs->version();
warn "assembly: $assembly\n";

# read SNP file, and foreach snp, get slice of snp coord, get transcripts
# overlapping that slice, and foreach of those transcripts, calculate the
# consequences for the snp on that transcript

my ($curr_chr, $chrslice);
open (my $SNPS, "$snpfile") or die "Can't open snpfile : $!\n";
while (<$SNPS>){
    chomp;
    my $line = $_;
    my ($chr, $pos, $ref, $snp, @rest)=split "\t", $line;
    
    unless ($chr eq $curr_chr){
	$curr_chr = $chr;
	$chrslice = $sa->fetch_by_region('chromosome',$chr);
    }

    # skip hets for now
    next unless $snp =~ /^[ATGC]$/;
    #next unless $qual >= $qual_cutoff;

    my $af = Bio::EnsEMBL::Variation::AlleleFeature->new (
						    -start   => $pos,
						    -end     => $pos,
						    -strand  => 1,
						    -slice   => $chrslice,
						    -allele_string => $snp,
						    -variation_name => 'test',
							);
    my %consequences;
    my $snpslice = $sa->fetch_by_region('chromosome',$chr, $pos, $pos);
    die "Reference mismatch\n" unless $ref eq $snpslice->seq;
    my @transcripts = @{$snpslice->get_all_Transcripts()};

    if (@transcripts){
	foreach my $trans (@transcripts){
	    # only call consequences on protein_coding transcripts
	    next unless $trans->biotype eq 'protein_coding';
	    $trans = $trans->transfer($chrslice);
	    my $transtrand = $trans->strand;
	    my $transid = $trans->stable_id;
	    my $gene = $ga->fetch_by_transcript_stable_id($transid);	
	    my $geneid = $gene->external_name||$gene->stable_id;
	    my $cons = get_all_ConsequenceType($trans, [$af]);
	    foreach my $ct (@{$cons}){
		foreach my $type (@{$ct->type}){
		    warn "No such type $type\n" unless $CONSEQUENCE_TYPES{$type};
		    my $ref_codon = "";
		    my $snp_codon = "";
		    my $codon_start = "";
		    my $codon_end = "";
		    my $codon_number="";
		    if ($ct->aa_alleles){   # coding
			$codon_number = $ct->aa_start;
			die "SNP $codon_number in $transid spans codons\n" unless $ct->aa_start eq $ct->aa_end;
			my @codoncoords = $trans->pep2genomic($codon_number,$codon_number,$transtrand);
			my @sortcodons = sort {$a->start <=> $b->start} @codoncoords;
			my $codonseq;
			foreach (@codoncoords){
			    my $codonslice = $sa->fetch_by_region('chromosome',$chr, $_->start, $_->end, $transtrand);
			    $codonseq .= $codonslice->seq;
			}
			$codon_start = $sortcodons[0]->start;
			$codon_end = $sortcodons[-1]->end;

			$ref_codon = $codonseq;
			$snp_codon = $ct->codon;
		    }

		    print join "\t", ($line, $type, $geneid, $transid, $transtrand, $codon_start, $codon_end, $codon_number, $ref_codon, $snp_codon);
		    print "\n";
		}
	    }
	}
    }
    else {
	print join "\t", ($line, "INTERGENIC", "", "", "", "", "","","","");
	print "\n";
    }
}
