#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use Getopt::Std;
#use Set::IntSpan;
use Set::IntSpan::Fast;
use Cwd;

# Author: Kim Wong kw10@sanger.ac.uk

# Parse and filter VCF, reformat to tabular format
# filter by: position, biotype, VEP consequence, PASS
#  VCF must be annoated with SO terms, VEP VCF output with consequences
# separated by "|". Optionally, the VCF has DP4 from MuTect extended output
# which has been added to INFO as DP4T (for DP4 Tumour).

##http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences

my @accepted_biotypes = (
	'3prime_overlapping_ncrna',
	'antisense',
	'IG_C_gene',
	'IG_D_gene',
	'IG_J_gene',
	'IG_LV_gene',
	'IG_V_gene',
	'IG_V_pseudogene',
	'lincRNA',
	'miRNA',
	'misc_RNA',
	'Mt_rRNA',
	'Mt_tRNA',
	'nonsense_mediated_decay',
	'non_stop_decay',
	'polymorphic_pseudogene',
	'processed_pseudogene',
	'processed_transcript',
	'protein_coding',
	'pseudogene',
	'retained_intron',
	'rRNA',
	'sense_intronic',
	'sense_overlapping',
	'snoRNA',
	'snRNA',
	'transcribed_processed_pseudogene',
	'transcribed_unprocessed_pseudogene',
	'translated_processed_pseudogene',
	'translated_unprocessed_pseudogene',
	'TR_V_gene',
	'TR_V_pseudogene',
	'unitary_pseudogene',
	'unprocessed_pseudogene',
); 

my @all_cons = (
	'3_prime_UTR_variant',
	'5_prime_UTR_variant',
	'coding_sequence_variant',
	'downstream_gene_variant',
	'feature_elongation',
	'feature_truncation',
	'frameshift_variant',
	'incomplete_terminal_codon_variant',
	'inframe_deletion',
	'inframe_insertion',
	'initiator_codon_variant',
	'intergenic_variant',
	'intron_variant',
	'mature_miRNA_variant',
	'missense_variant',
	'NMD_transcript_variant',
	'nc_transcript_variant',
	'non_coding_transcript_exon_variant',
	'non_coding_transcript_variant',
	'regulatory_region_ablation',
	'regulatory_region_amplification',
	'regulatory_region_variant',
	'splice_acceptor_variant',
	'splice_donor_variant',
	'splice_region_variant',
	'stop_gained',
	'stop_lost',
	'stop_retained_variant',
	'synonymous_variant',
	'TF_binding_site_variant',
	'TFBS_ablation',
	'TFBS_amplification',
	'transcript_ablation',
	'transcript_amplification',
	'upstream_gene_variant',
	'start_lost', # new; replaces initiator codin variant
	'protein_altering_variant', # new
);

my %opts;
getopts('f:m:b:o:t:pd:r:', \%opts);

# required params:
if (!$opts{f} || !$opts{m} || !$opts{o} ) {
	my $help = <<END;

  This is a tool to parse and filter VCF, and reformat to tabular format.
  Filter by: position, biotype, VEP consequence, PASS.
  VCF must be annoated with Ensembl's Variant Effect Predictor, using
  SO terms, and with consequences separated by "|". Optionally,
  the VCF has DP4 from MuTect 'extended output' which has
  been added to INFO as DP4T (for DP4 Tumour), eg: from
  run_mutect_and_annotate.pl

  The default for '-f' option is the first 11 consequences in this table:
    http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences

  Output is written to \$OUTDIR/\$sample-vs-\$ref-MuTect.vcf/txt where \$sample and \$ref are
  tumour and normal sample names from the original VCF or names replaced using the -r option.

  Usage: zcat in.vcf.gz | reformatVCF.pl -f default|all|consequencelist 
	-m default|extended -o vcf|table|both [-b biotypeslist] [-t targetregions] [-p]

	Options:

	-f all|default|list; filter VEP consequences
	   'all' -keep all
	   list is a "|" separated list, eg: "missense|stop_gained|stop_lost"
	   'default' -list is "/transcript_ablation|stop_g|stop_l|missense|splice_d|splice_a|
	   initiator_c|transcript_amp|frameshift_v|inframe_/"
	   REQUIRED

	-m extended|default; add or exclude DP4T in table format
	   'extended' -Mutect was run with --extended_output option and 
	               there is DP4T in the INFO column
	   'default' -Do not add DP4T info or DPT4 not available in VCF
	   REQUIRED

	-o vcf|table|both; output file format(s)
	   'vcf' -output filtered VCF only
	   'table' -output table format only
	   'both' -output VCF and table
	   REQUIRED

	-b list; filter by ENSEMBL transcript biotypes
	   list is a "|" separated list, eg: "protein_coding|non_coding"
	   OPTIONAL

	-t file; BED file (CHR,START,END; start is zero-based)
	   with regions to keep SNV calls.  Can be text file or
	   bgzip.
	    OPTIONAL

	-d path; optional output directory, otherwise prints to working dir

	-p print out PASS calls only

	-r file containing sample names to replace names found in
	   BAM/VCF header. The VCF created by VEP uses the sample
	   name in the BAM header. To replace tumour and/or normal
	   names use this option. Input file format: A[tab]B,
	   where A is the name in the BAM, B is the name you want
	   to replace with)

END
	die $help;
}	
my $prefix = $opts{n};
# check consequences
my $which = lc $opts{f};
my @conseq = split /\|/, $which;
if ($which ne 'all' && $which ne 'default') {
	foreach my $type (@conseq) {
		my ($match) = grep { $_ eq $type}  @all_cons;
		if (!$match) {
			die "Invalid consequence type: $type\n";
		}
	}

}
my $consequences;

if ($which eq 'default') {
	$consequences = "transcript_ablation|stop_g|stop_l|missense|splice_d|splice_a|initiator_c|transcript_amp|frameshift_v|inframe_|start_lost";
}
elsif ($which ne 'all') {
	$consequences = $which;
}

# check biotypes:

my $biotype = $opts{b} if $opts{b};
my @biotype = split /\|/, $biotype if $biotype;
foreach my $type (@biotype) {
	my ($match) = grep { $_ eq $type}  @accepted_biotypes;
	if (!$match) {
		die "Invalid transcript biotype: $type\n";
	}
}

my $mutect = lc $opts{m};
if ($mutect ne 'default' && $mutect ne 'extended') {
	die "choose -m default or extended (prints out DP4 of tumour)\n";
}

my $outformat = lc $opts{o};
if ($outformat ne 'vcf' && $outformat ne 'table' && $outformat ne 'both') {
	die "Choose vcf or table with -o option\n";
}


# check targets file:
my @targets;
my $sets;	
if ($opts{t}) {
	print STDERR "loading target regsions in $opts{t}\n";
	if ($opts{t} =~ /\.gz$/) {
		open T, "zcat $opts{t} |" || die "Can't open file $opts{t}\n";
	}
	else {
		open T, "<$opts{t}" || die "Can't open file $opts{t}\n";
	}
	my %spans;
	while (<T>) {
		chomp;
		push @{$spans{$1}},"$2-$3" if /^(\S+)\s+(\S+)\s+(\S+)/;
	}
	close T;
	for (sort keys %spans) {
		$sets->{$_} = Set::IntSpan::Fast->new( @{$spans{$_}} );
		print STDERR "adding $_\n";
	}
#	while (<T>) {
#		chomp;
#		my ($chrom,$s,$e) = split "\t", $_;
#		$s++; # bed file is zero basd
#		if (!$sets->{$chrom}) {
#			$sets->{$chrom} = new Set::IntSpan "$s-$e"
#		}
#		else {
#			$sets->{$chrom}+="$s-$e";
#			print STDERR "adding $chrom $s-$e\n";
#		}
#	}
#	close T;
	print STDERR "Done loading target regsions in $opts{t}\n";
}
	
my $dir = $opts{d} ?  $opts{d} : cwd();
mkdir($dir) if $opts{d} && ! -d $opts{d};
if (! -d -r -w -x $dir ) {
	die "No such directory or no permissions to read/write $dir\n";
}
# make a hash with old and new sample names
my %names;
if ($opts{r}) {
	open R, "<$opts{r}" or die "Can't open file $opts{r}\n";
	while (<R>) {
		my ($bamname, $newname) = ($1,$2) if /^(\S+)\s+(\S+)/;
		die ("Problem reading file $opts{r}\n") if !$bamname || !$newname;
		$names{$bamname}=$newname;
	}
	close R; 
} 

my $reference;
my $vcfheader;
my $headerline;
# read in VCF lines:
my $linecount=0;
my @csqorder;
my %csqindex;

while (<>) {
	chomp;
	my $line = $_;
	if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+).+CSQ=([^;]+)\S*\s+\S+\s+(\S+)\s+(\S+)/) {
		my ($chr,$pos) = ($1,$2);
		my $pass = $6;
		my $print = "$1\t$2\t$3\t$4\t$5\t$6\t"; 
		my $csq = $7; 
		my $gt1 = $8;
		my $gt2 = $9;
		my $dp4t;
		if ($mutect eq 'extended' && $line=~/DP4T=(\d+,\d+,\d+,\d+)/) {
			$dp4t=$1;
		}
		$linecount++;
		if ($linecount==1) {
			$headerline=~s/QUAL\t//;
			$headerline=~s/FILTER\t/MUTECT_FILTER\t/;
			$headerline=~s/INFO\t/DP4T\t/ if $mutect eq 'extended' && $dp4t;
			$headerline=~s/INFO\t// if $mutect ne 'extended' || !$dp4t;
			$headerline=~s/FORMAT/TumourAltBases (Ref,Alt)\tNormalAltBases (Ref,Alt)/;
			my ($header,$n1,$n2) = ($1,$2,$3) if $headerline=~/(.+\s+)(\S+)\s+(\S+)$/;
			my $sample = $gt1=~/^0:/ ? $n2 : $n1;
			my $ref = $gt1=~/^0:/ ? $n1 : $n2;
			if ($opts{r}) {
				my $sample1 = $names{$sample} if $names{$ref};
				my $ref1 = $names{$ref} if $names{$ref};
				if ($sample1 ne $sample) {
					print STDERR "Replacing $sample with $sample1\n";
					$sample = $sample1;
				}
				if ($ref1 ne $ref) {
					print STDERR "Replacing $ref with $ref1\n";
					$ref = $ref1;
				}
			}
			# Table output table or both:
			if ($outformat ne 'vcf') {
				open F, ">$dir/$sample-vs-$ref-MuTect.txt" or die "Can't open $dir/$sample-vs-$ref-MuTect.txt" ;
				print F "##Reformatted MuTect calls with ENSEMBL VEP consequence annotations\n";
				print F "##Filtered for conseqeuences: transcript_ablation|stop_lost|stop_gain|missense|splice_donor|splice_acceptor|initiator_codon|start_lost|frameshift_variant|inframe_insertion|inframe_deletion\n" if $which eq 'limit';
				print F "##Excludes consequences affecting NMD and non-coding transcripts\n" if $which eq 'limit';
				print F "##No consequence filtering\n" if $which eq 'all';
				print F $reference."\n";
				print F "##DPT4 is the number of reads in the tumour sample with forward ref, reverse REF, forward ALT and reverse ALT bases\n" if $mutect eq 'extended';
				print F "##SIFT score: deleterious is <=0.05; tolerated is >0.05; most severe score is reported only\n";
				print F "##TumourAltBases is the fraction of reads with ALT bases (Ref counts, Alt counts) in the Tumour sample\n";
				print F "##NormalAltBases is the fraction of reads with ALT bases (Ref counts, Alt counts) in the Normal sample\n";
				print F "##Consequence information see:http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences\n";
				print F "##First 11 most severe consequnences in the table from the above link are retained below\n";
				print F "##A site is listed multiple times if there is more than one consequence type, AA change type or gene affected\n";
				print F $header;
				print F "$sample\t$ref\t";
				print STDERR "PICK $csqindex{PICK}\n";
				print F "ENS_GENEID\tGENE_SYMBOL\tCONSEQUENCE\tCDS_POS\tPROT_POS\tAA_CHANGE\tTRANSCRIPTS\tTRANS_BIOTYPE\tSIFT\n" if !$csqindex{PICK};
				print F "ENS_GENEID\tGENE_SYMBOL\tCONSEQUENCE\tCDS_POS\tPROT_POS\tAA_CHANGE\tPICK\tTRANSCRIPTS\tTRANS_BIOTYPE\tSIFT\n"  if $csqindex{PICK};
			}
			if ($outformat ne 'table') {
				open V, ">$dir/$sample-vs-$ref-MuTect.vcf" or die "Can't open $dir/$sample-vs-$ref-MuTect.vcf" ;
				print V $vcfheader;
			}
		}
		# conseq and biotypes filtering:
		next if $opts{p}  && $pass ne 'PASS';
		next if $which ne 'all' && $csq !~ /$consequences/; # a list of consequences
		next if $biotype &&  $csq !~ /$biotype/;
		# position filter
		if ($opts{t}) {
			my $set = $sets->{$chr};
			print STDERR "skipping position $chr $pos\n" if ! $set->contains($pos);
			next unless $set->contains($pos);
		}
		if ($outformat eq 'table' || $outformat eq 'both') {
			my $newline;
			$newline .= $dp4t."\t" if $dp4t;
			my ($t_frac,$n_frac,$t_info,$n_info);
			if ($gt1=~/:2$/ || $gt1=~/\/1:/) {
				# GT:AD:BQ:DP:FA:SS
				while($gt1=~/^([^:]+):([^:]+):[^:]+:[^:]+:([^:]+)/g) {
					($t_frac,$t_info) = ("$3 ($2)", $1);
				};
				while($gt2=~/^([^:]+):([^:]+):[^:]+:[^:]+:([^:]+)/g) {
					($n_frac,$n_info) = ("$3 ($2)", $1);
				};
				$newline .= "$t_frac\t$n_frac\t$t_info\t$n_info\n";
				
			}
			elsif ($gt2=~/:2$/ || $gt2=~/\/1:/) {
				# GT:AD:BQ:DP:FA:SS
				while($gt2=~/^([^:]+):([^:]+):[^:]+:[^:]+:([^:]+)/g) {
					($t_frac,$t_info) = ("$3 ($2)", $1);
				};
				while($gt1=~/^([^:]+):([^:]+):[^:]+:[^:]+:([^:]+)/g) {
					($n_frac,$n_info) = ("$3 ($2)", $1);
				};
				$newline .= "$t_frac\t$n_frac\t$t_info\t$n_info\n";
				
			}
			$newline=~s/0\/0/HOM_REF/g; 
			$newline=~s/\t0$/\tHOM_REF/g; 
			$newline=~s/0\/1/HET/g; 
			$newline=~s/1\/1/HOM_ALT/g;
			my @checkcontrols = split /\t/, $newline;
			$print .= $newline;
			$csq =~ s/\|/\t/g;
			my @c = split /,/, $csq; 
			chomp @c;
			# Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL|SYMBOL_SOURCE|BIOTYPE"
			# human, vep81;  Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL_SOURCE|HGNC_ID|SIFT"
			my @tmp;
			for (@c) {
				if ($which eq 'default') {
					if ($biotype) {
						push @tmp, $_ if $_=~/$consequences/ && $_ !~/NMD_|nc_trans|non_coding/ && $_ =~/$biotype/;
					}
					else {
						push @tmp, $_ if $_=~/$consequences/ && $_ !~/NMD_|nc_trans|non_coding/;
					}
				}
				elsif ($biotype && $consequences) {
					push @tmp, $_ if $_=~/$consequences/ && $_=~/$biotype/;
				}
				elsif ($consequences) {
					push @tmp, $_ if $_=~/$consequences/;
				}
				elsif ($biotype) {
					push @tmp, $_ if $_=~/$biotype/;
				}
				else {
					push @tmp, $_;
				}
			}
			@c=@tmp; #if @tmp;
			@c=&formatgenes(@c);
			for (@c) {
				chomp $print;
				print F "$print\t$_\n";
			}
		}
		 #print original line to VCF too:
		if ($outformat eq 'both' || $outformat eq 'vcf') {
			print V $line."\n";
		}
	}
	elsif ($line =~ /#CHR/ ) {
		$headerline = $line;
		if ($outformat ne 'table') {
			# replace sample names if requested
			if ($opts{r}) {
				chomp $line;
				my @c = split "\t", $line;
				my $addline = "##SampleNames changed: ";
				foreach my $i (9..$#c) {
					if ($names{$c[$i]}) {
						$addline .= " $c[$i]:$names{$c[$i]};"; 
						$c[$i] = $names{$c[$i]};
					}
				}
				$line = $addline ? "$addline\n".join("\t", @c) :  join ("\t", @c);
			}
			$vcfheader .= "$line\n";
		}
	}
	elsif ($line =~ /##reference/) {
		$reference = $_;
		if ($outformat ne 'table') {
			$vcfheader .= "$line\n";
		}
	}
	elsif ($line=~/^##INFO=<ID=CSQ/) {
		my $format = $1 if $line=~/Format:\s+(\S+)"/;
		@csqorder = split /\|/, $format;
		chomp @csqorder;
		# record the index position of each annotation
		foreach my $i (0..$#csqorder) {
			$csqindex{$csqorder[$i]}=$i;
		}
		if ($outformat ne 'table') {
			$vcfheader .= "$line\n";
		}
		
	}
	elsif ($line=~/^#/ && $outformat ne 'table') {
		$vcfheader .= "$line\n";
	}
}
close V;
close F;


sub formatgenes {
	my @list = @_;
	# check for same gene name and consequence; put transcripts with same annotation on same line
	my %info;
 	#Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL|SYMBOL_SOURCE|BIOTYPE"
	foreach my $line (@list) {
		my @c = split "\t", $line;
		my $inf;
		# Genes, conseq, cds pos, aa pos, aa change
		#for (13,4,6,7,8) {
		#for ('Gene','SYMBOL','Consequence','CDS_position','Protein_position','Amino_acids') {
		for ('SYMBOL','Consequence','CDS_position','Protein_position','Amino_acids') {
			my $index = $csqindex{$_};
			$inf.= $c[$index] ? $c[$index] : '-';
			$inf .= "\t" unless $_ eq 'Amino_acids';
		}
		$csqindex{Gene} = $csqindex{Gene}; 
		$csqindex{Gene} = '-' if !$csqindex{Gene}; 
		#$c[1] = '-' if !$c[1];
		# push trans and biotypes for each gene/cons combination 
		push @{$info{$c[$csqindex{Gene}]}{"$inf"}{trans}}, $c[$csqindex{Feature}] || '-';
		push @{$info{$c[$csqindex{Gene}]}{"$inf"}{biotype}}, $c[$csqindex{BIOTYPE}] || '-';
		$info{$c[4]}{"$inf"}{pick}="PICK" if $c[$csqindex{PICK}] && $c[$csqindex{PICK}]==1;
		# order sift scores by most severe first
		##print STDERR "LINE $line\n";
		##print STDERR "sift is $c[$csqindex{SIFT}]\n";
		if ($c[$csqindex{SIFT}] && $c[$csqindex{SIFT}] ne "" && $c[$csqindex{SIFT}] >= 0) {
			if ( !defined $info{$c[$csqindex{Gene}]}{"$inf"}{sift}[0]  || ( $c[$csqindex{SIFT}] < $info{$c[$csqindex{Gene}]}{"$inf"}{sift}[0] )) {
				unshift (@{$info{$c[$csqindex{Gene}]}{"$inf"}{sift}},$c[$csqindex{SIFT}]);
			#	print STDERR "unshift $c[$csqindex{SIFT}]\n";
			}
			else {
				#push (@{$info{$c[$csqindex{Gene}]}{"$inf"}{sift}},$c[$csqindex{SIFT}]);
			#	print STDERR "add $c[$csqindex{SIFT}]\n";
			}
		}
	}
	my @res;
	# format line with gene info and list of transcripts and biotypes, SIFT, domains:
	foreach my $gene (keys %info) {
		foreach my $inf (keys %{$info{$gene}}) {
			my $line = "$gene\t$inf\t";
			$line .= $info{$gene}{$inf}{pick} || "-";
			$line .= "\t";
			$line .= join (",", @{$info{$gene}{$inf}{trans}});
			$line .= "\t";
			$line .= join (",", @{$info{$gene}{$inf}{biotype}});
			$line .= "\t";
			$line .= defined $info{$gene}{$inf}{sift}[0] && $info{$gene}{$inf}{sift}[0] >= 0 ?  $info{$gene}{$inf}{sift}[0]: '-' ; # most severe only
			push @res, $line;			
		}
	}
	return @res;
}
