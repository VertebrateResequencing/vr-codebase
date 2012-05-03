=head1 NAME

VertRes::Utils::GTypeCheckMpileup - module for genotype checking using samtools mpileup

=head1 SYNOPSIS

use VertRes::Utils::GTypeCheckMpileup;
my $o = VertRes::Utils::GTypeCheckMpileup->new();
my %hits = $o->check_genotype($bam, $snp_vcf, $snp_ped, $ref_fa, $outfile);

=head1 DESCRIPTION

Simple genotype checker.
Calls genotypes at specific sites in the genome, given in a vcf file,
using samtools mpileup.  Finds the 'best' match to the samples in the ped file.

=cut

package VertRes::Utils::GTypeCheckMpileup;
use base qw(VertRes::Base);

use strict;
use warnings;
use Carp;
use LSF;
use File::Spec;
use Utils;
use VertRes::Wrapper::samtools;
use VertRes::Utils::FileSystem;

=head2 new

  Title   : new
  Usage   : my $obj = VertRes::Utils::GTypeCheckMpileup();
  Function: Create a new VertRes::Utils::GTypeCheckMpileup object.
  Returns : VertRes::Utils::GTypeCheckMpileup object
  Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    return $self;
}

=head2 check_genotype

  Title   : check_genotype
  Usage   : my %results = VertRes::Utils::GTypeCheckMpileup();
  Function: Create a new VertRes::Utils::GTypeCheckMpileup object.
  Returns : hash. sample_id         => sample id from bam file ...
                  sample_score      => match score of the sample
                  best_hits         => hash of best hits.
                                       key = sample name
                                       value = [#matches, #sites compared]
                  second_hits       => hash, same as best_hits, but hits
                                       with the second best score.
                  same_score_hits   => hash, same as best_hits, but
                                       hits with the same score as the sample
                  best_match        => 1, iff the right sample was one of the
                                       best hits
                  best_match_unique => 1, iff the right sample was the only
                                       sanple with the highest score
                  best_score        => score of best match(es)
                  second_score      => score of second best match(es)

  Args    : bam     .. bam file
            snp_vcf .. vcf file containing snps.  Must have IDs which
                       match the IDs in the ped file
            snp_ped .. ped file contianing known genotypes
            ref_fa  .. reference fasta file

=cut

sub check_genotype {
    my ($self, $bam, $snp_vcf, $snp_ped, $ref_fa) = @_;

    my %location2snp; # "chr\tpos" => rs_id
    my %snps_on_rev_strand;
    my @regions;
    my %sample_snps;  # rs_id => "XX", where XX are the base calls
    my $input_sample_id;

    my $vfs = VertRes::Utils::FileSystem->new();
    my $tempdir = $vfs->tempdir();

    # load snp location info from the SNP vcf file
    open my $fh, $snp_vcf or $self->throw("Error opening file $snp_vcf");

    while (<$fh>) {
        next if /^#/;
        chomp;
        my @a = split /\t/;
        $location2snp{"$a[0]\t$a[1]"} = $a[2];
        push @regions,  "$a[0]:$a[1]-$a[1]";
        my @info = split /;/, $a[7];
        my $hits = grep /^RV$/, @info;

        if ($hits > 0) {
            $snps_on_rev_strand{$a[2]} = 1;
        }
    }

    close $fh;

    # make a bam of just the SNP regions, for faster SNP calling
    my $sw = VertRes::Wrapper::samtools->new();
    my $outbam_prefix = File::Spec->catfile($tempdir, 'subset');
    $sw->view($bam, "$outbam_prefix.bam", regions => [@regions], (h=>1, b=>1));
    $sw->sort("$outbam_prefix.bam", "$outbam_prefix.subset.sort");
    $sw->index("$outbam_prefix.subset.sort.bam", "$outbam_prefix.subset.sort.bam.bai");

    # make file of SNP sites for samtools
    my $snp_sites = File::Spec->catfile($tempdir, 'snp_sites');
    open $fh, '>', $snp_sites or $self->throw("Error opening $snp_sites");
    foreach (keys %location2snp) {
        print $fh "$_\n";
    }
    close $fh;

    # get genotype calls from input BAM
    my $snp_calls = File::Spec->catfile($tempdir, 'snps.vcf.gz');
    my $query_vcf = File::Spec->catfile($tempdir, 'snps.vcf.query');
    my $query_bcf =  File::Spec->catfile($tempdir, 'snps.cf.query');
    my $cmd = "samtools mpileup -C50 -Iug -l $snp_sites -f $ref_fa $outbam_prefix.subset.sort.bam > $query_bcf &&  bcftools view -gcv $query_bcf | bgzip -c > $snp_calls && tabix -p vcf $snp_calls && " . q/ query-vcf -f '%CHROM\t%POS[\t%GT]\n' / . "$snp_calls > $query_vcf";
    Utils::CMD($cmd);

    # get the sample ID from the vcf file we just made
    $cmd = qq[zcat $snp_calls | awk '/^#CHROM/ {print \$NF; exit}'];
    ($input_sample_id) = qx[$cmd];
    chomp $input_sample_id;

    # get genotype calls from the vcf file
    open $fh, $query_vcf or $self->throw("Error opening $query_vcf");
    while (<$fh>){
        chomp;
        my @a = split /\t/;
        my $id;

        if ($location2snp{"$a[0]\t$a[1]"}){
            $id = $location2snp{"$a[0]\t$a[1]"};
        }
        else {
            $self->throw("SNP called here: $a[0] $a[1], but position not in input SNP list");
        }

        my @tmp = split /\//, $a[2];
        @tmp = sort @tmp;
        $sample_snps{$id} = join('', @tmp);
    }

    close $fh;

    # samtools only reports true SNPs, so the above misses out on where the
    # sample agrees with the reference.  Get these from the bcf file
    open $fh, "bcftools view -g $query_bcf |" or die "error $!";

    while (<$fh>) {
        my @a = split /\t/;
        next if ($a[0] =~ /^#/ || $a[4] ne "." || $a[5] < 100);
        my $id =  $location2snp{"$a[0]\t$a[1]"};
        $sample_snps{$id} = "$a[3]$a[3]";
    }

    close $fh;

    # compare called SNPs with the genotype.ped file
    my %rsid2ped_column;
    my %match_scores;  # score => [rsid1, rsid2, ...]
    my $input_sample_score;

    open $fh, $snp_ped or die "Error opening $snp_ped";

    while (<$fh>){
        chomp;
        my @fields = split /\t/;
        my $sample_name;

        if (/^#FamilyID/) {
            for my $i (6 .. $#fields){
                $rsid2ped_column{$fields[$i]} = $i;
            }
        }
        else {
            $sample_name = $fields[1];
            my $match_count = 0;
            my $no_call_count = 0;
            my $test_sites = 0;
            my $mismatch_count = 0;
            my @mismatch_ids = ();

            while (my ($rsid, $sample_call) = each (%sample_snps)) {
                if ($snps_on_rev_strand{$rsid}) {
                    $fields[$rsid2ped_column{$rsid}] =~ tr/ACGT/TGCA/;
                }

                my @tmp = split / /,  $fields[$rsid2ped_column{$rsid}];
                @tmp = sort @tmp;
                my $test_sample_call = join ('', @tmp);

                next if ($test_sample_call =~ /N/);

                if ($test_sample_call eq $sample_call) {
                    $match_count++;
                }
                else {
                    $mismatch_count++;
                    push @mismatch_ids, $rsid;
                }
                $test_sites++;

            }

            my $score = $test_sites == 0 ? 0 : $match_count / $test_sites;
            push @{$match_scores{$score}}, [$sample_name, $match_count, $test_sites];
            if ($sample_name eq $input_sample_id) {
                $input_sample_score = $score;
            }
        }        
    }

    close $fh;

    # report results
    my %results;
    $results{sample_id} = $input_sample_id;

    my ($best_score, $second_score) = reverse sort {$a <=> $b} keys %match_scores;

    $results{best_score} = $best_score;

    foreach my $a (@{$match_scores{$best_score}}) {
        $results{best_hits}{$a->[0]} = [$a->[1], $a->[2]];
    }


    if (defined $second_score) {
        $results{second_score} = $second_score;
        foreach my $a (@{$match_scores{$second_score}}) {
            $results{second_hits}{$a->[0]} = [$a->[1], $a->[2]];
        }
    }

    $results{sample_score} = $input_sample_score;

    foreach my $a (@{$match_scores{$input_sample_score}}) {
        $results{same_score_hits}{$a->[0]} = [$a->[1], $a->[2]];
    }

    $results{best_match} = 0;
    $results{best_match_unique} = 0;

    if ($results{best_hits}{$input_sample_id}) {
        $results{best_match} = 1;
        $results{best_match_unique} = 1 if (1 == scalar keys %{$results{best_hits}})
    }
    

    return %results;
}


1;

