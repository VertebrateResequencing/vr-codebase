#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use Utils;
use VertRes::Utils::FileSystem;
use Cwd;
use Data::Dumper;

my %options = (
    run_ge => '/software/varinf/bin/run_ge',
    glf_dir => '/nfs/users/nfs_m/mh12/bin/glftools/glfv3',
    config => '', 
    tempdir => ''
);

my $ops_ok = GetOptions(
    \%options,
    'run_ge=s',
    'config=s',
    'tempdir=s'
);


if ($#ARGV != 0 or !($ops_ok) or !-s $options{config}) {
    die qq/usage:$0 outfile_prefix
options:
  -config\tfile describing studies and plexes to use (see below
  -run_ge\trun_ge executable [$options{run_ge}]
  -glf_dir\tPath to glf directory [$options{glf_dir}]
  -tempdir\tspecify existing tempdir with run_ge output to re-do

Produces a bin and snp_sites file across all the studies specified.

The config file needs to be tab-delimited and have fields:
Study abbreviation (in variation informatics database), e.g. 10KCHD
Plex name, e.g. W30467
Sample name field (either 0 for family name or 1 for individual name)
Sample name prefix (to be prepended to the sample name field) e.g. UK10K_CHD

/;
}

my $outfile_prefix = $ARGV[0];
my $vfs = VertRes::Utils::FileSystem->new();
my $tempdir = $options{tempdir};
$tempdir ||= $vfs->tempdir();
unless (-d $tempdir){
    mkdir $tempdir;
}
print "tempdir:$tempdir\n";

open (my $CFG, $options{config}) or die "Can't open ".$options{config}."$!\n";

my %study_config;
while (<$CFG>){
    chomp;
    my ($study, $plex, $field, $prefix) = split "\t", $_;
    $study_config{$study} = {plex => $plex,
                             namefield => $field,
                             prefix => $prefix
                             };
}


# For each of the studies, do run_ge to get the data
# pool all the snps, then call genotypes over all the snps
foreach my $study_name (keys %study_config){
    # Do we need to run run_ge, or are we using existing output?
    my $glf_out = File::Spec->catfile($tempdir, "run_ge_$study_name");
    if (-s "$glf_out.map" && -s "$glf_out.ped"){
        # existing output - don't re-run
    }
    else {
        # run run_ge
        my $plex = $study_config{$study_name}{plex};
        my $cmd = "$options{run_ge} -o $glf_out -i -pS -b37 T$study_name:P$plex:20";
        Utils::CMD($cmd);
    }
}

# Now pool all the snps
my %snps;
foreach my $study_name (keys %study_config){
    # convert run_ge output to glftools text genotype output
    my $glf_out = File::Spec->catfile($tempdir, "run_ge_$study_name");
    run_ge_map2ordered_snps("$glf_out.map", \%snps);
}

# sort the snps by chr, pos
my @sort_snps = sort {$snps{$a}{chr} cmp $snps{$b}{chr} || $snps{$a}{pos} <=> $snps{$b}{pos}} keys %snps;

my $genotype_dir = File::Spec->catdir($tempdir, 'Genotypes');
if (-d $genotype_dir){
    $vfs->move($genotype_dir,"${genotype_dir}_orig");
}
mkdir $genotype_dir or die "error making directory $genotype_dir";


# genotypes over all snps
foreach my $study_name (keys %study_config){
    my $glf_out = File::Spec->catfile($tempdir, "run_ge_$study_name");
    run_ge_ped2genotype_files("$glf_out.ped", $genotype_dir, $study_config{$study_name}, \@sort_snps,\%snps);
}

# write file of snp positions
my $snp_sites = "$outfile_prefix.snp_sites";
open my $f, '>', $snp_sites or die "Error opening $snp_sites";
foreach my $snp (@sort_snps) {
    print $f $snps{$snp}{chr},"\t",$snps{$snp}{pos},"\n";
}
close $f;

# run hapmap2bin
my $hapmap2bin = File::Spec->catfile($options{glf_dir}, 'hapmap2bin');
my $cwd = getcwd();
chdir $genotype_dir or die "Error chdir $genotype_dir";
my $cmd = "$hapmap2bin -f genotypes.fofn > hapmap2bin.out";
Utils::CMD($cmd);
chdir $cwd;
$vfs->copy(File::Spec->catfile($genotype_dir, 'hapmap2bin.out'), "$outfile_prefix.bin");
#$vfs->copy("$glf_out.map", "$outfile_prefix.map");
#$vfs->copy("$glf_out.ped", "$outfile_prefix.ped");

###############################################################################
###############################################################################
###############################################################################


# takes .map file made by run_ge and returns an
# hash of SNP chr & position against ID.
sub run_ge_map2ordered_snps {
    my ($infile, $hashref) = @_;

    open my $fh, $infile or die "error opening file $infile";

    while (<$fh>) {
        next if (/^#/);
        chomp;
        my ($chr, $snp, $dist, $pos) = split;
        # bit of sanity checking
        if (exists ($hashref->{$snp})){
            unless ($hashref->{$snp}{chr} == $chr &&
                    $hashref->{$snp}{pos} == $pos){
                die "snp $snp is at $chr:$pos in $infile, but already seen at".$hashref->{$snp}{chr}.":".$hashref->{$snp}{pos}."\n";
            }
        }
        else {
            $hashref->{$snp}={  chr => $chr,
                                pos => $pos 
                            };
        }
    }

    close $fh;
   

}


# takes info from run_ge_map2ordered_snps, and
# ped file made by run_ge.  Makes a genotype file
# per sample in the ped file, and a file of filenames
# of the genotype files.
sub run_ge_ped2genotype_files {
    my ($infile, $outdir, $config, $snplist, $snpinfo) = @_;
    my @ids_line;
    my $first_snp_column = 6;
    my %outfiles;

    open my $ifh, $infile or die "error opening file $infile";

    while (<$ifh>) {
        chomp;

        if (/^#/) {
            @ids_line = split /\t/;
        }
        else {
            my @data = split /\t/;
            # make a snp id => genotype mapping
            my %snp2genotype;
            foreach my $i ($first_snp_column..$#data){
                my $genotype = $data[$i];
                $genotype =~ s/ //;
                $snp2genotype{$ids_line[$i]} = $genotype;
            }

            # some run_ge output needs to use the individual id (col 1) prefixed
            # by the outfile_prefix as the sample name, while some needs to use
            # the family id (col 0) instead.
            my $sample_name;
            if ($config->{namefield} == 0){
                $sample_name = $data[0];
                $sample_name =~ s/^IID//;   # trim off prefix - necesary?   
            }
            else {
                $sample_name = $data[1];
            }
            $sample_name = $config->{prefix}.$sample_name;
            if ($outfiles{$sample_name}){
                die "sample name $sample_name is in $infile and has been seen before\n";
            }
            # write the genotype file for this sample
            my $outfile = File::Spec->catfile($outdir, $sample_name);
            open my $ofh, '>', $outfile or die "Error opening $outfile";
            foreach my $snp (@$snplist) {
                print $ofh join "\t", ($snpinfo->{$snp}{chr}, $snpinfo->{$snp}{pos}, $snp2genotype{$snp});
                print $ofh "\n";
            }
            close $ofh;
            $outfiles{$sample_name}++;
        }
    }
    close $ifh;

    # write the file of genotype filenames
    my $outfile = File::Spec->catfile($outdir, 'genotypes.fofn');
    open my $ofh, '>>', $outfile or die "Error opening $outfile";
    print $ofh join("\n", keys %outfiles);
    print $ofh "\n";
    close $ofh;
}
