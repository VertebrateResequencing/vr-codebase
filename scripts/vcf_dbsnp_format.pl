#!/usr/bin/env perl
#  Reformats a vcf into dbsnp submission format (snp and indel only), adding required additional hdrs and VRT INFO tag

use strict;
use warnings;
use Carp;
use Vcf;
use Data::Dumper;

my $opts = parse_params();
convert_file($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg )
    {
        croak @msg;
    }
    die
        "About: Adds required additional hdrs and INFO to vcf for dbsnp submission\n",
        "Usage: cat in.vcf | vcf_dbsnp_format -h <handle> -r <ref> [-g <genotypes>] > out.vcf\n",
        "Options:\n",
        "   -h, --handle <string>         dbSNP submissionl handle name\n",
        "   -r, --reference <string>      The RefSeq Assembly accession.version\n",
        "   -g, --genotypes <string>      File containing source/target genotype names where these need to be modified\n",
        "   -h, -?, --help                This help message.\n",
        "\n";
}

sub parse_params
{
    while (my $arg=shift(@ARGV))
    {
        if ( $arg eq '-h' || $arg eq '--handle' ) { $$opts{handle}=shift(@ARGV); next; }
        if ( $arg eq '-r' || $arg eq '--reference' ) { $$opts{reference}=shift(@ARGV); next; }
        if ( $arg eq '-g' || $arg eq '--genotypes' ) { $$opts{genotypes}=shift(@ARGV); next; }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }

        if ( !$$opts{reference} || !$$opts{handle} || !$$opts{type} || $$opts{type} ne 'snp' || $$opts{type} ne 'indel') { error() }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    return $opts;
}

sub convert_file
{
    my ($opts) = @_;
    my $MAX_VAR_LEN=50;

    my %gtypes;
    if ($$opts{genotypes}) {
        open FILE, $$opts{genotypes} or die $!;
        while (<FILE>) {
            next if /^#/;
            chomp;
            my ($from,$to) = split;
            $gtypes{$from}=$to;
        }
        close(FILE);
    }

    my $vcf_in  = Vcf->new(fh=>\*STDIN);
    my $vcf_out = Vcf->new();

    $vcf_in->parse_header();

    die 'first vcf header not the fileformat' unless $$vcf_in{header_lines}[0]{key} eq 'fileformat';
    $vcf_out->add_header_line($$vcf_in{header_lines}[0]);

    my $handle=$$opts{handle};
    my $ref = $$opts{reference};
    my $dt = sprintf "%d%02d%02d", (localtime)[5] + 1900, (localtime)[4] + 1, (localtime)[3];
    my $batch="${handle}_${dt}";

    $vcf_out->add_header_line( {'key' => 'fileDate', 'value' => $dt} );
    $vcf_out->add_header_line( {'key' => 'handle', 'value' => $handle} );
    $vcf_out->add_header_line( {'key' => 'batch', 'value' => $batch} );
    $vcf_out->add_header_line( {'key' => 'reference', 'value' => $ref} );
    $vcf_out->add_header_line({key=>'INFO',ID=>'VRT',Number=>-1,Type=>'String',Description=>"Variation type,1 - SNV: single nucleotide variation,2 - DIV: deletion/insertion variation,3 - HETEROZYGOUS: variable, but undefined at nucleotide level,4 - STR: short tandem repeat (microsatellite) variation, 5 - NAMED: insertion/deletion variation of named repetitive element,6 - NO VARIATON: sequence scanned for variation, but none observed,7 - MIXED: cluster contains submissions from 2 or more allelic classes,8 - MNV: multiple nucleotide variation with alleles of common length greater than 1,9 - Exception"});

    # Header lines for the genotype FORMAT fields
    for (my $i=1; $i<@{$$vcf_in{header_lines}}; $i++) {
        my $hdr_key = $$vcf_in{header_lines}[$i]{key};
        if ($hdr_key eq 'contig' or $hdr_key eq 'FILTER') {
            $vcf_out->add_header_line($$vcf_in{header_lines}[$i]);
        }
        elsif ($hdr_key eq 'FORMAT')  {
            next if $$vcf_in{header_lines}[$i]{ID} eq 'FI'; 
            $vcf_out->add_header_line($$vcf_in{header_lines}[$i]);
        }
    }

    # add pop id header line
    my @gtypes;
    for (my $i=9; $i<@{$$vcf_in{columns}}; $i++) {
        my $gtype =  $$vcf_in{columns}[$i];
        if ($$opts{genotypes}) {
            my $gt = $gtypes{$gtype};
            $vcf_out->add_header_line({'key' => 'population_id', 'value' => $gt});
        }
        else {
            $vcf_out->add_header_line({'key' => 'population_id', 'value' => $gtype});
        }
        push @gtypes,$gtype;
    }
    print $vcf_out->format_header();

    if ($$opts{genotypes}) {
        my @gts;
        foreach my $gt (@gtypes) {
            push(@gts,$gtypes{$gt});
        }
        print join("\t",('#CHROM','POS','ID','REF','ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'),@gts), "\n";
    }
    else {
        print join("\t",('#CHROM','POS','ID','REF','ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'),@gtypes), "\n";
    }

    $vcf_out->add_columns(@{$$vcf_in{columns}});

    my ($rec_tot,$ref_only,$none_passed,$snp_out,$mixed_vars,$indel_out) = (0,0,0,0,0,0);

    while (my $rec=$vcf_in->next_data_hash()) {

        $rec_tot++;
        print $vcf_in->format_line($rec);

        # genotypes must have FORMAT tag FI=1 (pass filter)
        my $passed = 0;
        for (my $i=0; $i<@gtypes; $i++) {
            my $gtype = $gtypes[$i];
            if ($$rec{gtypes}{$gtype}{FI}) {
                $passed++;
            }
            else {
                delete $$rec{gtypes}{$gtype} unless $$rec{gtypes}{$gtype}{FI};
            }
        }
        if (!$passed) {
            $none_passed++;
            next;
        }

        $vcf_in->remove_format_field($rec,'FI');

        # Check ALT lengths and create GT translation hash if we need to delete any
        my (@alts,%gt_trans, @gt_deletion);
        my ($snps,$indels) = (0,0);
        my ($i,$gt_new) = (1,1);


        for my $alt (@{$$rec{ALT}}) {

            if (abs(length($alt) - length($$rec{REF})) > $MAX_VAR_LEN) {
                push (@gt_deletion, $i);
            }
            else {
                $gt_trans{$i} = $gt_new;
                $gt_new++;
                push (@alts,$alt);
            }
            $i++;

            # Check they are all either snps or dels
            my ($type,$len,$ht) = $vcf_in->event_type($rec,$alt);
            if  ($type eq 's') {
                $snps++;
            }
            elsif  ($type eq 'i') {
                $indels++;
            }
        }
        if ($snps > 0 && $indels > 0) {
            $mixed_vars++;
            next;
        }
print '@alts ', Dumper(@alts);
print '%gt_trans ', Dumper(%gt_trans);
print '@gt_deletion ', Dumper(@gt_deletion);

        delete $rec->{INFO};

        # Need to modify the ALT and each affected GT if we removed any ALTS > max
        if (@gt_deletion) {
            $rec->{ALT} = \@alts;

            foreach my $gtype (keys %{$rec->{'gtypes'}}) {
                my $gt = $rec->{gtypes}{$gtype}{GT};
                my ($gt1,$gt2);
                if ($gt =~ /\//) {
                    ($gt1,$gt2) = split (/\//,$gt);
                }
                else {
                    ($gt1,$gt2) = split (/\|/,$gt);
                }
                
                if ($gt_trans{$gt1} && $gt_trans{$gt2}) {  # translate the GT numbers
                    $gt = $gt_trans{$gt1} . substr($gt,1,1) . $gt_trans{$gt2};
                    $rec->{gtypes}{$gtype}{GT} = $gt;

                    # Remove scores from PL field for deleted genotypes
                    my $pl = $rec->{'gtypes'}{$gtype}{'PL'};

                    my @dels;   # produce a zerobased index delete list 
                    foreach my $i(@gt_deletion) {
                        push (@dels,$i-1);
                    }
                    $pl=delete_gt($pl,\@dels);
                    $rec->{gtypes}{$gtype}{PL} = $pl;
                }
                else { # get rid of GT if it contains a deleted ALT
                    delete $rec->{gtypes}{$gtype};
                }
            }
        }

        # Check that we do not have just ref/ref gtypes
        my $variant=0;
        foreach my $gtype (keys %{$rec->{'gtypes'}}) {
            my $gt = $rec->{gtypes}{$gtype}{GT};
            $variant++ if $gt && $gt ne '0/0' && $gt ne '0|0';
        }
        if ($variant == 0) {
            $ref_only++;
            next;
        }
        $rec->{ID}=sprintf("%s.%d",$batch,$rec_tot);

        if ($snps > 0) {
            $snp_out++;
            $rec->{INFO}{VRT} = 1;
            print $vcf_out->format_line($rec);
        }
        elsif ($indels > 0) {
            $indel_out++;
            $rec->{INFO}{VRT} = 2;
            print $vcf_out->format_line($rec);
        }
    }
    print STDERR "Input VCF:$rec_tot Rejected mixed variants:$mixed_vars No GTs passed:$none_passed Ref Only:$ref_only SNPs:$snp_out Indels:$indel_out\n";
}

sub delete_gt {
# delete entries from a GT matrix for a given GT index
# eg delete genotype C from AA,AB,BB,AC,BC,CC,AD,BD,CD,DD  ->  AA,AB,BB,AD,BD,DD

    my ($gts,$idxs) = @_;
    my @gt=split(/,/,$gts);

    # parse GT data in canonical order into a matrix so we can manipulate the values
    # Given the alleles A,B,C,D, and genotypes AA,AB,BB,AC,BC,CC,AD,BD,CD,DD generate:
    # AA
    # AB BB
    # AC BC CC
    # AD BD CD DD

    my @matrix;
    my $i=0;
    my $j=0;
    for (my $c=0;$c<@gt;$c++) {
        $matrix[$i][$j] = $gt[$c];

        $j++;
        if ($j > $i) {
            $j=0;
            $i++;
        }
    }

    # delete each GT in the idx list
    foreach my $idx(reverse sort { $a <=> $b } @$idxs) {
        splice(@matrix,$idx,1);

        for (my $i=$idx;$i<@matrix;$i++) {
            splice(@{$matrix[$i]},$idx,1);
        }
    }

    # parse GT matrix back into string
    my $gt = '';

    for (my $i=0;$i<@matrix;$i++) {
        $gt .= ',' if $gt or $gt eq '0';
        $gt .= join(',',@{$matrix[$i]});
    }
    return $gt;
}
