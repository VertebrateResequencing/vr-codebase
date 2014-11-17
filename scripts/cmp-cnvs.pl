#!/usr/bin/env perl
#
# Author: petr.danecek@sanger
#

use strict;
use warnings;
use Carp;

my $opts = parse_params();
read_dat($opts);
unify_regions($opts);
report_stats($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print 
        "Usage: cmp-cnvs.pl [OPTIONS] <dir1> <dir2> ..\n",
        "Options:\n",
        "   -q, --quality <float>       Filter CNVs by quality\n",
        "   -h, -?, --help              This help message.\n",
        "\n";
    exit -1;
}


sub parse_params
{
    my $opts = { qual=>0 };
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( $arg eq '-q' || $arg eq '--quality' ) { $$opts{qual} = shift(@ARGV); next; }
        if ( -d $arg ) { push @{$$opts{dirs}},$arg; next; }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( !$$opts{dirs} ) { error(); }
    return $opts;
}

sub read_dat
{
    my ($opts) = @_;
    my $regs = {};
    for my $dir (@{$$opts{dirs}})
    {
        open(my $fh,'<',"$dir/summary.tab") or error("$dir/summary.tab: $!");
        my $query;
        while (my $line=<$fh>)
        {
            if ( $line=~/^# RG/ )
            {
                chomp($line);
                my @items = split(/\t/,$line);
                $query   = $items[4];
                my $control = $items[5];
                $query   =~ s/^.*Copy number://;
                $control =~ s/^.*Copy number://;
                if ( exists($$opts{control}) && $$opts{control} ne $control ) { error("Different controls: $control vs $$opts{control}\n"); }
                if ( exists($$opts{query}{$query}) ) { error("Duplicate sample names? $query\n"); }
                $$opts{query}{$query} = 1;
                $$opts{control} = $control;
                next;
            }
            if ( $line=~/^#/ ) { next; }
            my ($rg,$chr,$start,$end,$cnq,$cnc,$qual) = split(/\t/,$line);
            chomp($qual);
            push @{$$opts{dat}{$query}{$chr}}, { start=>$start, end=>$end, cnc=>$cnc, cnq=>$cnq, qual=>$qual };
            $$regs{$chr}{$start} = 1;
        }
        close($fh);
        for my $chr (keys %{$$opts{dat}{$query}})
        {
            filter_by_quality($$opts{dat}{$query}{$chr}, $$opts{qual});
        }
    }
    $$opts{regs} = $regs;
}

sub filter_by_quality
{
    my ($regs,$qual) = @_;
    for my $reg (@$regs)
    {
        if ( $$reg{qual}>=$qual ) { next; }
        if ( $$reg{cnc} ne $$reg{cnq}  ) { $$reg{cnq} = $$reg{cnc}; }   # low quality
    }
}

sub unify_regions
{
    my ($opts) = @_;
    for my $chr (keys %{$$opts{regs}})
    {
        my $regs = $$opts{regs}{$chr};
        for my $pos (sort {$a<=>$b} keys %$regs)
        {
            for my $smpl (keys %{$$opts{query}})
            {
                my $dat = $$opts{dat}{$smpl}{$chr};
                my $i;
                for ($i=0; $i<@$dat; $i++)
                {
                    if ( $pos >= $$dat[$i]{start} && $pos<$$dat[$i]{end} ) { last; }
                    if ( $pos < $$dat[$i]{start} ) { last; }
                }
                if ( $i==0 && $pos < $$dat[$i]{start} )
                {
                    # one more block at the beginning
                    my $to = $$dat[0]{start};
                    splice(@$dat,0,0, { start=>$pos, end=>$to, cnc=>$$dat[0]{cnc}, cnq=>$$dat[0]{cnq}, qual=>$$dat[0]{qual} });
                    next;
                }

                if ( $i==@$dat ) { error("Uh: i=$i  pos=$pos  last=$$dat[$i-1]{start},$$dat[$i-1]{end}  $smpl,$chr\n"); }
                if ( $$dat[$i]{start}==$pos ) { next; }

                # one more block in the middle or at the end
                my $to = $$dat[$i]{end};
                $$dat[$i]{end} = $pos;
                splice(@$dat,$i+1,0, { start=>$pos, end=>$to, cnc=>$$dat[$i]{cnc}, cnq=>$$dat[$i]{cnq}, qual=>$$dat[$i]{qual} });
            }
        }
    }
}

sub is_continuation
{
    my ($regs,$i) = @_;
    if ( $i==0 ) { return 0; }
    if ( $$regs[$i-1]{cnc} ne $$regs[$i]{cnc} ) { return 0; }
    if ( $$regs[$i-1]{cnq} ne $$regs[$i]{cnq} ) { return 0; }
    if ( $$regs[$i-1]{end} ne $$regs[$i]{start} ) { return 0; }
    return 1;
}

sub is_shared
{
    my ($dat,$smpl,$chr,$i) = @_;
    for my $a (keys %$dat)
    {
        if ( $a eq $smpl ) { next; }
        if ( $$dat{$a}{$chr}[$i]{cnc}!=$$dat{$smpl}{$chr}[$i]{cnc} ) { next; }
        if ( $$dat{$a}{$chr}[$i]{cnq}!=$$dat{$smpl}{$chr}[$i]{cnq} ) { next; }
        return 1;
    }
    return 0;
}

sub report_stats
{
    my ($opts) = @_;
    my $tot_lens    = {};
    my $diff_lens   = {};
    my $shared_lens = {};
    my $ndiffs      = {};
    for my $smpl (keys %{$$opts{dat}})
    {
        $$ndiffs{$smpl} = 0;
        $$diff_lens{$smpl} = 0;
        $$shared_lens{$smpl} = 0;
        for my $chr (keys %{$$opts{dat}{$smpl}})
        {
            my $regs = $$opts{dat}{$smpl}{$chr};
            for (my $i=0; $i<@$regs; $i++)
            {
                my $reg = $$regs[$i];

                $$tot_lens{$smpl} += $$reg{end} - $$reg{start};
                if ( $$reg{cnc} ne $$reg{cnq} && $$reg{qual}>=$$opts{qual} ) 
                { 
                    $$diff_lens{$smpl} += $$reg{end} - $$reg{start}; 
                    if ( !is_continuation($regs,$i) ) { $$ndiffs{$smpl}++; }
                    if ( is_shared($$opts{dat},$smpl,$chr,$i) ) { $$shared_lens{$smpl} += $$reg{end} - $$reg{start}; }
                }
            }
        }
    }

    my @smpls = sort keys %{$$opts{query}};

    print "# SM, Samples\n";
    print "# ND, Number of different regions\n";
    print "# LD, Length of different regions (Mbp)\n";
    print "# SD, Length of shared differences (Mbp)\n";
    printf "# RG, List of regions, [2]Chromosome, [3]Start, [4]End, [5]Length, [6]Quality, [7-%d]CN of the samples\n", 7+@smpls;

    print "SM\t$$opts{control}\t" . join("\t",@smpls), "\n";
    print "ND\t0";
    for (my $i=0; $i<@smpls; $i++) { printf "\t%d", $$ndiffs{$smpls[$i]}; }
    print "\n";
    print "LD\t0";
    for (my $i=0; $i<@smpls; $i++) { printf "\t%.0f", $$diff_lens{$smpls[$i]}/1e6; }
    print "\n";
    print "SD\t0";
    for (my $i=0; $i<@smpls; $i++) { printf "\t%.0f", $$shared_lens{$smpls[$i]}/1e6; }
    print "\n";

    my @regs;
    my $smpl = $smpls[0];
    for my $chr (keys %{$$opts{dat}{$smpl}})
    {
        my $regs = $$opts{dat}{$smpl}{$chr};
        for (my $i=0; $i<@$regs; $i++)
        {
            my $reg = $$regs[$i];
            my $is_diff = 0;
            my $qual;
            for my $a (@smpls)
            {
                my $areg = $$opts{dat}{$a}{$chr}[$i];
                if ( $$areg{cnc} eq $$areg{cnq} ) { next; }
                $is_diff = 1;
                if ( !defined $qual or $qual<$$areg{qual} ) { $qual = $$areg{qual}; }
            }
            if ( !$is_diff ) { next; }
            if ( !defined $qual ) { $qual = $$reg{qual}; }
            if ( $qual < $$opts{qual} ) { next; }

            my $row = "RG\t$chr\t$$reg{start}\t$$reg{end}";
            $row .= sprintf "\t%.1f\t$qual\t$$reg{cnc}", ($$reg{end}-$$reg{start})/1e6;
            for my $smpl (@smpls)
            {
                $row .= "\t$$opts{dat}{$smpl}{$chr}[$i]{cnq}";
            }
            $row .= "\n";
            push @regs, [ $$reg{end}-$$reg{start}, $row];
        }
    }
    for my $rg (sort {$$b[0] <=> $$a[0]} @regs)
    {
        print $$rg[1];
    }
}



