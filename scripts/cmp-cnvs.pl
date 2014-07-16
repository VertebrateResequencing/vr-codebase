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
    if ( !$$opts{dirs} or scalar @{$$opts{dirs}}<2 ) { error(); }
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
            push @{$$opts{dat}{$query}{$chr}}, [$start,$end,$cnc,$cnq,$qual];
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
        if ( $$reg[4]>=$qual ) { next; }
        if ( $$reg[2] ne $$reg[3]  ) { $$reg[3] = $$reg[2]; }
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
                    if ( $pos >= $$dat[$i][0] && $pos<$$dat[$i][1] ) { last; }
                    if ( $pos < $$dat[$i][0] ) { last; }
                }
                if ( $i==0 && $pos < $$dat[$i][0] )
                {
                    # one more block at the beginning
                    my $to = $$dat[0][0];
                    splice(@$dat,0,0,[$pos,$to,$$dat[0][2],$$dat[0][3],$$dat[0][4]]);
                    next;
                }

                if ( $i==@$dat ) { error("Uh: i=$i  pos=$pos  last=$$dat[$i-1][0],$$dat[$i-1][1]  $smpl,$chr\n"); }
                if ( $$dat[$i][0]==$pos ) { next; }

                # one more block in the middle or at the end
                my $to = $$dat[$i][1];
                $$dat[$i][1] = $pos;
                splice(@$dat,$i+1,0,[$pos,$to,$$dat[$i][2],$$dat[$i][3],$$dat[$i][4]]);
            }
        }
    }
}

sub is_continuation
{
    my ($regs,$i) = @_;
    if ( $i==0 ) { return 0; }
    if ( $$regs[$i-1][2] ne $$regs[$i][2] ) { return 0; }
    if ( $$regs[$i-1][3] ne $$regs[$i][3] ) { return 0; }
    if ( $$regs[$i-1][1] ne $$regs[$i][0] ) { return 0; }
    return 1;
}

sub is_shared
{
    my ($dat,$smpl,$chr,$i) = @_;
    for my $a (keys %$dat)
    {
        if ( $a eq $smpl ) { next; }
        if ( $$dat{$a}{$chr}[$i][2]!=$$dat{$smpl}{$chr}[$i][2] ) { next; }
        if ( $$dat{$a}{$chr}[$i][3]!=$$dat{$smpl}{$chr}[$i][3] ) { next; }
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

                $$tot_lens{$smpl} += $$reg[1] - $$reg[0];
                if ( $$reg[2] ne $$reg[3] ) 
                { 
                    $$diff_lens{$smpl} += $$reg[1] - $$reg[0]; 
                    if ( !is_continuation($regs,$i) ) { $$ndiffs{$smpl}++; }
                    if ( is_shared($$opts{dat},$smpl,$chr,$i) ) { $$shared_lens{$smpl} += $$reg[1] - $$reg[0]; }
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
                if ( $$areg[2] eq $$areg[3] && $$areg[2] eq $$reg[2] ) { next; }
                $is_diff = 1;
                if ( !defined $qual or $qual<$$areg[4] ) { $qual = $$areg[4]; }
            }
            if ( !$is_diff ) { next; }
            if ( !defined $qual ) { $qual = $$reg[4]; }

            my $row = "RG\t$chr\t$$reg[0]\t$$reg[1]";
            $row .= sprintf "\t%.1f\t$qual\t$$reg[2]", ($$reg[1]-$$reg[0])/1e6;
            for my $smpl (@smpls)
            {
                $row .= "\t$$opts{dat}{$smpl}{$chr}[$i][3]";
            }
            $row .= "\n";
            push @regs, [ $$reg[1]-$$reg[0], $row];
        }
    }
    for my $rg (sort {$$b[0] <=> $$a[0]} @regs)
    {
        print $$rg[1];
    }
}



