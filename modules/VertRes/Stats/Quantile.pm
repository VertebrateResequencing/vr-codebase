# Author: petr.danecek@sanger
#

=head1 NAME

Quantile.pm.  Single pass quantile estimation from Numerical Recipes 8.5.2, p. 435

=head1 SYNOPSIS

use VertRes::Stats::Quantile;

my $qest = Quantile->new();
my @test = ();
for (my $i=0; $i<1_000_000; $i++)
{
    my $val = rand(1000);
    push @test, $val;
    $qest->add($val);
}
@test = sort { $a<=>$b } @test;
for my $p qw(0.1 1 25 50 75 99 99.9)
{
    my $q = $p/100.;
    my $exact  = $test[ int($q * @test) ];
    my $approx = $qest->report($q);
    printf "%.1f \t exact=%.3f \t approx=%.3f \t err=%.3f\n", $p,$exact,$approx,$exact-$approx;
}

=cut


package Quantile;

use strict;
use warnings;
use Carp;

sub new
{
    my ($class,@args) = @_;
    my $self = @args ? {@args} : {};
    bless $self, ref($class) || $class;

    # init
    my @pval;
    for (my $i=85; $i<166; $i++) { $pval[$i] = ($i-75)/100.; }
    for (my $i=84; $i>=0; $i--)
    {
        $pval[$i] = 0.87191909 * $pval[$i+1];
        $pval[250-$i] = 1 - $pval[$i];
    }
    $$self{pval} = \@pval;
    $$self{dbuf} = [];
    $$self{qile} = [];
    for (my $i=0; $i<@pval; $i++) { $$self{qile}[$i] = 0; }
    $$self{nbuf} = 1_000;
    $$self{nq}   = scalar @pval;
    $$self{nt}   = 0;

    return $self;
}

sub add
{
    my ($self,$val) = @_;

    # update extremes
    if ( !exists($$self{q0}) or $val < $$self{q0} ) { $$self{q0} = $val; }
    if ( !exists($$self{qm}) or $val > $$self{qm} ) { $$self{qm} = $val; }

    my $dbuf = $$self{dbuf};
    push @$dbuf, $val;
    if ( @$dbuf == $$self{nbuf} ) { $self->update(); }
}

sub report
{
    my ($self,$p) = @_;

    if ( @{$$self{dbuf}} ) { $self->update(); }
    my $qile = $$self{qile};
    my $pval = $$self{pval};
    my $nq = $$self{nq};
    my $jl = 0;
    my $jh = $$self{nq} - 1;
    while ( $jh - $jl > 1 )
    {
        my $j = ($jh + $jl) >> 1;
        if ( $p > $$pval[$j] ) { $jl = $j; }
        else { $jh = $j; }
    }
    my $j = $jl;

    # interpolate
    my $q = $$qile[$j] + ($$qile[$j+1]-$$qile[$j])*($p-$$pval[$j])/($$pval[$j+1]-$$pval[$j]);
    if ( $$qile[$nq-1] < $q ) { $q = $$qile[$nq-1]; }
    return $$qile[0] > $q ? $$qile[0] : $q;
}

sub update
{
    my ($self) = @_;

    my $nt = $$self{nt};
    my $nd = @{$$self{dbuf}};
    my $nq = $$self{nq};
    my $dbuf = $$self{dbuf};
    @$dbuf = sort { $a<=>$b } @$dbuf;

    my $qile = $$self{qile};
    my @newqile;   # new quantiles after update
    for (my $i=0; $i<$nq; $i++) { $newqile[$i] = 0; }
    $$qile[0]     = $$self{q0};
    $newqile[0]   = $$self{q0};
    $$qile[$nq-1] = $$self{qm};
    $newqile[-1]  = $$self{qm};
    my $qold      = $$self{q0};
    my $qnew      = $$self{q0};
    my $tnew      = 0;
    my $told      = 0;
    my $jd = 0;
    my $jq = 1;

    my $pval = $$self{pval};
    $$pval[0]     = 0.5/($nt+$nd) < 0.5*$$pval[1] ? 0.5/($nt+$nd) : 0.5*$$pval[1];
    $$pval[$nq-1] = 1-0.5/($nt+$nd) > 0.5*(1+$$pval[$nq-2]) ? 1-0.5/($nt+$nd) : 0.5*(1+$$pval[$nq-2]);

    # loop over target p-values for interpolation
    for (my $iq=1; $iq<$nq-1; $iq++)
    {
        my $tgt = ($nt+$nd)*$$pval[$iq];

        # Find a succession of abcissa-ordinate pairs (qnew,tnew) that
        # are the discontinuity of value or slope and break to perform
        # an interpolation as we cross each target
        if ( $tnew < $tgt )
        {
            while (1)
            {
                if ( $jq<$nq && ($jd >= $nd or $$qile[$jq] < $$dbuf[$jd]) )
                {
                    # found slope discontinuity from old cdf
                    $qnew = $$qile[$jq];
                    $tnew = $jd + $nt * $$pval[$jq++];
                    if ( $tnew >= $tgt ) { last; }
                }
                else
                {
                    # found value discontinuity from batch cdf
                    $qnew = $$dbuf[$jd];
                    $tnew = $told;

                    if ( $$qile[$jq] > $$qile[$jq-1] )
                    {
                        $tnew += $nt * ($$pval[$jq] - $$pval[$jq-1]) * ($qnew - $qold)/($$qile[$jq] - $$qile[$jq-1]);
                    }
                    $jd++;
                    if ( $tnew >= $tgt ) { last; }
                    $told = $tnew++;
                    $qold = $qnew;
                    if ( $tnew >= $tgt ) { last; }
                }
                $told = $tnew;
                $qold = $qnew;
            }
        }
        # new interpolation
        $newqile[$iq] = $tnew == $told ? 0.5* ($qold+$qnew) : $qold + ($qnew - $qold)*($tgt-$told)/($tnew-$told);
        $told = $tnew;
        $qold = $qnew;
    }
    @{$$self{qile}} = @newqile;
    $$self{nt} += $nd;
    @{$$self{dbuf}} = ();
}


1;
