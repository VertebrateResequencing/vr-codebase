# Author: petr.danecek@sanger
#

=head1 NAME

MannWhitney.pm.  Mann-Whitney U test

=head1 SYNOPSIS

use VertRes::Stats::MannWhitney;

my @a = qw(1 2 3 3 4);
my @b = qw(3 4 4 5 6 7);
my ($pval,$u1,$u2) = MannWhitney->Utest(\@a,\@b);

=cut


package MannWhitney;

use strict;
use warnings;
use Carp;

sub Utest
{
    my ($class,$vals1,$vals2) = @_;
    return calc_mwu($vals1,$vals2);
}

sub mann_whitney_1947
{
    my ($n, $m, $U) = @_;
    if ($U<0) { return 0; }
    if ($n==0||$m==0) { return $U==0 ? 1 : 0; }
    return $n/($n+$m)*mann_whitney_1947($n-1,$m,$U-$m) + $m/($n+$m)*mann_whitney_1947($n,$m-1,$U);
}

sub mann_whitney_1947_cdf
{
    my ($n, $m, $U) = @_;
    my $sum = 0;
    for (my $i=0; $i<=$U; $i++) 
    {
        $sum += mann_whitney_1947($n,$m,$i);
    }
    return $sum;
}

sub calc_mwu
{
    my ($avals, $bvals) = @_;

    my @avals = sort { $a<=>$b } @$avals;
    my @bvals = sort { $a<=>$b } @$bvals;

    my $na = scalar @avals;
    my $nb = scalar @bvals;
    my $ia = 0;
    my $ib = 0;
    my $Ua = 0;
    my $xa = 0;
    my $xb = 0;
    while ( $ia<$na && $ib<$nb )
    {
        my $va = $avals[$ia]; 
        if ( !$xa ) { $xa = 1; }
        while ( $ia+1<$na && $avals[$ia]==$avals[$ia+1] ) { $ia++; $xa++; }

        my $vb = $bvals[$ib]; 
        while ( $ib+1<$nb && $bvals[$ib]==$bvals[$ib+1] ) { $ib++; $xb++; }

        if ( $va==$vb ) { $ia++; $Ua += $xa*(($ib-$xb) + (1+$xb)*0.5); $ib++; $xa = 0; $xb = 0; }
        elsif ( $vb<$va ) { $ib++; $xb = 0; }
        else { $ia++; $Ua += $xa*($ib-$xb); $xa = 0; }
    }
    while ( $ia<$na ) { $Ua += $nb; $ia++; }

    my $Ub = ($na * $nb) - $Ua;
    my $U_min = $Ua < $Ub ? $Ua : $Ub;

    my $pval;
    if ( $na==1 ) 
    { 
        $pval = 2.0 * (floor($U_min)+1) / ($nb+1); 
    }
    elsif ( $nb==1 ) 
    { 
        $pval = 2.0 * (floor($U_min)+1) / ($na+1); 
    }
    elsif ( $na>=8 || $nb>=8 )
    {
        # Normal approximation, very good for na>=8 && nb>=8 and reasonable if na<8 or nb<8
        my $mean = ($na*$nb)*0.5;
        my $var2 = ($na*$nb)*($na+$nb+1)/12.0;
        my $z = ($U_min - $mean)/sqrt(2*$var2);   # z is N(0,1)
        $pval = 2.0-erfc($z);  # which is 1 + erf(z)
    }
    else
    {
        # Exact calculation
        $pval = 2*mann_whitney_1947_cdf($na,$nb,$U_min);
        if ( $pval>1 ) { $pval = 1; }
    }
    return ($pval,$Ua,$Ub);
}

# complementary error function
# \frac{2}{\sqrt{\pi}} \int_x^{\infty} e^{-t^2} dt
# AS66, 2nd algorithm, http://lib.stat.cmu.edu/apstat/66
# 
sub erfc
{
    my ($x) = @_;
    my $p0 = 220.2068679123761;
    my $p1 = 221.2135961699311;
    my $p2 = 112.0792914978709;
    my $p3 = 33.912866078383;
    my $p4 = 6.37396220353165;
    my $p5 = .7003830644436881;
    my $p6 = .03526249659989109;
    my $q0 = 440.4137358247522;
    my $q1 = 793.8265125199484;
    my $q2 = 637.3336333788311;
    my $q3 = 296.5642487796737;
    my $q4 = 86.78073220294608;
    my $q5 = 16.06417757920695;
    my $q6 = 1.755667163182642;
    my $q7 = .08838834764831844;
    my $z = abs($x) * 2**0.5;
    if ($z > 37.) { return $x > 0.? 0. : 2.; }
    my $expntl = exp($z**2* -.5);
    my $p;
    if ($z < 10. / 2**0.5)
    {
        # for small z
        $p = $expntl * (((((($p6 * $z + $p5) * $z + $p4) * $z + $p3) * $z + $p2) * $z + $p1) * $z + $p0)
            / ((((((($q7 * $z + $q6) * $z + $q5) * $z + $q4) * $z + $q3) * $z + $q2) * $z + $q1) * $z + $q0);
    }
    else 
    {
        $p = $expntl / 2.506628274631001 / ($z + 1. / ($z + 2. / ($z + 3. / ($z + 4. / ($z + .65)))));
    }
    return $x > 0.? 2. * $p : 2. * (1. - $p);
}



1;
