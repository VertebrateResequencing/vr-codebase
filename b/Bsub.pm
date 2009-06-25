package Bsub;

use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(expandAndSubmit bsubInvoke bsub isRunning jobSucceeded waitUntilDone unlinkLogs);

use Utility;
use File::Basename;
use Carp;
use Cwd;

use strict;

# See SGRP.pm for an example of how to call this 

sub expandAndSubmit {
    my ($pOpt,$pInfo,@args) = @_;
    my @expanded = splitArgsAndReport(@args);
    if (defined $pInfo->{expansionRoot}) {
	@expanded = expandFromRoot($pInfo,@args);
    }
    my $printCmd = (@expanded < 30);
    foreach my $pA (@expanded) {
	expandAndSubmit1($pOpt,$pInfo,$printCmd,@$pA);
    }
}

sub expandFromRoot {
    my ($pInfo,$cmd,$anaType,@args) = @_;
    my $eRoot = $pInfo->{expansionRoot};
    my $na = max($#args,($pInfo->{nFixed}{$cmd} || $pInfo->{nFixedDefault})-2_);
    my @ans;
    $anaType = 'default' if ($anaType eq '.' || !defined $anaType);
    foreach my $i (0 .. $na) {
	$args[$i] = '*' if ($i > $#args || $args[$i] eq '.');
    }
    my $pat = join("/",$eRoot,@args);
    foreach my $dir (glob($pat)) {
	my @comps = split(/\//,$dir);
	push @ans,[$cmd,$anaType,@comps[-scalar(@args)..-1]];
    }
    return @ans;
}

sub splitArgsAndReport {
    my (@args) = @_;
    return reportExpansions(splitArgs(@args));
}

# Print expansions if there were any:

sub reportExpansions {
    my (@expanded) = @_;
    if (@expanded > 1) {
	my $base = basename($0);
	foreach my $pA (@expanded) {
	    printf("Expansion: %s %s\n",$base,join(" ",@$pA));
	}
    }
    return @expanded;
}

# Characters "," and ":" inside an argument induce a split;
# we try the args resulting from each combination, except that
# we always choose the same position for "," args (and so each
# "," arg must have the same number of components). Example:
#
# sgrp.pl -b esnp cere,para meth1:meth2 CERE,PARA 0 1 2 3 4 5 6 7 8 9
# Expansion: cere meth1 CERE 0 1 2 3 4 5 6 7 8 9
# Expansion: cere meth2 CERE 0 1 2 3 4 5 6 7 8 9
# Expansion: para meth1 PARA 0 1 2 3 4 5 6 7 8 9
# Expansion: para meth2 PARA 0 1 2 3 4 5 6 7 8 9

sub splitArgs {
    my @args = @_;
    # First look for commas
    my ($nc,@args2);
    foreach my $a (@args) {
	my @aa = split(/,/,$a);
	push @args2,\@aa;
	my $n = scalar(@aa);
	if ($n > 1) {
	    if (defined $nc) {
		confess("Comma count mismatch at $a") unless ($n == $nc);
	    } else {
		$nc = $n;
	    }
	}
    }
    if (defined $nc) {
	my @ans;
	foreach my $k (0 .. $nc-1) {
	    my @ak;
	    foreach my $pA (@args2) {
		my $i = (@$pA > 1) ? $k : 0;
		push @ak,$pA->[$i];
	    }
	    push @ans,splitArgs1(0,@ak);
	}
	return @ans;
    } else {
	return splitArgs1(0,@args); # nothing with commas
    }
}

sub splitArgs1 {
    my ($kMin,@args) = @_;
    # Look for first position with a colon from kMin upwards.
    # Also look for ranges of numbers e.g. 3-12
    foreach my $k ($kMin .. $#args) {
	my $ak = $args[$k];
	my @aa;
	if ($ak =~ /^(\d+)-(\d+)$/ && $1 <= $2) {
	    @aa = ($1 .. $2);
	} else {
	    @aa = split(/:/,$ak);
	}
	if (@aa > 1) {
	    my @ans;
	    foreach my $i (0 .. $#aa) {
		my @ak;
		foreach my $j (0 .. $#args) {
		    push @ak,($j==$k) ? $aa[$i] : $args[$j];
		}
		push @ans,splitArgs1($k+1,@ak);
	    }
	    return @ans;
	}
    }
    return (\@args);
}

sub expandAndSubmit1 {
    my ($pOpt,$pInfo,$printCmd,$cmd,@args) = @_;
    if ($pOpt->{b}) {
	my $nf = $pInfo->{nFixed}{$cmd};
	$nf = $pInfo->{nFixedDefault} unless (defined $nf);
	my $force = $pOpt->{f};
	my $nSubs = $pOpt->{n} || 1;
	my $cwd = cwd();
	chdir($pInfo->{rootDir});
	my ($fcp,$queue,$req) = $pInfo->{decide}($cmd,@args);
	bsubInvoke($force,$nSubs,$cmd,$nf,$fcp,$queue,$req,@args);
	chdir($cwd);
    } elsif (!$cmd) {
	print "Usage: $0 [-b] [-r rootDir] ...\n";
	my $pm = $pInfo->{instructionModule};
	if ($pm) {
	    my $fh = openToRead($pm);
	    while (<$fh>) {
		if (/^##/) {
		    print substr($_,2);
		}
	    }
	}
	die("");
    } elsif (defined $pInfo->{command}{$cmd}) {
	my $shCmd = $pInfo->{shellCommand};
	print ">> $shCmd $cmd @args\n" if ($printCmd);
	&{$pInfo->{command}{$cmd}}(@args);
    } elsif (defined $pInfo->{commandSet}{$cmd}) {
	foreach my $cmd2 (split(/ /,$pInfo->{commandSet}{$cmd})) {
	    $pInfo->{callback}($pOpt,$cmd2,@args);
	}
    } else {
	confess("unrecognized command: $cmd");
    }
}

# Interface for SGRP-like calls. Note we expect to have chdir's to
# the right place BEFORE calling this.

sub bsubInvoke {
    my ($force,$nSubs,$cmd,$nf,$fixedCmdPfx,$queue,$req,@args) = @_;
    $nf = $#args if ($nf < 0); # -1: ensure one job exactly
    my @ids;
    my @fixed = @args[0 .. $nf-1];
    my $fixedCmd = join(" ",($fixedCmdPfx,@fixed));
    my $fixedKey = join(".",($cmd,@fixed));
    my $bSubDir = sprintf("bout/%s/%s",$fixed[2] || 'ANY',$cmd);
    mkdir_p($bSubDir);
  ARG:
    foreach my $k ($nf .. $#args) {
	my $key = "$fixedKey.$args[$k]";
	my $stem = "$bSubDir/$key";
	my $tf = "$stem.touch";
	my $fullCmd = "$fixedCmd $args[$k]; rm -f $tf";
	foreach my $suf (qw(bout berr touch)) {
	    my $f = "$stem.$suf";
	    if (-f $f) {
		if ($force) {
		    moveAside($f);
		} else {
		    print "$f already exists\n";
		    next ARG;
		}
	    }
	}
	touch($tf);
	for (my $n=0; $n<$nSubs; $n++) {
	    my $id = bsub($stem,$fullCmd,undef,$queue,$key,$req);
	    push @ids,$id;
	    my $msg = "$key: submitted job $id";
	    if ($nSubs > 1) {
		$msg .= sprintf(" (%d of %d)",$n+1,$nSubs);
	    }
	    print "$msg\n";
	}
    }
    if ($cmd eq 'bugb') {
	print "When done, run: sgrp.pl cmgb ...\n";
    }
    return @ids;
}

# base is base for .bout and .berr files. Returns job number.

sub bsub {
    my ($base,$cmd,$pri,$q,$name,$req) = @_;
    die("Bad base, containing space: $base\n") if ($base =~ /\s/);
    my $sp = "";
    $sp = "-sp $pri" if ($pri);
    $sp .= " -q $q" if ($q);
    $sp .= " -J$name" if ($name);
    $sp .= " -R'$req'" if ($req);
    my $fh = openFromPipe("bsub $sp -o $base.bout -e $base.berr '$cmd'");
    while (<$fh>) {
	chomp;
	if (/^Job <(\d+)> /) {
	    $fh->close();
	    return $1;
	}
    }
    die("Cannot find job number for $base\n");
}

# Returns 1 if ANY of the jobs in list are running.
sub isRunning {
    my (@nList) = @_;
    my %toFind;
    foreach my $n (@nList) {
	$toFind{$n} = 1 if ($n >= 0);
    }
    return 0 unless (keys %toFind);
    my $fh = openFromPipe("bjobs 2>&1");
    my @found;
    while (<$fh>) {
	if (/^No unfinished/) {
	    return;
	} elsif (/^\s*(\d+)\s/) {
   	    push @found,$1 if (defined $toFind{$1});
        }
    }
    return \@found if (@found);
}

sub jobSucceeded {
    my ($base) = @_;
    my $f = "$base.bout";
    while (! -s $f) {
	sleep(10);
    }
    my $fh = openToRead($f);
    <$fh>; # Sender line
    my $line = <$fh>; # Subject line
    if ($line =~ /^Subject/) {
	return ($line =~ / Done$/); # as opposed to Exited
    } else {
	die("jobSucceeded: cannot parse $f\n");
    }
}

sub waitUntilDone {
    my (@nList) = @_;
    return unless (@nList); # there may be nothing to wait for...
    sleep(20); # to ensure all have registered...
    while (isRunning(@nList)) {
	sleep(10);
    }
}

sub unlinkLogs {
    my (@bases) = @_;
    foreach my $base (@bases) {
	foreach my $f ("$base.bout","$base.berr") {
	    unlink($f) if (-f $f);
	}
    }
}

1;
