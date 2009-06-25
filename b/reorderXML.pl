#!/usr/bin/env perl

# Args: file containing ordering specs of the form
#       PARENT CHILD1 CHILD2 ...
# and one or more XML files.

use FileHandle;
use strict;

my $specFile = shift @ARGV;
my %ORDER;
my $fh = new FileHandle($specFile);
while (<$fh>) {
    chomp;
    tr/,/ /;
    my ($parent,@children) = split;
    foreach my $i (0 .. $#children) {
	$ORDER{$parent}{$children[$i]} = $i;
    }
}

my @lines;
while (<>) {
    chomp;
    my ($spaces,$content) = (/^(\s*)(.*)$/);
    $content =~ s/\s+$//;
    my $type = ($content =~ /^<\/([-A-Z0-9_]+)/) ? 'END' :
	($content =~ /^<([-A-Z0-9_]+).*\/>$/) ? 'WHOLE' :
	($content =~ /^<([-A-Z0-9_]+).*<\/.*>$/) ? 'WHOLE' :
	($content =~ /^<([-A-Z0-9_]+)/) ? 'BEGIN' : 'TEXT';
    my $tag = $1;
    push @lines,[length($spaces),$type,$tag,$content];
}
foreach my $pTree (readTrees(\@lines,0,$#lines)) {
    printTree($pTree);
}

sub readTrees {
    my ($pLines,$min,$max) = @_;
    my @trees;
    my $lo=$min;
    while ($lo <= $max) {
	my ($ind,$type,$tag,$cont) = @{$pLines->[$lo]};
	if ($type eq 'TEXT' || $type eq 'WHOLE') {
	    push @trees,$pLines->[$lo];
	    $lo ++;
	} elsif ($type eq 'END') {
	    die("");
	} else {
	    my $hi=$lo+1;
	    while ($hi <= $max) {
		if ($pLines->[$hi][0] == $ind) {
		    if ($pLines->[$hi][1] eq 'END' && $pLines->[$hi][2] eq $tag) {
			my @subs = readTrees($pLines,$lo+1,$hi-1);
			push @trees,[@{$pLines->[$lo]},reorder($tag,@subs)];
			$lo = $hi+1;
			last;
		    } else {
			die();
		    }
		}
		$hi ++;
	    }
	}
    }
    return @trees;
}

sub reorder {
    my ($parent,@trees) = @_;
    return sort {$ORDER{$parent}{$a->[2]} <=> $ORDER{$parent}{$b->[2]}} @trees;
}

sub printTree {
    my ($pTree) = @_;
    my ($ind,$type,$tag,$content,@children) = @$pTree;
    print ' ' x $ind;
    print "$content\n";
    if ($type eq 'BEGIN') {
	foreach my $pChild (@children) {
	    printTree($pChild);
	}
	print ' ' x $ind;
	print "</$tag>\n";
    }
}
