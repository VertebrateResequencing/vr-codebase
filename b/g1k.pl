#!/usr/bin/env perl

=head1 NAME

g1k.pl

=head1 SYNOPSIS

g1k.pl action cond proj indiv lib lane

=head1 DESCRIPTION

Thousand genomes (g1k) project analysis pipeline invocation script.
Call without arguments for help string.

=head1 AUTHOR

David Carter, C<dmc@sanger.ac.uk>

=head1 MAINTAINER

Jim Stalker<jws@sanger.ac.uk>
Thomas Keane<tk2@sanger.ac.uk>

=cut

# G1K/*.pm and other *.pm files are in here

use File::Basename;

use lib dirname($0);

use G1K::G1K;
use Getopt::Std;

use strict;

# This pipeline is derived from the one used in SGRP; we don't
# really use these switches yet...
my %OPT;
getopts('bfn:r:i:',\%OPT);

# In G1K/G1K.pm:
g1k(\%OPT,@ARGV);
