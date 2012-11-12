#!/usr/bin/env perl

use strict;
use warnings;

my $base = '/software/vertres';

# update the checkouts that PATH and PERL5LIB point to
foreach my $repo ("$base/bin-external", "$base/codebase", "$base/vrpipe/master") {
    foreach my $server ('uk10k-login', 'precise-dev64') {
        warn "\nupdating $repo on $server\n";
        system(qq[ssh $server "umask 002; cd $repo; git checkout .; git pull"]);
    }
}

exit;

