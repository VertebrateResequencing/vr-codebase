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

# rsync things from 'main' /software/vertres to precise-only /software/vertres
my @to_rsync = qw(/software/vertres/etc/ /software/vertres/vrpipe/archive_disc_pool);
foreach my $path (@to_rsync) {
    warn "\nrsyncing $path\n";
    system("rsync -rlptOD -z --delete $path precise-dev64:$path");
}

exit;

