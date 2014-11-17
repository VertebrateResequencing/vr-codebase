#!/usr/bin/env perl

use strict;
use warnings;

my $base = '/software/vertres';

# update the checkouts that PATH and PERL5LIB point to
foreach my $repo ("$base/bin-external", "$base/codebase", "$base/vrpipe/master", "$base/update_pipeline", "$base/conf") {
    foreach my $server ('vr-login') {
        warn "\nupdating $repo on $server\n";
        system(qq[ssh -A $server "umask 002; cd $repo; git checkout .; git pull"]);
    }
}

# rsync things from precise /software/vertres to lenny-only /software/vertres
#my @to_rsync = qw(/software/vertres/etc/ /software/vertres/vrpipe/archive_disc_pool /software/vertres/vrpipe/scratch107_disc_pool /software/vertres/vrpipe/siteconfig/);
#foreach my $path (@to_rsync) {
#    warn "\nrsyncing $path\n";
#    system("rsync -rlptOD -z --delete $path uk10k-login:$path");
#}

exit;

