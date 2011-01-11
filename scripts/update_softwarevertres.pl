#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;

$ENV{PERL5LIB} = '';

my $cwd = getcwd;

my $base = '/software/vertres';

# update the checkouts that PATH and PERL5LIB point to
foreach my $repo ("$base/external-bin", "$base/codebase") {
    chdir($repo);
    warn "\nupdating $repo\n";
    system("git pull");
}

chdir($cwd);

exit;
