#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;

$ENV{PERL5LIB} = '';

my $cwd = getcwd;

# keep svn and git in sync
my $base = '/software/vertres';
my $origin_base = "$base/git/origin";
foreach my $repo ("$origin_base/vr-bin-external", "$origin_base/vr-codebase-svn") {
    chdir($repo);
    warn "syncing $repo\n";
    system("git svn rebase; git svn dcommit");
}

# update the checkouts that PATH and PERL5LIB point to
foreach my $repo ("$base/bin-git", "$base/codebase") {
    chdir($repo);
    warn "updating $repo\n";
    system("git pull");
}

# for now, we'll also update the defunct svn checkouts that nothing should be
# using
foreach my $repo ("$base/bin", "$base/scripts", "$base/modules") {
    chdir($repo);
    warn "updating $repo\n";
    system("svn up");
}

chdir($cwd);

exit;
