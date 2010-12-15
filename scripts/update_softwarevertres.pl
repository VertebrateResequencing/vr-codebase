#!/usr/bin/env perl

use strict;
use warnings;

#open(my $pipe, "| ssh etch-dev64") || die "failed to open pipe\n";

system("svn update /software/vertres/bin");
system("svn update /software/vertres/scripts");
system("svn update /software/vertres/modules");
#print $pipe "exit\n";

#close($pipe) || die "failed to close pipe\n";

exit;
