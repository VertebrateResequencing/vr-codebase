#!/usr/bin/env perl

use strict;
use warnings;

open(my $pipe, "| ssh etch-dev64") || die "failed to open pipe\n";

print $pipe "svn update /software/vertres/bin\n";
print $pipe "svn update /software/vertres/modules\n";
print $pipe "exit\n";

close($pipe) || die "failed to close pipe\n";

exit;
