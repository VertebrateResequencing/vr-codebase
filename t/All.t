#!/usr/bin/perl -w
# test basics of all perl modules and scripts in the distribution;
use strict;
use warnings;

use Test::Strict;

all_perl_files_ok('modules', 'bin');

exit;
