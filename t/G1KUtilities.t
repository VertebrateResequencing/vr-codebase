#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
	use Test::Most tests => 5;
	
	use_ok('G1KUtilities');
	use_ok('Cwd');
}

is G1KUtilities::path2Gender, 'unknown', 'when not in the hierarchy at all, path2Gender gives us unknown';

my $cwd = getcwd();

chdir('/nfs/sf8/G1K/DATA/Trio-CEU/NA12878');
is G1KUtilities::path2Gender, 'female', 'path2Gender tells us NA12878 is female - this is UNCONFIRMED!';
chdir('/nfs/sf8/G1K/DATA/Trio-CEU/NA12891');
is G1KUtilities::path2Gender, 'male', 'path2Gender tells us NA12891 is male - this is UNCONFIRMED!';

chdir($cwd);

exit;