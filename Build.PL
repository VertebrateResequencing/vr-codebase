#!/usr/bin/env perl
# Build script mainly used to run tests.
#
# Do this once, first:
# perl Build.PL
# Then before each commit, run the test for the module
# you altered:
# ./Build test --verbose --test_files t/appropriatetestscript.t
# And if you did something major, run the whole test suite:
# ./Build test
#
# less often, when interested:
# ./Build testcover
# ./Build testpod
# ./Build testpodcoverage
#
# setenv DEVEL_COVER_OPTIONS '+ignore,5.8.8'
# beforehand if coverage includes extraneous stuff
#

use strict;
use warnings;
use lib 'build_lib';
use MyBuild;


my $build = MyBuild->new(
    module_name          => 'VertRes',
    dist_version         => 0,
    dist_author          => 'Vertebrate Resequencing group at the Sanger Institute',
    dist_abstract        => 'A collection of modules and scripts for processing and analysing large quantities of sequencing data',
    license              => 'perl',
    recursive_test_files => 1,
    pm_files             => get_pm_files(),
    script_files         => 'scripts'
);

$build->create_vrconfig_file();
$build->create_build_script;

exit;

sub get_pm_files {
	my %pm_files;
	foreach my $module (check_dir('modules')) {
		my $in_lib = $module;
		$in_lib =~ s/^modules/lib/;
		$pm_files{$module} = $in_lib;
	}
	return \%pm_files;
}

sub check_dir {
        my $dir = shift;
        opendir(my $dir_handle, $dir);
	
	my @pm_files;
        foreach my $thing (readdir($dir_handle)) {
                if ($thing =~ /^\.+$/) { next; }
                $thing = $dir."/".$thing;

                if (-d $thing) {
                        push(@pm_files, check_dir($thing));
                        next;
                }
		
		if ($thing =~ /\.pm$/) {
			push(@pm_files, $thing);
		}
        }
	return @pm_files;
}
