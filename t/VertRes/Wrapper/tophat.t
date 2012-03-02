#!/usr/bin/env perl
use strict;
use warnings;
use File::Copy;
use File::Spec;

BEGIN {
    use Test::Most tests => 36;
    
    use_ok('VertRes::Wrapper::tophat');
    use_ok('VertRes::Utils::Mappers::tophat');
}


my $tophat = VertRes::Wrapper::tophat->new(quiet => 1);
isa_ok $tophat, 'VertRes::Wrapper::WrapperI';
is $tophat->quiet, 1, 'quiet set via new';

is $tophat->exe, 'tophat', 'exe ok';
like $tophat->version, qr/\d+\.\d+/, 'version ok';

