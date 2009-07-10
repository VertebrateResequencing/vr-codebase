#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 3;
    
    use_ok('VertRes::Wrapper::WrapperI');
}

my $wi = VertRes::Wrapper::WrapperI->new();
isa_ok $wi, 'VertRes::Base';

# quick basic test
package VertRes::Wrapper::Head;
use base qw(VertRes::Wrapper::WrapperI);
sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args,
                                  exe      => 'head',
                                  params   => [qw(n lines)],
                                  switches => [qw(silent help)]);
    $self->run_method('open');
    return $self;
}
package main;
my $wrapper = VertRes::Wrapper::Head->new(n => 5, quiet => 1);
my $fh = $wrapper->run($0);
my @lines = <$fh>;
is @lines, 5, 'looks like head -n 5 worked';

# needs a zillion more tests!...

exit;
