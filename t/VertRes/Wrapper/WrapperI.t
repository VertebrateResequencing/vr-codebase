#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 7;
    
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

# arg handling tests
package VertRes::Wrapper::Foo;
use base qw(VertRes::Wrapper::WrapperI);
sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args,
                                  exe => 'foo');
    return $self;
}
sub action1 {
    my ($self, @args) = @_;
    
    $self->switches([qw(baz naz)]);
    $self->params([qw(foo bar)]);
    
    $self->_set_params_and_switches_from_args(@args);
}
sub action2 {
    my ($self, @args) = @_;
    
    $self->switches([qw(caz maz)]);
    $self->params([qw(goo car)]);
    
    $self->_set_params_and_switches_from_args(@args);
    $self->_set_params_string(join => '=');
}
package main;
$wrapper = VertRes::Wrapper::Foo->new();
$wrapper->action1(baz => 1, foo => 'boo');
is $wrapper->_get_params_string(), ' -foo boo -baz', '_get_params_string default test';
$wrapper->action2(maz => 1, car => 'far');
is $wrapper->_get_params_string(), ' car=far maz', '_get_params_string non-default test';

$wrapper->extras('--unsupported woo', '-new goo');
is $wrapper->_get_params_string(), ' -car far -maz --unsupported woo -new goo', '_get_params_string extras test, reverts to default behaviour too';

# needs a zillion more tests!...
TODO: {
    local $TODO = "Currently little more than a stub";
    ok 0, 'needs more tests...';
}

exit;
