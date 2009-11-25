# This example pipeline contains two actions named "hello" and "world".
# Note: This code is not maintained and may not work with latest version 
#   of Pipeline.pm.
#

package VertRes::Pipelines::TrackDummy;
use base qw(VertRes::Pipeline);

use strict;
use warnings;

our @actions =
(
    {
        'name'     => 'hello',
        'action'   => \&hello,
        'requires' => \&hello_requires, 
        'provides' => \&hello_provides,
    },
    {
        'name'     => 'world',
        'action'   => \&world,
        'requires' => \&world_requires, 
        'provides' => \&world_provides,
    },
);

our $options = 
{
    'Hello' => 'Hello',
    'World' => 'World',
};


# --------- OO stuff --------------

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,@args);
    return $self;
}


# --------- hello --------------

sub hello_requires
{
    my ($self,$lane) = @_;
    my @requires = ();
    return \@requires;
}

sub hello_provides
{
    my ($self,$lane) = @_;
    my @provides = ('hello.txt');
    return \@provides;
}

sub hello
{
    my ($self,$lane_path,$action_lock) = @_;
    open(my $fh, '>', "$lane_path/hello.txt") or Utils::error("$lane_path/hello.txt: $!");
    print $fh $$self{'Hello'}, "\n";
    close $fh;

    return $$self{'Yes'};
}


# --------- world --------------

sub world_requires
{
    my ($self,$lane) = @_;
    my @requires = ('hello.txt');
    return \@requires;
}

sub world_provides
{
    my ($self,$lane) = @_;
    my @provides = ('world.txt');
    return \@provides;
}

sub world
{
    my ($self,$lane_path,$action_lock) = @_;
    open(my $fh, '>', "$lane_path/world.txt") or Utils::error("$lane_path/world.txt: $!");
    print $fh $$self{'World'}, "\n";
    close $fh;

    return $$self{'Yes'};
}


1;

