=head1 NAME

VertRes::Utils::Scaffold - wrapper for calling scaffolding

=head1 SYNOPSIS

=head1 DESCRIPTION

Provides a uniform interface for scaffolding so you dont need to worry about the difference in the underlying applications

=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Utils::Scaffold;

use strict;
use warnings;
use VertRes::IO;
use File::Basename;
use VertRes::Utils::Hierarchy;
use VertRes::Utils::FileSystem;
use VRTrack::Lane;

use base qw(VertRes::Base);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Scaffold->new();
 Function: Create a new VertRes::Utils::Scaffold object.
 Returns : VertRes::Utils::Scaffold object
 Args    : scaffolder => 'sspace'

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 find_module

 Title   : find_module
 Usage   : my $module = $obj->find_module();
 Function: Find out what assembly utility module to use
 Returns : class string (call new on it)

=cut

sub find_module {
    my ($self) = @_;
    return "VertRes::Utils::Scaffolders::".$self->{scaffolder};
}








1;
