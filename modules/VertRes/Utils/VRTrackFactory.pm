=head1 NAME

VertRes::Utils::VRTrackFactory - factory class for getting VRTrack objects

=head1 SYNOPSIS

use VertRes::Utils::VRTrackFactory;

my $fsu = VertRes::Utils::VRTrackFactory->instanstiate('mouse_reseq_track', 'r');

=head1 DESCRIPTION

A simple factory class that returns VRTrack objects to centralise the database connection information

=head1 AUTHOR

Thomas Keane tk2@sanger.ac.uk

=cut

package VertRes::Utils::VRTrackFactory;
use base qw(VertRes::Base);

use strict;
use warnings;
use Carp;

use VRTrack::VRTrack;

my $HOST = 'mcs4a';
my $PORT = 3306;
my $READ_USER = 'vreseq_ro';
my $WRITE_USER = 'vreseq_rw';
my $WRITE_PASS = 't3aml3ss';

=head2 new

 Title   : instantiate
 Usage   : my $vrtrack = VRTrackFactory->instantiate( 'mouse', 'r' );
 Function: Ask the factory to return a VRTrack object to the database specified
 Returns : VRTrack object
 Args    : database name: a valid VRTrack database, mode: either 'r' or 'rw' connection

=cut
sub instantiate
{
	my $database = shift;
	my $mode = lc( shift );
	
	croak "Invalid connection mode (r or rw valid): $mode\n" unless $mode =~ /^r|rw$/;
	
	my $user = $mode =~ /^r$/ ? $READ_USER : $WRITE_USER;
	my $pass = $mode =~ /^r$/ ? '' : $WRITE_PASS;
	
	my $vrtrack = VRTrack::VRTrack->new({ host => 'mcs4a',
                                    port => 3306,
                                    user => $user,
                                    password => $pass,
                                    database => $database,
                                   });
    return $vrtrack;
}

1;
