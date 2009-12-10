=head1 NAME

VertRes::Utils::SRA - SRA submission utility functions

=head1 SYNOPSIS

use VertRes::Utils::SRA;

my $dbsnp_util = VertRes::Utils::SRA->new();

=head1 DESCRIPTION

General utility functions for making dbSNP submission files - 1 method per dbSNP submission section

=head1 AUTHOR

Thomas Keane: thomas.keane@sanger.ac.uk

=cut

package VertRes::Utils::SRA;

use strict;
use warnings;

use VRTrack::VRTrack;

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Submissions::SRA->new(database=>'database_name',lanes=>'lanes file', holduntil=>'YYYY-MM-DD', project=>'project');
 Function: Create a new VertRes::Submissions::SRA object.
 Returns : VertRes::Submissions::SRA object
 Args    : n/a

=cut

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
	
	if ( !exists($$self{'database'}) ) { $self->throw("Expected a database name to be provided!");}
	if ( exists($$self{'lanes'}) && ! -f $$self{'lanes'} ) {$self->throw("Cant find file of lane names: ".$$self{'lanes'})}
	
	if( exists($$self{'lanes'}) )
	{
		_validateLanes();
	}
	else
	{
		_gatherAllUnsubmittedLanes($$self{'project'});
	}
	
    return $self;
}

sub _validateLanes
{
	#read in the lane names
	my @laneNames;
	open( my $fh, $$self{'lanes'} ) or $self->throw('Cant open file of lane names: '.$$self{'lanes'});
	while( <$fh> )
	{
		chomp;
		
		#small sanity check
		$self->throw("Invalid lane name: ".$_) unless $_ !~ /^\d+_\d+$/;
		
		push( @laneNames, $_ );
	}
	close( $fh );
	
	#connect to the database
	my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate($$self{'database'}, 'r');
	
	my @lanes;
	foreach( @laneNames )
	{
		my $lane = VRTrack::Lane->new_by_name( $vrtrack, $_ );
		
		if( defined( $lane->submission_id() ) && $lane->submission_id() > 0 )
		{
			push( @lanes, $lane );
		}
	}
	
	$$self{'laneObjects'} = \@lanes;
}

sub _gatherAllUnsubmittedLanes
{
	
}

1;
