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

#private function to validate the lanes are eligible for submission
sub _validateLanes
{
	my $self = shift;
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
	
	my $originally = @laneNames;
	
	#connect to the database
	my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate($$self{'database'}, 'r');
	
	for(my $i = 0; $i < @laneNames; $i ++ )
	{
		my $lane = VRTrack::Lane->new_by_name( $vrtrack, $_ );
		
		if( ( defined( $lane->submission_id() ) && $lane->submission_id() > 0 ) || ( $lane->qc_status() ne 'passed' ) )
		{
			delete( $laneNames[ $i ] );
			$i --;
		}
	}
	
	print scalar( @laneNames ).' out of '.$originally.' lanes are eligible for submission\n';
	
	$$self{'lanes'} = \@laneNames;
}

sub _gatherAllUnsubmittedLanes
{
	my ($self, $project) = @_;
	
	my @projects;
	
	#connect to the database
	my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate($$self{'database'}, 'r');
	
	if( defined( $project ) )
	{
		my $p = VRTrack::Project->new_by_name($vrtrack, $project);
		$self->throw("Cant find project: $project") unless defined( $p );
		push( @projects, $p );
	}
	else
	{
		@projects = @{ $vrtrack->projects() };
	}
	
	my @laneNames;
	foreach( @projects )
	{
		my $samples = $project->samples();
		foreach( @{$samples} )
		{
			my $sample = $_;
			my $libraries = $sample->libraries();
			
			foreach( @{$libraries} )
			{
				my $library = $_;
				my $lanes = $library->lanes();
				
				foreach( @{$lanes} )
				{
					my $lane = $_;
					if( ! defined( $lane->submission_id() ) && $lane->qc_status() eq 'passed' )
					{
						push( @laneNames, $lane->name() );
					}
				}
			}
		}
	}
	
	$$self{'lanes'} = \@laneNames;
}

=head2 printLanes

 Title   : printLanes
 Usage   : my $obj = VertRes::Submissions::SRA->printLanes();
 Function: Print the list of lanes to be submitted
 Returns : n/a
 Args    : n/a

=cut
sub printLanes()
{
	my $self = shift;
	
	my @lanes = @{ $$self{'lanes'} };
	foreach( @lanes )
	{
		print $_."\n";
	}
}

1;
