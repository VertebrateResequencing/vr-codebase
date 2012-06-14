package Pathogens::Reports::Mapping::Report ;

use Moose;
use Pathogens::Reports::Mapping::Row ;
use Pathogens::Reports::Mapping::Spreadsheet ;
use VRTrack::VRTrack;

has 'filehandle' => ( is => 'ro', isa => 'FileHandle',                                 required   => 1); # output filehandle 
has 'vrtrack'    => ( is => 'ro', isa => 'VRTrack::VRTrack',                           required   => 1); # database object
has 'lanes'      => ( is => 'ro', isa => 'ArrayRef[VRTrack::Lane]',                    required   => 1); # lanes array
has '_rows'      => ( is => 'ro', isa => 'ArrayRef[Pathogens::Reports::Mapping::Row]', lazy_build => 1); # spreadsheet rows

sub _build__rows
{
    my($self)= @_;

    my @rows;

    # sort lanes
    my @sorted_lanes = sort _sort_lanes @{$self->lanes};

    # make row objects for each lane
    for my $lane (@sorted_lanes)
    {
	# get mapstats for lane sorted by id
	my @sorted_mapstats  = sort {$a->row_id <=> $b->row_id} @{$lane->mappings};

	# create rows
	my $qc_row; # zero or one qc row 
	my @lane_rows = (); # zero or more mapping rows
	for my $map (@sorted_mapstats)
	{
	    my $row = Pathogens::Reports::Mapping::Row->new(vrtrack => $self->vrtrack, lane => $lane, mapstats => $map);	    
	    if($row->is_qc_mapstats)
	    {
		$qc_row = $row; # select qc row with the highest mapstats_id
	    }
	    else
	    {
		push @lane_rows, $row; # mapping rows
	    } 
	}
	
	# add qc fields to mapping data
	if(defined $qc_row)
	{
	    for my $row (@lane_rows)
	    {
		$row->transfer_qc_values($qc_row);
	    }
	}

	# push rows
	push @rows, $qc_row if defined $qc_row; # single qc row
	push @rows, @lane_rows if scalar @lane_rows; # multiple mapping rows
    }
   
    return \@rows;
}

# Sort routine for multiplexed lane names (eg 1234_5#6)
# Runs are sorted in decending order (ie Latest first.)
# Lane and Tag are sorted in ascending order.
sub _sort_lanes
{
    my @a = split(/\_|\#/,$a->name());
    my @b = split(/\_|\#/,$b->name());

    $b[0]<=>$a[0] || $a[1]<=>$b[1] || $a[2]<=>$b[2];
}

sub output_csv
{
    my ($self) = @_;

    # Output data.
    my $spreadsheet = Pathogens::Reports::Mapping::Spreadsheet->new( filehandle => $self->filehandle, rows => $self->_rows );
    $spreadsheet->output_csv;
    return 1;
}

1;
