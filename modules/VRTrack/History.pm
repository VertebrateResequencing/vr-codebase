package VRTrack::History;

=head1 NAME

VRTrack::History - Sequence Tracking History object

=head1 SYNOPSIS

use VRTrack::History;

my $hist = VRTrack::History->new();

# create an instance of one of the Core_obj-inheriting classes as normal, eg:
my $vrlane = VRTrack::Lane->new_by_name($vrtrack, 'foo');

# supply it to one of the methods of this class. Eg. to get the date the lane's
# library was last changed:
my $datetime = $hist->state_change($vrlane, 'library_id');

# or the date the lane's processed flag was last set to 'swapped'
$datetime = $hist->was_processed($vrlane, 'swapped');

# Now make all Core_obj-inheriting classes return the most recent versions that
# are older than our datetime (ie. view the database as it was immediately
# prior to the event associated with our datetime):
$hist->time_travel($datetime);

# ... do stuff with new Core_obj instances, which represent the old db state

# Revert back to normal behaviour (latest version):
$hist->time_travel('latest');

=head1 DESCRIPTION

This module lets you choose a particular point in time based on a desired
state of a particular database entry, then sets the Core_obj api to return
the latest version of objects not younger than that time point.

Eg. if at some point a lane swapped libraries, you can effectively return to
the state of the database just before it was swapped, to discover what the
original library was.

=head1 AUTHOR

sb10@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);
use VRTrack::Core_obj;


=head2 new

  Example    : my $obj = $class->new()
  Description: Make a new VRTrack::History object
  Returntype : $class object

=cut

sub new {
    my $class = shift;
    my $self = {};
    bless $self, ref($class) || $class;
    return $self;
}


=head2 historical_objects

  Arg [1]    : VRTrack::Core_obj-inheriting instance
  Example    : my @objs = $class->historical_objects($vrlane);
  Description: Get a list of instances of all versions of your Core_obj
               throughout history.
  Returntype : list of Core_obj-inheriting instances (one of which will
               correspond to Arg [1], ordered oldest to latest)

=cut

sub historical_objects {
    my ($self, $core_obj) = @_;
    ($core_obj && ref($core_obj) && $core_obj->isa('VRTrack::Core_obj')) || confess "A VRTrack::Core_obj is required";
    my $vrtrack = $core_obj->vrtrack;
    my $id = $core_obj->id;
    
    my @row_ids = $core_obj->row_ids;
    my @objs;
    foreach my $row_id (@row_ids) {
	push(@objs, $core_obj->new($vrtrack, $id, $row_id));
    }
    
    return @objs;
}


=head2 state_change

  Arg [1]    : VRTrack::Core_obj-inheriting instance
  Arg [2]    : name of a method of Arg[1] - the state to look for being changed
  Example    : my $datetime = $class->state_change($vrlane, 'library_id');
  Description: Get the datetime corresponding to the last time that a certain
               state changed.
  Returntype : datetime formatted string (or 'latest' if the state never
               changed)

=cut

sub state_change {
    my ($self, $core_obj, $method) = @_;
    confess "A Core_obj and one of its methods must be supplied" unless $core_obj && $method && $core_obj->can($method);
    
    my @objs = reverse($self->historical_objects($core_obj));
    my $latest = shift(@objs);
    my $latest_state = $latest->$method;
    my $changed = 'latest';
    foreach my $obj (@objs) {
	my $state = $obj->$method;
	if ("$state" ne "$latest_state") {
	    $changed = $obj->changed;
	    last;
	}
    }
    
    return $changed;
}


=head2 was_processed

  Arg [1]    : VRTrack::Core_obj-inheriting instance
  Arg [2]    : name of an allowed processed flag
  Example    : my $datetime = $class->was_processed($vrlane, 'swapped');
  Description: Get the datetime corresponding to the last time that a certain
               processed flag was set. Arg[1] must have an is_processed()
	       method.
  Returntype : datetime formatted string (or 'latest' if the flag was never
               set)

=cut

sub was_processed {
    my ($self, $core_obj, $flag) = @_;
    confess "A Core_obj that supports is_processed() is required" unless $core_obj && $core_obj->can('is_processed');
    
    my @objs = reverse($self->historical_objects($core_obj));
    my $changed = 'latest';
    foreach my $obj (@objs) {
	my $processed = $obj->is_processed($flag);
	if ($processed) {
	    $changed = $obj->changed;
	    last;
	}
    }
    
    return $changed;
}


=head2 time_travel

  Arg [1]    : datetime string|'latest'
  Example    : $obj = $class->time_travel('2010-01-04 10:49:10');
  Description: Set all Core_obj inheriting classes to return new instances as if
               we had travelled back in time to immediately prior to the given
	       datetime. Revert back to normal behaviour by setting the 'latest'
	       keyword.
  Returntype : n/a

=cut

sub time_travel {
    my $self = shift;
    VRTrack::Core_obj->global_history_date(@_);
}

1;
