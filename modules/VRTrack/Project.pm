package VRTrack::Project;

=head1 NAME

VRTrack::Project - Sequence Tracking Project object

=head1 SYNOPSIS
    my $proj = VRTrack::Project->new($vrtrack, $project_id);

    #get arrayref of sample objects in a project
    my $samples = $project->samples();
    
    my $id = $project->id();
    my $name = $project->name();

=head1 DESCRIPTION

An object describing the tracked properties of a project.

=head1 AUTHOR

jws@sanger.ac.uk (author)

=head1 METHODS

=cut

use strict;
use warnings;
use Carp qw(cluck confess);
use VRTrack::Sample;
use VRTrack::Study;

use base qw(VRTrack::Core_obj
            VRTrack::Hierarchy_obj
	    VRTrack::Named_obj
	    VRTrack::SequenceScape_obj);


=head2 fields_dispatch

  Arg [1]    : none
  Example    : my $fieldsref = $proj->fields_dispatch();
  Description: Returns hashref dispatch table keyed on database field
               Used internally for new and update methods
  Returntype : hashref

=cut

sub fields_dispatch {
    my $self = shift;
    
    my %fields = %{$self->SUPER::fields_dispatch()};
    %fields = (%fields,
               project_id        => sub { $self->id(@_)},
               ssid              => sub { $self->ssid(@_)},
               hierarchy_name    => sub { $self->hierarchy_name(@_)},
               study_id          => sub { $self->study_id(@_)},
			   data_access_group => sub { $self->data_access_group(@_)},
	       name              => sub { $self->name(@_)});
    
    return \%fields;
}

###############################################################################
# Class methods
###############################################################################


=head2 new_by_name

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : project name
  Example    : my $project = VRTrack::Project->new_by_name($vrtrack, $name);
  Description: Class method. Returns latest Project object by name and project_id.  If no such name is in the database, returns undef
  Returntype : VRTrack::Project object

=cut


=head2 new_by_hierarchy_name

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : project hierarchy_name
  Example    : my $project = VRTrack::Project->new_by_hierarchy_name($vrtrack, $hierarchy_name)
  Description: Class method. Returns latest Project object by hierarchy_name.  If no such hierarchy_name is in the database, returns undef.  Dies if multiple hierarchy_names match.
  Returntype : VRTrack::Project object

=cut


=head2 new_by_ssid

  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : project sequencescape id
  Example    : my $project = VRTrack::Project->new_by_ssid($vrtrack, $ssid);
  Description: Class method. Returns latest Project object by ssid.  If no such ssid is in the database, returns undef
  Returntype : VRTrack::Project object

=cut


=head2 create
  
  Arg [1]    : vrtrack handle to seqtracking database
  Arg [2]    : name
  Example    : my $file = VRTrack::Project->create($vrtrack, $name)
  Description: Class method.  Creates new Project object in the database.
  Returntype : VRTrack::Project object
   
=cut


=head2 is_name_in_database

  Arg [1]    : project name
  Arg [2]    : hierarchy name
  Example    : if(VRTrack::Project->is_name_in_database($vrtrack, $name, $hname)
  Description: Class method. Checks to see if a name or hierarchy name is already used in the project table.
  Returntype : boolean

=cut


###############################################################################
# Object methods
###############################################################################

=head2 id

  Arg [1]    : id (optional)
  Example    : my $id = $proj->id();
               $proj->id('104');
  Description: Get/Set for ID of a project
  Returntype : Internal ID integer

=cut


=head2 hierarchy_name

  Arg [1]    : directory name (optional)
  Example    : my $hname = $project->hierarchy_name();
  Description: Get/set project hierarchy name.  This is the directory name (without path) that the project will be named in a file hierarchy.
  Returntype : string

=cut


=head2 name

  Arg [1]    : name (optional)
  Example    : my $name = $proj->name();
               $proj->name('1000Genomes-B1-TOS');
  Description: Get/Set for project name
  Returntype : string

=cut


=head2 ssid

  Arg [1]    : ssid (optional)
  Example    : my $ssid = $proj->ssid();
               $proj->ssid(104);
  Description: Get/Set for project SequenceScape ID
  Returntype : string

=cut


=head2 changed

  Arg [1]    : changed (optional)
  Example    : my $changed = $project->changed();
               $project->changed('20080810123000');
  Description: Get/Set for project changed
  Returntype : string

=cut


=head2 samples

  Arg [1]    : None
  Example    : my $samples = $project->samples();
  Description: Returns a ref to an array of the sample objects that are associated with this project
  Returntype : ref to array of VRTrack::Sample objects

=cut

sub samples {
    my $self = shift;
    return $self->_get_child_objects('VRTrack::Sample');
}


=head2 sample_ids

  Arg [1]    : None
  Example    : my $sample_ids = $project->sample_ids();
  Description: Returns a ref to an array of the sample IDs that are associated with this project
  Returntype : ref to array of integer sample IDs

=cut

sub sample_ids {
    my $self = shift;
    return $self->_get_child_ids('VRTrack::Sample');
}


=head2 add_sample

  Arg [1]    : sample name
  Example    : my $newproj = $track->add_sample('NOD mouse 1');
  Description: create a new sample, and if successful, return the object
  Returntype : VRTrack::Sample object

=cut

sub add_sample {
    my ($self, $sname) = @_;
    # TODO: if ssid is defined, then it should also not be added twice
    return $self->_add_child_object('new_by_name_project', 'VRTrack::Sample', $sname, $self->id);
}


=head2 get_sample_by_name

  Arg [1]    : sample name
  Example    : my $sample = $track->get_sample_by_name('My sample');
  Description: retrieve sample object by name
  Returntype : VRTrack::Sample object

=cut

sub get_sample_by_name {
    my $self = shift;
    return $self->_get_child_by_field_value('samples', 'name', @_);
}


=head2 get_sample_by_id

  Arg [1]    : sample id 
  Example    : my $sample = $proj->get_sample_by_id(1154);
  Description: retrieve sample object by internal id
  Returntype : VRTrack::Sample object

=cut

sub get_sample_by_id {
    my $self = shift;
    return $self->_get_child_by_field_value('samples', 'id', @_);
}


=head2 get_sample_by_ssid

  Arg [1]    : sample sequencescape id
  Example    : my $sample = $proj->get_sample_by_ssid(1154);
  Description: retrieve sample object by sequencescape id
  Returntype : VRTrack::Sample object

=cut

sub get_sample_by_ssid {
    my $self = shift;
    return $self->_get_child_by_field_value('samples', 'ssid', @_);
}


=head2 study_id

  Arg [1]    : study_id (optional)
  Example    : my $study_id = $proj->study_id();
               $proj->study_id(1);
  Description: Get/Set for project internal study_id
  Returntype : integer

=cut

sub study_id {
    my $self = shift;
    return $self->_get_set('study_id', 'number', @_);
}


=head2 study

  Arg [1]    : study accession (optional)
  Example    : my $study = $proj->study();
               $proj->study('SRP000031');
  Description: Get/Set for project study.  Lazy-loads study object from $self->study_id.  If a study accession is supplied, then study_id is set to the corresponding study in the database. If no such study exists, returns undef.  Use add_study to add a study in this case.
  Returntype : VRTrack::Study object

=cut

sub study {
    my $self = shift;
    return $self->_get_set_child_object('get_study_by_acc', 'VRTrack::Study', @_);
}


=head2 add_study

  Arg [1]    : study acc
  Example    : my $ind = $proj->add_study('NA19820');
  Description: create a new study, and if successful, return the object
  Returntype : VRTrack::Study object

=cut

sub add_study {
    my $self = shift;
    return $self->_create_child_object('get_study_by_acc', 'VRTrack::Study', @_);
}


=head2 get_study_by_acc

  Arg [1]    : study_name
  Example    : my $ind = $proj->get_study_by_acc('NA19820');
  Description: Retrieve a VRTrack::Study object by name
  Returntype : VRTrack::Study object

=cut

sub get_study_by_acc {
    my ($self, $acc) = @_;
    return VRTrack::Study->new_by_acc($self->{vrtrack}, $acc);
}


=head2 descendants

  Arg [1]    : none
  Example    : my $desc_objs = $obj->descendants();
  Description: Returns a ref to an array of all objects that are descendants of this object
  Returntype : arrayref of objects

=cut

sub _get_child_methods {
    return qw(samples);
}

=head2 data_access_group

  Arg [1]    : data_access_group
  Example    : my $num_bp = $file->data_access_group();
	       $file->data_access_group('abc');
  Description: Get/Set for data_access_group
  Returntype : string

=cut

sub data_access_group {
    my $self = shift;
    return $self->_get_set('data_access_group', 'string', @_);
}

=head2 main_data_access_group

  Arg [1]    : main_data_access_group
  Example    : my $num_bp = $file->main_data_access_group();
  Description: Get main_data_access_group
  Returntype : string

=cut

sub main_data_access_group {
    my $self = shift;
	return undef unless(defined($self->data_access_group()));
	my @data_access_groups = split(/\s/,$self->data_access_group);
	return $data_access_groups[0];
}

1;
