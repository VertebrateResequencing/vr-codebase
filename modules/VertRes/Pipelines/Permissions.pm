
=head1 NAME

VertRes::Pipelines::Permissions -  Pipeline for setting file permissions by brute force

=head1 DESCRIPTION

Pipeline for setting file permissions by brute force. Very disk intensive

=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Pipelines::Permissions;

use strict;
use warnings;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::Library;
use VRTrack::Sample;
use VertRes::LSF;
use base qw(VertRes::Pipeline);
use VertRes::Utils::FileSystem;
use VertRes::Utils::Assembly;
use File::Spec;
use Utils;

our $actions = [
    {
        name     => 'permissions',
        action   => \&permissions,
        requires => \&permissions_requires,
        provides => \&permissions_provides
    }
];

our %options = ( bsub_opts => '' );

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new( %options, actions => $actions, @args );
    return $self;
}


=head2 permissions_requires

 Title   : permissions_requires
 Usage   : my $required_files = $obj->permissions_requires('/path/to/lane');
 Function: Find out what files the permissions action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub permissions_requires {
  my ($self) = @_;
  return [];
}

=head2 permissions_provides

 Title   : permissions_provides
 Usage   : my $provided_files = $obj->permissions_provides('/path/to/lane');
 Function: Find out what files the permissions action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub permissions_provides {
    my ($self) = @_;
    return ['file_which_will_never_exist'];
}

=head2 permissions

 Title   : permissions
 Usage   : $obj->permissions('/path/to/lane', 'lock_filename');
 Function: Unlink all the pipeline-related files (_*) in a lane, as well
           as the split directory.
 Returns : $VertRes::Pipeline::Yes
 Args    : lane path, name of lock file to use

=cut

sub permissions {
  my ($self, $lane_path, $action_lock) = @_;
  
  return unless(defined($$self{octal_permissions}));

  if(defined($$self{unix_group}) )
  {
    my $change_permissions_obj = Bio::VertRes::Permissions::ModifyPermissions->new(
        input_directories => [$lane_path],
        group             => $$self{unix_group},
  	threads           => 0,
        octal_permissions => $$self{octal_permissions});
    $change_permissions_obj->update_permissions;
  }	
  
  return $self->{Yes};
}

1;

