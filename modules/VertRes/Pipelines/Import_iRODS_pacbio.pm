=head1 NAME

VertRes::Pipelines::Import_iRODS_pacbio - pipeline for importing pacbio files from iRODS


=head1 EXAMPLE CONFIG FILES:

=head2 pipeline.conf

__VRTrack_Import_Pacbio__ importirods.conf

=head2 import_irods.conf

root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::Import_iRODS_pacbio',
prefix  => '_',

db =>
   {
        database => 'pathogen_example_track',
        host     => 'web-mii-shap',
        port     => 3303,
        user     => 'pathpipe_ro',
        password => '',
   },
data => 
   {  # Pipeline-specific data
        db  => 
            {
               database => 'pathogen_example_track',
               host     => 'web-mii-shap',
               port     => 3303,
               user     => 'pathpipe_rw',
               password => 'xxx',       
            },
   }, 


=head1 DESCRIPTION

A module for importing pacbio files into a tracking database. 


=head1 AUTHOR

path-help@sanger.ac.uk

=cut


package VertRes::Pipelines::Import_iRODS_pacbio;
use VertRes::Pipeline;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use Data::Dumper;
use VertRes::Wrapper::iRODS;
use VertRes::LSF;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::File;
use VertRes::Utils::FileSystem;

our $actions =
[
    # Create the hierarchy path, download and bamcheck the bam files.
    {
        'name'     => 'get_files',
        'action'   => \&get_files,
        'requires' => \&get_files_requires, 
        'provides' => \&get_files_provides,
    },

];
our %options = ( bsub_opts => '' );


sub new 
{
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new( %options, actions => $actions, @args );
  if ( defined( $self->{db} ) ) {
      $self->{vrtrack} = VRTrack::VRTrack->new( $self->{db} ) or $self->throw("Could not connect to the database\n");
  }
  $self->{fsu} = VertRes::Utils::FileSystem->new;

  return $self;
}

#---------- get_files ---------------------

# Requires nothing
sub get_files_requires
{
    my ($self) = @_;
    return [];
}

sub get_files_provides
{
    my ($self, $lane_path) = @_;
    return $$self{files};
}

sub get_files
{
      my ($self,$lane_path,$lock_file) = @_;
      
      my $prefix   = $$self{prefix};
      my $work_dir = $lane_path;
      my $opts = $self->dump_opts(qw(files));
      `mkdir -p $lane_path`;

      # Create a script to be run on LSF.
      open(my $fh,'>', "$work_dir/${prefix}import_files.pl") or $self->throw("$work_dir/${prefix}import_files.pl: $!");
      print $fh qq[
  use strict;
  use warnings;
  use VertRes::Pipelines::Import_iRODS_pacbio;

  my $opts

  my \$import = VertRes::Pipelines::Import_iRODS_pacbio->new(%\$opts);
  \$import->download_files();

  ];

      close($fh);
      VertRes::LSF::run($lock_file,$work_dir,"${prefix}import_files",$self,qq[perl -w ${prefix}import_files.pl]);

      return $$self{No};
}



sub download_files
{
    my ($self) = @_;

    my $irods = VertRes::Wrapper::iRODS->new();

    for my $ifile (@{$$self{files}})
    {
        next if($ifile =~ /\.bam/);
        if ( !defined $ifile ) { $self->warn("No such file in iRODS? [$ifile]\n"); next; }

        if ( !($ifile=~m{([^/]+)$}) ) { $self->throw("FIXME: [$ifile]"); }        
        my $outfile = $1;
        if ( -e $outfile ) { next; }

        $irods->get_file($ifile,"$outfile.tmp");
        chmod 0664,"$outfile.bai";

        # Get the md5sum and check
        my $md5 = $irods->get_file_md5($ifile);
        open(my $fh,'>',"$outfile.md5") or $self->throw("$outfile.md5: $!");
        print $fh "$md5  $outfile.tmp\n";
        close($fh);

        Utils::CMD(qq[md5sum --status -c $outfile.md5]);

        # Recreate the checksum file to contain the correct file name
        open($fh,'>',"$outfile.md5") or $self->throw("$outfile.md5: $!");
        print $fh "$md5  $outfile\n";
        close($fh);
    }
}

sub dump_opts
{
    my ($self,@keys) = @_;
    my %opts;
    for my $key (@keys)
    {
        $opts{$key} = exists($$self{$key}) ? $$self{$key} : undef;
    }
    return Data::Dumper->Dump([\%opts],["opts"]);
}

1;

