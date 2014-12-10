
=head1 NAME

VertRes::Pipelines::Import_iRODS_cram - pipeline for importing cram files from iRODS


=head1 EXAMPLE CONFIG FILES:

=head2 pipeline.conf

__VRTrack_Import_cram__ importirods.conf

=head2 import_irods.conf

root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::Import_iRODS_cram',
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

A module for importing cram files into a tracking database. 


=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Pipelines::Import_iRODS_cram;
use VertRes::Pipeline;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use Data::Dumper;
use File::Copy;
use VertRes::Wrapper::iRODS;
use VertRes::LSF;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::File;
use VertRes::Utils::FileSystem;
use VertRes::Wrapper::fastqcheck;
use Pathogens::Import::ValidateFastqConversion;
use File::Basename;

our $actions = [

    # Create the hierarchy path, download the files.
    {
        'name'     => 'get_files',
        'action'   => \&get_files,
        'requires' => \&get_files_requires,
        'provides' => \&get_files_provides,
    },
    {
        'name'     => 'convert_to_fastq',
        'action'   => \&convert_to_fastq,
        'requires' => \&convert_to_fastq_requires,
        'provides' => \&convert_to_fastq_provides,
    },

    # If all files were downloaded OK, update the VRTrack database.
    {
        'name'     => 'update_db',
        'action'   => \&update_db,
        'requires' => \&update_db_requires,
        'provides' => \&update_db_provides,
    },

];
our %options = (
    bsub_opts  => '',
    cramtools_jar  => '/software/pathogen/external/apps/usr/share/java/cramtools-2.1.jar',
    cramtools_java => '/software/jdk1.8.0_11/bin/java',
    samtools_exec  => '/software/pathogen/external/apps/usr/local/samtools-1.1/samtools'
);

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new( %options, actions => $actions, @args );
    if ( defined( $self->{db} ) ) {
        $self->{vrtrack} = VRTrack::VRTrack->new( $self->{db} ) or $self->throw("Could not connect to the database\n");
    }
    $self->{fsu} = VertRes::Utils::FileSystem->new;

    return $self;
}

#---------- Convert to FASTQ ----------------

sub convert_to_fastq_requires {
    my ($self) = @_;
    return $self->get_files_provides;
}

sub convert_to_fastq_provides {
    my ( $self, $lane_path ) = @_;

    my @fastqs ;
    for my $file (@{$self->{files}})
    {
      for my $f (@{$self->_fastq_from_cram($file)})
      {
        push(@fastqs,$f );
      }
    }

    return \@fastqs;
}

sub convert_to_fastq {
    my ( $self, $lane_path, $lock_file ) = @_;
    my $memory_in_mb = 1000;

    my $prefix   = $$self{prefix};
    my $work_dir = $lane_path;
    
    my $file = $lane_path.'/'.$self->{files}->[0];
    
    open( my $fh, '>', "$work_dir/${prefix}convert_to_fastq.pl" ) or $self->throw("$work_dir/${prefix}convert_to_fastq.pl: $!");
    print $fh qq[
  use strict;
  use warnings;
  use VertRes::Pipelines::Import_iRODS_cram;
  
  my \$import = VertRes::Pipelines::Import_iRODS_cram->new();
  \$import->cram_to_fastq(qq[$file]);
  system('touch ${prefix}convert_to_fastq_done');
  ];

    close($fh);
    VertRes::LSF::run(
        $lock_file, $work_dir, "${prefix}convert_to_fastq",
        { bsub_opts => " -M${memory_in_mb} -R 'select[mem>$memory_in_mb] rusage[mem=$memory_in_mb]' " },
        qq[perl -w ${prefix}convert_to_fastq.pl]
    );

    return $$self{No};
}



sub _fastq_from_cram
{
  my ( $self, $cram ) = @_;
  my ( $filename, $dirs, $suffix ) = fileparse( $cram, '.cram' );
  return [$dirs.$filename.'_1.fastq.gz', $dirs.$filename.'_2.fastq.gz' ];
}


sub cram_to_fastq {
    my ( $self, $file ) = @_;

    my ( $filename, $dirs, $suffix ) = fileparse( $file, '.cram' );
    my $fastq_base = $dirs.$filename ;
    
    # remove supplementary and secondary alignments
    system($self->{samtools_exec}." view -F 0x900 -h -C -o output.cram $file");
    system("mv output.cram $file");
    
    my $cmd = $self->{cramtools_java} .' -jar '.$self->{cramtools_jar}. ' fastq --gzip -I '.$file .' --fastq-base-name '. $fastq_base ;

    system($cmd );
    my $fastqcheck = VertRes::Wrapper::fastqcheck->new();
    
    my @fastqcheck_files;
    for my $fastq (@{$self->_fastq_from_cram($file)})
    {
      push(@fastqcheck_files,$fastq . '.fastqcheck');
      $fastqcheck->run( $fastq, $fastq . '.fastqcheck' );
    }
    
    my $filename_without_path = $filename. $suffix;
    my $validate = Pathogens::Import::ValidateFastqConversion->new(
         fastqcheck_filenames => \@fastqcheck_files,
         irods_filename       => $filename_without_path
        );
    $self->throw("The number of reads in the FASTQ files doesnt match the number if reads in iRODS") unless($validate->is_total_reads_valid());
}

#---------- get_files ---------------------

# Requires nothing
sub get_files_requires {
    my ($self) = @_;
    return [];
}

sub get_files_provides {
    my ($self, $lane_path) = @_;
    my @provides;
    foreach my $file (@{$$self{files}}) {
        push @provides, File::Basename::basename($file);
    }
    @provides || $self->throw("Something went wrong; we don't seem to provide any bams!");
    return \@provides;
}

sub get_files {
    my ( $self, $lane_path, $lock_file ) = @_;

    my $prefix   = $$self{prefix};
    my $work_dir = $lane_path;
    my $opts     = $self->dump_opts(qw(files));
    `mkdir -p $lane_path`;
    
    open( my $fh, '>', "$work_dir/${prefix}import_files.pl" ) or $self->throw("$work_dir/${prefix}import_files.pl: $!");
    print $fh qq[
  use strict;
  use warnings;
  use VertRes::Pipelines::Import_iRODS_cram;

  my $opts

  my \$import = VertRes::Pipelines::Import_iRODS_cram->new(%\$opts);
  \$import->download_files();

  ];

    close($fh);
    VertRes::LSF::run( $lock_file, $work_dir, "${prefix}import_files", $self, qq[perl -w ${prefix}import_files.pl] );

    return $$self{No};
}

sub download_files {
    my ($self) = @_;

    my $irods = VertRes::Wrapper::iRODS->new();

    for my $file ( @{ $$self{files} } ) {
        next unless ( $file =~ /\.cram/ );
        my $ifile = $irods->find_file_by_name($file);
        
        if ( !defined $ifile ) { $self->warn("No such file in iRODS? [$ifile]\n"); next; }
        my ( $filename, $dirs, $suffix ) = fileparse($ifile);
        next if ( -e $filename );

        if ( !( $ifile =~ m{([^/]+)$} ) ) { $self->throw("FIXME: [$ifile]"); }
        my $outfile = $1;
        if ( -e $outfile ) { next; }

        $irods->get_file( $ifile, "$outfile.tmp" );

        # Get the md5sum and check
        my $md5 = $irods->get_file_md5($ifile);
        open( my $fh, '>', "$outfile.md5" ) or $self->throw("$outfile.md5: $!");
        print $fh "$md5  $outfile.tmp\n";
        close($fh);

        Utils::CMD(qq[md5sum --status -c $outfile.md5]);

        # Recreate the checksum file to contain the correct file name
        open( $fh, '>', "$outfile.md5" ) or $self->throw("$outfile.md5: $!");
        print $fh "$md5  $outfile\n";
        close($fh);

        move( $outfile . '.tmp', $outfile );
        chmod 0664, $outfile;

    }
}

# Requires the gzipped fastq files. How many? Find out how many .md5 files there are.
sub update_db_requires {
    my ( $self, $lane_path ) = @_;
    return $self->convert_to_fastq_provides;
}

# This subroutine will check existence of the key 'db'. If present, it is assumed
#   that Import should write the stats and status into the VRTrack database. In this
#   case, 0 is returned, meaning that the task must be run. The task will change the
#   QC status from NULL to pending, therefore we will not be called again.
#
#   If the key 'db' is absent, the empty list is returned and the database will not
#   be written.
#
sub update_db_provides {
    my ($self) = @_;
    if ( exists( $$self{db} ) ) { return 0; }
    my @provides = ();
    return \@provides;
}

sub update_db {
    my ( $self, $lane_path, $lock_file ) = @_;

    if ( !$$self{db} ) { $self->throw("Expected the db key.\n"); }

    my $prefix = $self->{prefix};

    # remove job files
    foreach my $file (qw(import_files convert_to_fastq )) {
        foreach my $suffix (qw(o e pl)) {
            unlink( $self->{fsu}->catfile( $lane_path, $prefix . $file . '.' . $suffix ) );
        }
    }
    
      # Remove Large Files
      my @cram_suffix   = ('cram','cram.md5');

      my $bam = $self->{files}->[0];

      for my $suffix (@cram_suffix)
      {
  	# Remove crams
    	Utils::CMD(qq[rm $lane_path/$$self{lane}.$suffix]) if(-e qq[$lane_path/$$self{lane}.$suffix]);
      }
    

    my $vrtrack = VRTrack::VRTrack->new( $$self{db} ) or $self->throw("Could not connect to the database\n");
    my $vrlane = VRTrack::Lane->new_by_name( $vrtrack, $$self{lane} ) or $self->throw("No such lane in the DB: [$$self{lane}]\n");

    $vrtrack->transaction_start();

    my $rawreads = 0;
    my $rawbases = 0;
    for my $fastq ( @{ $self->convert_to_fastq_provides } ) {
        my $parser = VertRes::Parser::fastqcheck->new( file => $lane_path.'/'.$fastq . '.fastqcheck' );
        $rawreads += $parser->num_sequences() || 0;
        $rawbases += $parser->total_length()  || 0;
    }

    $vrlane->is_processed( 'import', 1 );
    $vrlane->is_withdrawn(0);
    $vrlane->raw_reads($rawreads);
    $vrlane->raw_bases($rawbases);
    $vrlane->update();
    $vrtrack->transaction_commit();
    return $$self{Yes};
}

sub dump_opts {
    my ( $self, @keys ) = @_;
    my %opts;
    for my $key (@keys) {
        $opts{$key} = exists( $$self{$key} ) ? $$self{$key} : undef;
    }
    return Data::Dumper->Dump( [ \%opts ], ["opts"] );
}

1;

