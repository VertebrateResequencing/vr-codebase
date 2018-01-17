
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
use File::Copy;
use VertRes::Wrapper::iRODS;
use VertRes::LSF;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::File;
use VertRes::Utils::FileSystem;
use VertRes::Wrapper::fastqcheck;
use File::Basename;
use File::Path qw(remove_tree);

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
    bsub_opts  => ''
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
    return $self->_bas_h5_filenames;
}

sub convert_to_fastq_provides {
    my ( $self, $lane_path ) = @_;

    my @output_files ;
    for my $file (@{$self->_output_fastq_filenames})
    {
      push(@output_files,$file);
      push(@output_files,$file.'.fastqcheck');
    }
	for my $file(@{$self->_output_bam_filenames})
	{
		push(@output_files,$file);
	}
    push(@output_files,'_convert_to_fastq_done');
    return \@output_files;
}

sub _output_fastq_filenames {
    my ($self) = @_;
    my @output_files;
    for my $file ( @{ $self->_bas_h5_filenames } ) {
        my ( $filename, $dirs, $suffix ) = fileparse( $file, '.bas.h5' );

        my $fastq_filename = $filename . '.fastq.gz';
        push( @output_files, $fastq_filename );
    }
    return \@output_files;
}

sub _output_bam_filenames {
    my ($self) = @_;
    my @output_files;
    for my $file ( @{ $self->_bas_h5_filenames } ) {
        my ( $filename, $dirs, $suffix ) = fileparse( $file, '.bas.h5' );
        push( @output_files, $filename . '.subreads.bam' );
		push( @output_files, $filename . '.scraps.bam' );
    }
    return \@output_files;
}

sub _bas_h5_filenames {
    my ($self) = @_;
    my @output_files;
    my $file_regex = $self->{lane_path}."/".'*.bas.h5';
    my @files = glob( $file_regex);
	
    for my $file ( @files ) {
        my ( $filename, $dirs, $suffix ) = fileparse($file);
        push( @output_files, $filename );
    }
    return \@output_files;
}

sub _h5_filenames {
    my ($self) = @_;
    my @output_files;
	
    my $file_regex = $self->{lane_path}."/".'*.h5';
    my @files = glob( $file_regex);
	
    for my $file ( @files ) {
        my ( $filename, $dirs, $suffix ) = fileparse($file);
        push( @output_files, $filename );
    }
    return \@output_files;
}

sub convert_to_fastq {
    my ( $self, $lane_path, $lock_file ) = @_;
    my $memory_in_mb = 2500;
	my $threads = 8;

    my $prefix   = $$self{prefix};
    my $work_dir = $lane_path;
    my $opts     = Data::Dumper->Dump( [ \@{ $self->_bas_h5_filenames } ], ["bas_files"] );

    open( my $fh, '>', "$work_dir/${prefix}convert_to_fastq.pl" ) or $self->throw("$work_dir/${prefix}convert_to_fastq.pl: $!");
    print $fh qq[
  use strict;
  use warnings;
  use VertRes::Pipelines::Import_iRODS_pacbio;

  my $opts

  my \$import = VertRes::Pipelines::Import_iRODS_pacbio->new();
  \$import->convert_bax_to_fastq(\$bas_files);
  system('touch _convert_to_fastq_done');
  ];

    close($fh);
    VertRes::LSF::run(
        $lock_file, $work_dir, "${prefix}convert_to_fastq",
        { bsub_opts => "-M${memory_in_mb} -R 'select[mem>$memory_in_mb] rusage[mem=$memory_in_mb]' -n $threads -R 'span[hosts=1]' " },
        qq[perl -w ${prefix}convert_to_fastq.pl]
    );

    return $$self{No};
}

sub convert_bax_to_fastq {
    my ( $self, $bas_files ) = @_;

    for my $bas_file ( @{$bas_files} ) {
        my ( $filename, $dirs, $suffix ) = fileparse( $bas_file, '.bas.h5' );
        my $fastq = $filename . '.fastq.gz';
		my $subreads_bam = $filename . '.subreads.bam';
		my $scraps_bam = $filename . '.scraps.bam';
		
        #next if ( -e $fastq);
		if(-d 'bax2fastq')
		{
			remove_tree('bax2fastq');
		}
		
		system("pacbio_smrtpipe -t 8 -o bax2fastq bax2fastq $filename*.bax.h5");
		Utils::CMD(qq[mv bax2fastq/reads.fastq.gz $fastq]);
		
		if ( -e 'bax2fastq/subreads.bam')
		{
			Utils::CMD(qq[mv bax2fastq/subreads.bam $subreads_bam]);
			Utils::CMD(qq[mv bax2fastq/subreads.bam.pbi ${subreads_bam}.pbi]);
		}
		if ( -e 'bax2fastq/scraps.bam')
		{
			Utils::CMD(qq[mv bax2fastq/scraps.bam $scraps_bam]);
			Utils::CMD(qq[mv bax2fastq/scraps.bam.pbi ${scraps_bam}.pbi]);
		}
		if(-d 'bax2fastq')
		{
			remove_tree('bax2fastq');
		}

        my $fastqcheck = VertRes::Wrapper::fastqcheck->new();
        $fastqcheck->run( $fastq, $fastq . '.fastqcheck' );
    }
}

#---------- get_files ---------------------

# Requires nothing
sub get_files_requires {
    my ($self) = @_;
    return [];
}

sub get_files_provides {
    my ( $self, $lane_path ) = @_;
    return $self->_bas_h5_filenames;
}

sub get_files {
    my ( $self, $lane_path, $lock_file ) = @_;

    my $prefix   = $$self{prefix};
    my $work_dir = $lane_path;
    my $opts     = $self->dump_opts(qw(files));
    `mkdir -p $lane_path`;

    # Create a script to be run on LSF.
    open( my $fh, '>', "$work_dir/${prefix}import_files.pl" ) or $self->throw("$work_dir/${prefix}import_files.pl: $!");
    print $fh qq[
  use strict;
  use warnings;
  use VertRes::Pipelines::Import_iRODS_pacbio;

  my $opts

  my \$import = VertRes::Pipelines::Import_iRODS_pacbio->new(%\$opts);
  \$import->download_files();

  ];

    close($fh);
    VertRes::LSF::run( $lock_file, $work_dir, "${prefix}import_files", $self, qq[perl -w ${prefix}import_files.pl] );

    return $$self{No};
}

sub download_files {
    my ($self) = @_;

    my $irods = VertRes::Wrapper::iRODS->new();

    for my $ifile ( @{ $$self{files} } ) {
        next if ( $ifile =~ /\.bam/ );
        if ( !defined $ifile ) { $self->warn("No such file in iRODS? [$ifile]\n"); next; }
        my ( $filename, $dirs, $suffix ) = fileparse($ifile);
        next if ( -e $filename );

        if ( !( $ifile =~ m{([^/]+)$} ) ) { $self->throw("FIXME: [$ifile]"); }
        my $outfile = $1;
        if ( -e $outfile ) { next; }

        $irods->get_file( $ifile, "$outfile.tmp" );
        chmod 0664, "$outfile.bai";

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

    my $vrtrack = VRTrack::VRTrack->new( $$self{db} ) or $self->throw("Could not connect to the database\n");
    my $vrlane = VRTrack::Lane->new_by_name( $vrtrack, $$self{lane} ) or $self->throw("No such lane in the DB: [$$self{lane}]\n");

    $vrtrack->transaction_start();

    my $rawreads = 0;
    my $rawbases = 0;
    for my $fastq ( @{ $self->_output_fastq_filenames } ) {
        my $parser = VertRes::Parser::fastqcheck->new( file => $lane_path.'/'.$fastq . '.fastqcheck' );
        $rawreads += $parser->num_sequences() || 0;
        $rawbases += $parser->total_length()  || 0;
    }
	
	# Add the BAMs and FASTQ files to the file table	
    for my $file ((@{$self->_output_bam_filenames()},@{$self->_output_fastq_filenames()}))
    {
		my($filename_base, $dirs, $suffix) = fileparse($file);
		
		Utils::CMD(qq[md5sum $file > $file.md5]);
		my ($md5) = Utils::CMD(qq[awk '{printf "%s",\$1}' $file]);
		Utils::CMD(qq[rm -rf $file.md5]);
		
	    # The file may be absent from the database
	    my $vrfile = $vrlane->get_file_by_name($filename_base);
	    if ( !$vrfile ) 
	    { 
	        $vrfile = $vrlane->add_file($filename_base); 
	        $vrfile->hierarchy_name($filename_base);
	    }
		$vrfile->md5($md5);
        $vrfile->is_processed('import',1);
        $vrfile->update();
    }

	# Delete the h5 files
	for my $file (@{$self->_h5_filenames()})
	{
		if(defined($file) && -e $file )
		{
			unlink( $file );
		}
	}
	
	# Remove the h5 files from the file table
    for my $ifile ( @{ $$self{files} } ) {
        next if ( $ifile =~ /\.bam/ );
        my ( $filename, $dirs, $suffix ) = fileparse($ifile);

        if ( $ifile =~ m!^/seq/! ) { 
			my $vrfile = $vrlane->get_file_by_name($ifile);
	        $vrfile->is_latest(0);
	        $vrfile->update();
		}
	}
	
    $vrlane->is_processed( 'import', 1 );
    $vrlane->is_withdrawn(0);
    $vrlane->raw_reads($rawreads);
    $vrlane->raw_bases($rawbases);
    $vrlane->update();
    $vrtrack->transaction_commit();
	$self->update_file_permissions($lane_path);
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

