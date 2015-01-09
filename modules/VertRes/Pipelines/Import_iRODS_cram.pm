
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
use Cwd 'abs_path';

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
    samtools_exec  => '/software/pathogen/external/apps/usr/bin/samtools-1.1.30'
);

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new( %options, actions => $actions, @args );
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

    my $fastqs ;
    $fastqs = $self->_list_of_fastq_files($lane_path);
    
    push(@{$fastqs},'_cram_to_fastq_done');

    return $fastqs;
}


sub _list_of_fastq_files {
    my ( $self, $lane_path ) = @_;

    my @fastqs ;
    for my $file (@{$self->{files}})
    {
      for my $f (@{$self->_fastq_from_cram($lane_path.'/'.$file,$self->{vrlane}->{is_paired})})
      {
        push(@fastqs,$f );
      }
    }
    return \@fastqs;
}


sub convert_to_fastq {
    my ( $self, $lane_path, $lock_file ) = @_;
    my $memory_in_mb = 3000;

    my $prefix   = $$self{prefix};
    my $work_dir = $lane_path;
    
    my $file = $lane_path.'/'.$self->{files}->[0];
    my $is_paired = $self->{vrlane}->{is_paired};
    
    open( my $fh, '>', "$work_dir/${prefix}convert_to_fastq.pl" ) or $self->throw("$work_dir/${prefix}convert_to_fastq.pl: $!");
    print $fh qq[
  use strict;
  use warnings;
  use VertRes::Pipelines::Import_iRODS_cram;
  
  my \$import = VertRes::Pipelines::Import_iRODS_cram->new();
  \$import->cram_to_fastq(qq[$file],$is_paired);
  system('touch _cram_to_fastq_done');
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
  my ( $self, $cram,$is_paired, ) = @_;
  my $full_cram_path = abs_path($cram);
  my ( $filename, $dirs, $suffix ) = fileparse( $full_cram_path, '.cram' );
  
  if($is_paired == 1)
  {
    return [$dirs.$filename.'_1.fastq.gz', $dirs.$filename.'_2.fastq.gz' ];
  }
  else
  {
    return [$dirs.$filename.'_1.fastq.gz' ];
  }
}


sub cram_to_fastq {
    my ( $self, $file,$is_paired ) = @_;

    my ( $filename, $dirs, $suffix ) = fileparse( $file, '.cram' );
    my $fastq_base = $dirs.$filename ;
    
    my @bamtofastq_command = ( 
      'bamtofastq',
      'collate=1',
      'inputformat=cram',
      "F=".$filename."_1.fastq.gz",
      "O=".$filename."_1.fastq_orphan.gz",
      "F2=".$filename."_2.fastq.gz",
      "O2=".$filename."_2.fastq_orphan.gz",
      "S=".$filename.".fastq.gz",
      'gz=1',
      'exclude=SECONDARY,SUPPLEMENTARY',
      '<',
      $file,
    );

    system(join(" ", @bamtofastq_command) );
    unless($is_paired)
    {
      system("mv ${filename}.fastq.gz  ${filename}_1.fastq.gz ");
    }
    unlink($filename."_1.fastq_orphan.gz") if(-e $filename."_1.fastq_orphan.gz");
    unlink($filename."_2.fastq_orphan.gz") if(-e $filename."_2.fastq_orphan.gz");
    unlink($filename.".fastq.gz")          if(-e $filename.".fastq.gz");
    
    my $fastqcheck = VertRes::Wrapper::fastqcheck->new();
    
    my @fastqcheck_files;
    for my $fastq (@{$self->_fastq_from_cram($file,$is_paired)})
    {
      push(@fastqcheck_files,$fastq . '.fastqcheck');
      $fastqcheck->run( $fastq, $fastq . '.fastqcheck' );
      Utils::CMD(qq[md5sum $fastq > $fastq.md5]);

      my $fastq_without_gz = $fastq;
      $fastq_without_gz  =~ s!\.gz!!i;
      Utils::CMD(qq[gunzip -c $fastq > $fastq_without_gz]);
      Utils::CMD(qq[md5sum $fastq_without_gz > $fastq_without_gz.md5]);
      unlink($fastq_without_gz);
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
    
    return $self->convert_to_fastq_provides($lane_path);
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
    my $vrtrack = VRTrack::VRTrack->new( $$self{db} ) or $self->throw("Could not connect to the database\n");
    $self->update_db_master($lane_path,$lock_file,$vrtrack);
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
    

   
    my $vrlane = VRTrack::Lane->new_by_name( $vrtrack, $$self{lane} ) or $self->throw("No such lane in the DB: [$$self{lane}]\n");

    $vrtrack->transaction_start();

    my $rawreads = 0;
    my $rawbases = 0;
    for my $fastq ( @{ $self->_list_of_fastq_files($lane_path) } ) {
        $self->throw("Cannot find ". $fastq . '.fastqcheck'."\n") unless(-e $fastq . '.fastqcheck');
        my $parser = VertRes::Parser::fastqcheck->new( file => $fastq . '.fastqcheck' );
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


sub update_db_master
{
    my ($self,$lane_path,$lock_file,$vrtrack) = @_;

    if ( !$$self{db} ) { $self->throw("Expected the db key.\n"); }

    my $vrlane  = VRTrack::Lane->new_by_name($vrtrack,$$self{lane}) or $self->throw("No such lane in the DB: [$$self{lane}]\n");

    $vrtrack->transaction_start();

    # To determine file types, a simple heuristic is used: When there are two files 
    #   lane_1.fastq.gz and lane_2.fastq.gz, fwd (1) and rev (2) will be set for file type.
    #   When only one file is present, single-end (0) will be set.
    my $i=0;
    my %processed_files;
    while (1)
    {
        $i++;
        my $name = "$$self{lane}_$i.fastq";

        # Check what fastq files actually exist in the hierarchy
        if ( ! -e "$lane_path/$name.gz" ) { last; }

        $processed_files{$name} = $i;
    }

    # do the same for every _nonhuman.fastq files
    $i=0;
    while (1)
    {
	$i++;
	my $name = "$$self{lane}_$i"."_nonhuman.fastq";

	# Check what fastq files actually exist in the hierarchy
	if ( ! -e "$lane_path/$name.gz" ) { last; }

	$processed_files{$name} = $i;
    }

    my $nfiles = scalar keys %processed_files;
    while (my ($name,$type) = each %processed_files)
    {
        # The file may be absent from the database, if it was created by splitting the _s_ fastq.
        my $vrfile = $vrlane->get_file_by_name($name);
        if ( !$vrfile ) 
        { 
            $vrfile = $vrlane->add_file($name); 
            $vrfile->hierarchy_name($name);
        }
        $vrfile->md5(`awk '{printf "%s",\$1}' $lane_path/$name.md5`);

        # Hm, this must be evaled, otherwise it dies without rollback
        my ($avg_len,$tot_len,$num_seq,$avg_qual);
        eval {
            my $fastq = VertRes::Parser::fastqcheck->new(file => "$lane_path/$name.gz.fastqcheck");
            $avg_len  = $fastq->avg_length();
            $tot_len  = $fastq->total_length();
            $num_seq  = $fastq->num_sequences();
            $avg_qual = $fastq->avg_qual();
        };
        if ( $@ )
        {
            $vrtrack->transaction_rollback();
            $self->throw("Problem reading the fastqcheck file: $lane_path/$name.gz.fastqcheck\n");
        }
        $vrfile->read_len($avg_len);
        $vrfile->raw_bases($tot_len);
        $vrfile->raw_reads($num_seq);
        $vrfile->mean_q($avg_qual); 
        $vrfile->is_processed('import',1);

        # Only the gzipped variant will carry the latest flag
        $vrfile->is_latest(0);  

        # If there is only one file, it is single-end file
        if ( $nfiles == 1 ) { $type = 0; }
        $vrfile->type($type);
        $vrfile->update();

        # Now add also the fastq.gz file into the File table. (Requested for the Mapping pipeline.)
        $self->fix_lane_file_if_exists("$name.gz",$vrlane, $vrtrack);
        my $vrfile_gz = $vrlane->get_file_by_name("$name.gz");
        if ( !$vrfile_gz )
        {
            $vrfile_gz = $vrlane->add_file("$name.gz");
            vrtrack_copy_fields($vrfile,$vrfile_gz,[qw(file_id name hierarchy_name)]);
            $vrfile_gz->hierarchy_name("$name.gz");
            $vrfile_gz->md5(`awk '{printf "%s",\$1}' $lane_path/$name.gz.md5`);
            $vrfile_gz->is_processed('import',1);
            $vrfile_gz->update();
        }
    }

    # Update also the status of the _s_ file, if any. Loop over the
    #   files passed to the pipeline and filter out those which were not
    #   processed above.
    for my $file (@{$$self{files}})
    {
        if ( exists($processed_files{$file}) ) { next; }

        my $vrfile = $vrlane->get_file_by_name($file);
        if ( !$vrfile ) { $self->throw("FIXME: no such file [$file] for the lane [$lane_path]."); }
        $vrfile->is_processed('import',1);
        $vrfile->is_latest(0);
        $vrfile->update();
    }


    # Update the lane stats
    my $read_len=0;
    my $raw_reads=0;
    my $raw_bases=0;
    my $vfiles=$vrlane->files;
    for my $vfile(@$vfiles){
	$raw_reads+=$vfile->raw_reads;
	$raw_bases+=$vfile->raw_bases;
	$read_len=$vfile->read_len;
    }
    $vrlane->raw_reads($raw_reads);
    $vrlane->raw_bases($raw_bases);
    $vrlane->read_len($read_len);


    # Finally, change the import status of the lane, so that it will not be picked up again
    #   by the run-pipeline script.
    $vrlane->is_processed('import',1);
    $vrlane->update();
    $vrtrack->transaction_commit();

    return $$self{Yes};
}

=head2 fix_lane_file_if_exists

        Example : fix_lane_file_if_exists("my_file_name", $vlane);
        Args    : filename to check
                : lane object which will in future have the file object

  This will see if a file row already exists for a lane. If its attached to a different lane, then delete it so that it doesnt cause
  issues later.
=cut

sub fix_lane_file_if_exists
{
  my ($self,$file_name, $vrlane, $vrtrack) = @_;
  
  if(VRTrack::File->is_name_in_database($vrtrack,$file_name,$file_name) == 1)
  {
    my $direct_file_object = VRTrack::File->new_by_name($vrtrack, $file_name);
    my $lane_file_object = $vrlane->get_file_by_name($file_name);
    
    if(! defined($lane_file_object)  || ( defined($lane_file_object) &&  $direct_file_object->row_id() != $lane_file_object->row_id()  ) )
    {
      # file exists already but its not attached to our lane, delete it
      $direct_file_object->delete();
    }
  }
}

sub vrtrack_copy_fields
{
    my ($src,$dst,$except) = @_;

    my %omit = $except ? map { $_=>1 } @$except : ();
    my $src_fields = $src->fields_dispatch();
    my $dst_fields = $dst->fields_dispatch();

    while (my ($key,$handler) = each %$src_fields)
    {
        if ( exists($omit{$key}) ) { next; }

        my $value = &{$handler}();
        &{$$dst_fields{$key}}($value);
    }
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

