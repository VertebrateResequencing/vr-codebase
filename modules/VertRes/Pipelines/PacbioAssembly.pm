
=head1 NAME

VertRes::Pipelines::PacbioAssembly -  Pipeline for annotating an assembly

=head1 SYNOPSIS

# make the config files, which specifies the details for connecting to the
# VRTrack g1k-meta database and the data roots:
#echo '__VRTrack_PacbioAssembly__ annotation.conf' > pipeline.config
# where annotation.conf contains:
root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::PacbioAssembly',
prefix  => '_',
 
limit => 50,
db  => {
    database => 'g1k_meta',
    host     => 'mcs4a',
    port     => 3306,
    user     => 'vreseq_rw',
    password => 'xxxxxxx',
}

data => {

}

# __VRTrack_PacbioAssembly__ 

# run the pipeline:
run-pipeline -c pipeline.config

=head1 DESCRIPTION

Pipeline for assembling pacbio data
=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Pipelines::PacbioAssembly;

use strict;
use warnings;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::Library;
use VRTrack::Sample;
use VertRes::LSF;
use base qw(VertRes::Pipeline);
use VertRes::Utils::FileSystem;
use File::Spec;
use Utils;
use Bio::AssemblyImprovement::Circlator::Main;
use Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs;

our $actions = [
    {
        name     => 'pacbio_assembly',
        action   => \&pacbio_assembly,
        requires => \&pacbio_assembly_requires,
        provides => \&pacbio_assembly_provides
    },
    {
        name     => 'update_db',
        action   => \&update_db,
        requires => \&update_db_requires,
        provides => \&update_db_provides
    }
];

our %options = ( bsub_opts => '' ,
				 circularise => 0,
				 circlator_exec => '/software/pathogen/external/bin/circlator', # move to config files?
				 quiver_exec => '/software/pathogen/internal/prod/bin/pacbio_smrtanalysis');

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new( %options, actions => $actions, @args );
    if ( defined( $self->{db} ) ) {
        $self->{vrtrack} = VRTrack::VRTrack->new( $self->{db} ) or $self->throw("Could not connect to the database\n");
    }
    $self->{fsu} = VertRes::Utils::FileSystem->new;

    return $self;
}

sub pacbio_assembly_provides {
    my ($self) = @_;
    return [$self->{lane_path}."/pacbio_assembly/contigs.fa", $self->{lane_path}."/".$self->{prefix}."assembly_done"];
}

sub pacbio_assembly_requires {
    my ($self) = @_;
    my $file_regex = $self->{lane_path}."/".'*.bax.h5';
    my @files = glob( $file_regex);
    die 'no files to assemble' if(@files == 0);
    
    return \@files;
}

=head2 pacbio_assembly

 Title   : pacbio_assembly
 Usage   : $obj->pacbio('/path/to/lane', 'lock_filename');
 Function: Take an assembly from the assembly pipeline and automatically pacbio it.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub pacbio_assembly {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $memory_in_mb = $self->{memory} || 100000;
    my $threads = $self->{threads} || 16;
    my $genome_size_estimate = $self->{genome_size} || 4000000;
    my $files = join(' ', @{$self->pacbio_assembly_requires()});
    my $output_dir= $self->{lane_path}."/pacbio_assembly";
    my $queue = $self->{queue}|| "normal";
    my $pipeline_version = $self->{pipeline_version} || '7.0';
    my $target_coverage = $self->{target_coverage} || 30;
    my $umask    = $self->umask_str;
    my $lane_name = $self->{vrlane}->name;
    my $contigs_base_name = $self->generate_contig_base_name($lane_name);
    
    my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."assembly.pl");
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
  use strict;
  use Bio::AssemblyImprovement::Circlator::Main;
  use Bio::AssemblyImprovement::Quiver::Main;
  use Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs;
  $umask

  # ~~~~~ Basic HGAP assembly ~~~~~~
  system("rm -rf $output_dir");
  system("pacbio_assemble_smrtanalysis --no_bsub --target_coverage $target_coverage $genome_size_estimate $output_dir $files");
  die "No assembly produced\n" unless( -e qq[$output_dir/assembly.fasta]);  
  system("mv $output_dir/assembly.fasta $output_dir/contigs.fa");
  system("sed -i -e 's/|quiver\\\$//' $output_dir/contigs.fa"); # remove the |quiver from end of contig names
  system("mv $output_dir/All_output/data/aligned_reads.bam $output_dir/contigs.mapped.sorted.bam");
  system("mv $output_dir/All_output/data/corrected.fastq $self->{lane_path}/$lane_name.corrected.fastq");
  system("gzip -f -9 $self->{lane_path}/$lane_name.corrected.fastq"); 
  system("rm -rf $output_dir/All_output/");
  
  # ~~~~~~ Circlator ~~~~~~~
  if(defined($self->{circularise}) && $self->{circularise} == 1) {
        # rename contigs here so that circlator logs have final contig names
	Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs->new(
                                        input_assembly => qq[$output_dir/contigs.fa],
                                        base_contig_name => qq[$contigs_base_name])->run();
  
  	my \$circlator = Bio::AssemblyImprovement::Circlator::Main->new(
    			'assembly'	      => qq[$output_dir/contigs.fa],
    			'corrected_reads'     => qq[$self->{lane_path}].'/$lane_name.corrected.fastq.gz',
    			'circlator_exec'      => qq[$self->{circlator_exec}],
    			'working_directory'   => qq[$output_dir],
    			);
       \$circlator->run();
       my \$circlator_final_file = qq[$output_dir/circularised/circlator.final.fasta];
    
    # ~~~~~~ Quiver ~~~~~~~~~~
    if(-e \$circlator_final_file) {
	system("touch $self->{prefix}circularisation_done");
        my \$quiver = Bio::AssemblyImprovement::Quiver::Main->new(
    			'reference'           => \$circlator_final_file,
    			'bax_files'           => qq[$self->{lane_path}].'/*.bax.h5',
    			'working_directory'   => qq[$output_dir/circularised],	
			'quiver_exec'         => qq[$self->{quiver_exec}],
    			);
    	\$quiver->run();
    	
    	my \$quiver_final_file = qq[$output_dir/circularised/quiver/quiver.final.fasta];
    	my \$quiver_bam_file = qq[$output_dir/circularised/quiver/quiver.aligned_reads.bam];
    	my \$quiver_bai_file = qq[$output_dir/circularised/quiver/quiver.aligned_reads.bam.bai];
    	
    	if(-e \$quiver_final_file){
		system("touch $self->{prefix}quiver_done");
    		system(qq[mv \$quiver_final_file $output_dir/circularised/quiver/hgap.circlator.quiver.contigs.fa]); # rename final assembly
                system("sed -i -e 's/|quiver\\\$//' $output_dir/circularised/quiver/hgap.circlator.quiver.contigs.fa"); # remove the |quiver from ends of contig names
    		system(qq[mv $output_dir/contigs.fa $output_dir/hgap.contigs.fa]); # rename original hgap assembly
		system(qq[mv \$quiver_bam_file $output_dir/contigs.mapped.sorted.bam]); # copy over bam file
		system(qq[rm -f \$quiver_bai_file]); # clean this file up - we generate our own later. Should we use it?
		# create a symlink called contigs.fa pointing to the new quiverised fasta file so that *find scripts continue to work
    		system(qq[ln -s $output_dir/circularised/quiver/hgap.circlator.quiver.contigs.fa $output_dir/contigs.fa]);
    	}#end quiver success   
    }#end circlator success
  }#end if circularised
  
  # ~~~~~~~ Bamcheck, assemblystats, cleanup ~~~~~~~~~
  system("samtools index $output_dir/contigs.mapped.sorted.bam");
  system("bamcheck -c 1,20000,5 -r $output_dir/contigs.fa $output_dir/contigs.mapped.sorted.bam > $output_dir/contigs.mapped.sorted.bam.bc");
  system("assembly_stats $output_dir/contigs.fa  > $output_dir/contigs.fa.stats");
  system("rm -f $output_dir/contigs.mapped.sorted.bam"); # delete BAM and BAI files to save space
  system("rm -f $output_dir/contigs.mapped.sorted.bam.bai");
  
  if(not defined($self->{circularise}) || $self->{circularise} == 0) {
	# If circularisation has not been done, rename here (i.e. after bamcheck so that there are no problems with contig
	# names not matching the names in the BAM file generated by hgap)
  	  Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs->new(
  					input_assembly => qq[$output_dir/contigs.fa],
  					base_contig_name => qq[$contigs_base_name])->run();    		      
  }
  
  # touch done file and pipeline_version_7
  system("touch $output_dir/pipeline_version_$pipeline_version");
  system("touch $self->{prefix}assembly_done");
  exit;
  };
  close $scriptfh;
 
  my $job_name = $self->{prefix}.'pacbio_assembly';
      
    $self->archive_bsub_files($lane_path, $job_name);
    VertRes::LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-n$threads -q $queue -M${memory_in_mb} -R 'select[mem>$memory_in_mb] rusage[mem=$memory_in_mb] span[hosts=1]'"}, qq{perl -w $script_name});

    return $self->{No};
}

=head2 generate_contig_base_name

 Title   : generate_contig_base_name
 Usage   : 
 Function: Generates a name that can be used as a prefix for contig names in assembly (usually the lane's accession number + lane name)
 Returns : base name
 Args    : lane path

=cut

sub generate_contig_base_name
{
  my ($self, $lane_name) = @_;
  my $vrlane  = VRTrack::Lane->new_by_name($self->{vrtrack}, $lane_name) or $self->throw("No such lane in the DB: [".$lane_name."]");
  if(defined($vrlane->acc())) {
      return join('.',($vrlane->acc(),$lane_name));
    
  }
  return join('.',('',$lane_name));
}


=head2 update_db_requires

 Title   : update_db_requires
 Usage   : my $required_files = $obj->update_db_requires('/path/to/lane');
 Function: Find out what files the update_db action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub update_db_requires {
    my ($self, $lane_path) = @_;
    return $self->pacbio_assembly_provides();
}

=head2 update_db_provides

 Title   : update_db_provides
 Usage   : my $provided_files = $obj->update_db_provides('/path/to/lane');
 Function: Find out what files the update_db action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub update_db_provides {
   my ($self) = @_;
    return [ $self->{lane_path}."/".$self->{prefix}."pacbio_assembly_update_db_done"];
}

=head2 update_db

 Title   : update_db
 Usage   : $obj->update_db('/path/to/lane', 'lock_filename');
 Function: Records in the database that the lane has been improved.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut
sub update_db {
    my ($self, $lane_path, $action_lock) = @_;

    # all the done files are there, so just need to update the processed
    # flag in the database (if not already updated)
    my $vrlane = $self->{vrlane};
    return $$self{'Yes'} unless(defined($vrlane));
    
    my $vrtrack = $vrlane->vrtrack;
    
    return $$self{'Yes'} if $vrlane->is_processed('assembled');

    unless($vrlane->is_processed('assembled')){
      $vrtrack->transaction_start();
      $vrlane->is_processed('assembled',1);
      $vrlane->update() || $self->throw("Unable to set assembled status on lane $lane_path");
      $vrtrack->transaction_commit();
    }

    
    my $job_status =  File::Spec->catfile($lane_path, $self->{prefix} . 'job_status');
    Utils::CMD("rm $job_status") if (-e $job_status);
    
    
    my $prefix = $self->{prefix};
    # remove job files
    foreach my $file (qw(pacbio_assembly )) 
      {
        foreach my $suffix (qw(o e pl)) 
        {
          unlink($self->{fsu}->catfile($lane_path, $prefix.$file.'.'.$suffix));
        }
    }
    
    Utils::CMD("touch ".$self->{fsu}->catfile($lane_path,"$self->{prefix}pacbio_assembly_update_db_done")   );  
    $self->update_file_permissions($lane_path);
    return $$self{'Yes'};
}

1;

