
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
        name     => 'correct_reads',
        action   => \&correct_reads,
        requires => \&correct_reads_requires,
        provides => \&correct_reads_provides
    },
    {
        name     => 'hgap_assembly',
        action   => \&hgap_assembly,
        requires => \&hgap_assembly_requires,
        provides => \&hgap_assembly_provides
    },
    {
        name     => 'canu_assembly',
        action   => \&canu_assembly,
        requires => \&canu_assembly_requires,
        provides => \&canu_assembly_provides
    },
    {
        name     => 'modification',
        action   => \&modification,
        requires => \&modification_requires,
        provides => \&modification_provides
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
				 samtools_htslib_exec => '/software/pathogen/external/apps/usr/bin/samtools-1.6',
				 circlator_exec => '/software/pathogen/external/bin/circlator', # move to config files?
				 minimap2_exec => '/software/pathogen/external/apps/usr/bin/minimap2',
				 hgap_version_str => 'hgap_4_0',
				 canu_version_str => 'canu_1_6',
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

sub correct_reads_provides {
    my ($self) = @_;
    return [$self->{lane_path}.'/'.$self->{vrlane}->name.'.corrected.fasta.gz'];
}

sub correct_reads_requires {
    my ($self) = @_;
    return $self->get_subreads_bams();
}

# Generate corrected reads with CANU. HGAP corrected reads dont work with circlator
sub correct_reads {
    my ($self, $lane_path, $action_lock) = @_;

    my $queue = $self->{queue}|| "normal";
    my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."correct_reads.pl");

    my $umask    = $self->umask_str;
    my $lane_name = $self->{vrlane}->name;
	
    my @files = @{$self->get_subreads_bams()};
	my $uncorrected_fastq = $self->generate_fastq_filename_from_subreads_filename($files[0]);
	
	my $threads = 4;
	my $memory_in_gb = 10;
	my $memory_in_mb = $memory_in_gb*1000;
	# This can be hardcoded
	my $genome_size_estimate_kb = 8000;
	my $correction_output_dir = $self->{lane_path}.'/tmp_correction';
	
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
  use strict;
  $umask
  
  system("rm -rf $correction_output_dir");
  system("canu -correct -p corrected -d $correction_output_dir genomeSize=${genome_size_estimate_kb}k maxMemory=${memory_in_gb}g maxThreads=${threads} -pacbio-raw $uncorrected_fastq");
  system("mv $correction_output_dir/corrected.correctedReads.fasta.gz $self->{lane_path}/$lane_name.corrected.fasta.gz");
  system("rm -rf $correction_output_dir");
    
  exit;
  };
  close $scriptfh;
 
  my $job_name = $self->{prefix}.'correct_reads';
      
    $self->delete_bsub_files($lane_path, $job_name);
    VertRes::LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-n$threads -q $queue -M${memory_in_mb} -R 'select[mem>$memory_in_mb] rusage[mem=$memory_in_mb] span[hosts=1]'"}, qq{perl -w $script_name});

    return $self->{No};
}

sub canu_assembly_provides {
    my ($self) = @_;
     return [$self->{lane_path}."/".$self->{canu_version_str}."_assembly/contigs.fa", $self->{lane_path}."/".$self->{prefix}."canu_assembly_done"];
}

sub canu_assembly_requires {
    my ($self) = @_;
	my @files = @{$self->get_subreads_bams()};
	push(@files, $self->{lane_path}.'/'.$self->{vrlane}->name.'.corrected.fasta.gz');
	
    return \@files;
}

# Assemble with CANU
sub canu_assembly {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $memory_in_mb = 30000;
	my $memory_in_gb = $memory_in_mb/1000;
    my $threads = 8;
    my $genome_size_estimate = 8000000;
	my $genome_size_estimate_kb = $genome_size_estimate/1000;
	my $canu_output_dir = $self->{lane_path}.'/tmp_canu';

    my $files = join(' ', @{$self->get_subreads_bams()});
	
    my $output_dir= $self->{lane_path}.'/canu_1_6_assembly';
	my $resequence_output_dir = $output_dir."/resequence";
    my $queue = $self->{queue}|| "normal";
    my $pipeline_version = '9.0';
    my $target_coverage = $self->{target_coverage} || 25;
    my $umask    = $self->umask_str;
    my $lane_name = $self->{vrlane}->name;
    my $contigs_base_name = $self->generate_contig_base_name($lane_name);
	my $correction_output_dir = $output_dir.'/correction';
	my $corrected_reads = $self->{lane_path}."/$lane_name.corrected.fasta.gz";
	
    my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."canu_assembly.pl");
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
  use strict;
  use Bio::AssemblyImprovement::Circlator::Main;
  use Bio::AssemblyImprovement::Quiver::Main;
  use Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs;
  $umask

  # ~~~~~ Basic HGAP assembly ~~~~~~
  system("rm -rf $canu_output_dir");
  system("canu -p canu -d $canu_output_dir genomeSize=${genome_size_estimate_kb}k maxMemory=${memory_in_gb}g maxThreads=${threads} -pacbio-corrected $corrected_reads");
  system("mv $canu_output_dir/canu.unitigs.fasta $self->{lane_path}/$output_dir/contigs.fa");
  system("rm -rf $canu_output_dir");
  
  # ~~~~~~ Circlator ~~~~~~~
  if(defined($self->{circularise}) && $self->{circularise} == 1) {
        # rename contigs here so that circlator logs have final contig names
	Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs->new(
                                        input_assembly => qq[$output_dir/contigs.fa],
                                        base_contig_name => qq[$contigs_base_name])->run();
  
       system(qq[$self->{circlator_exec} --assembler canu $output_dir/contigs.fa $corrected_reads $output_dir/circularised]);
       my \$circlator_final_file = qq[$output_dir/circularised/06.fixstart.fasta];
    
	# ~~~~~~ Quiver/Resequencing ~~~~~~~~~~
	if(-e \$circlator_final_file) {
		system("rm -f $resequence_output_dir");
		system("pacbio_smrtpipe -t $threads -r \$circlator_final_file -o $resequence_output_dir resequence $files");
		system("sed -i -e 's/|quiver\\\$//' $resequence_output_dir/contigs.fa");
		system(qq[mv $output_dir/contigs.fa $output_dir/canu.contigs.fa]);
		system(qq[mv $resequence_output_dir/contigs.fa $output_dir/contigs.fa]);
		system("rm -rf $resequence_output_dir");
	}
  }#end if circularised
  
  # map corrected reads to assembly
  system("$self->{minimap2_exec} -ax map-pb -t $threads $output_dir/contigs.fa $lane_name.corrected.fasta.gz | $self->{samtools_htslib_exec} sort -o $output_dir/contigs.mapped.sorted.bam -");
  system("$self->{samtools_htslib_exec} index $output_dir/contigs.mapped.sorted.bam");
  system("$self->{samtools_htslib_exec} stats -r $output_dir/contigs.fa $output_dir/contigs.mapped.sorted.bam > $output_dir/contigs.mapped.sorted.bam.bc");

  system("assembly-stats $output_dir/contigs.fa  > $output_dir/contigs.fa.stats");
  system("rm -f $output_dir/contigs.mapped.sorted.bam"); # delete BAM and BAI files to save space
  system("rm -f $output_dir/contigs.mapped.sorted.bam.bai");
  
  if(not defined($self->{circularise}) || $self->{circularise} == 0) {
	# If circularisation has not been done, rename here
  	  Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs->new(
  					input_assembly => qq[$output_dir/contigs.fa],
  					base_contig_name => qq[$contigs_base_name])->run();    		      
  }
  
  # touch done file and pipeline_version_9
  system("touch $output_dir/pipeline_version_$pipeline_version");
  system("touch $self->{prefix}canu_assembly_done");
  exit;
  };
  close $scriptfh;
 
  my $job_name = $self->{prefix}.'canu_assembly';
      
    $self->delete_bsub_files($lane_path, $job_name);
    VertRes::LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-n$threads -q $queue -M${memory_in_mb} -R 'select[mem>$memory_in_mb] rusage[mem=$memory_in_mb] span[hosts=1]'"}, qq{perl -w $script_name});

    return $self->{No};
}


sub hgap_assembly_provides {
    my ($self) = @_;
    return [$self->{lane_path}."/".$self->{hgap_version_str}."_assembly/contigs.fa", $self->{lane_path}."/".$self->{prefix}."hgap_assembly_done"];
}

sub hgap_assembly_requires {
    my ($self) = @_;
	my @files = @{$self->get_subreads_bams()};
	push(@files, $self->correct_reads_provides()->[0]);
    return \@files;
}

# TODO: refactor to remove all the logic from the script, split into multiple tasks
sub hgap_assembly {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $memory_in_mb = 40000;
	my $memory_in_gb = $memory_in_mb/1000;
    my $threads = 12;
    my $genome_size_estimate = $self->{genome_size} || 8000000;
	my $genome_size_estimate_kb = $genome_size_estimate/1000;

    my $files = join(' ', @{$self->get_subreads_bams()});
	
    my $output_dir= $self->{lane_path}."/".$self->{hgap_version_str}."_assembly";
	my $resequence_output_dir = $output_dir."/resequence";
	my $modification_output_dir = $output_dir."/modification";
    my $queue = $self->{queue}|| "normal";
    my $pipeline_version = '8.0';
    my $target_coverage = $self->{target_coverage} || 25;
    my $umask    = $self->umask_str;
    my $lane_name = $self->{vrlane}->name;
    my $contigs_base_name = $self->generate_contig_base_name($lane_name);
	my $correction_output_dir = $output_dir.'/correction';
	my $uncorrected_fastq = $self->generate_fastq_filename_from_subreads_filename($self->hgap_assembly_requires()->[0]);
	my $corrected_reads = $self->{lane_path}."/$lane_name.corrected.fasta.gz";
	
    my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."hgap_assembly.pl");
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
  use strict;
  use Bio::AssemblyImprovement::Circlator::Main;
  use Bio::AssemblyImprovement::Quiver::Main;
  use Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs;
  $umask

  # ~~~~~ Basic HGAP assembly ~~~~~~
  system("rm -rf $output_dir");
  system("pacbio_smrtpipe -t $threads --genome_size $genome_size_estimate --min_coverage $target_coverage -o $output_dir assembly $files");
  
  die "No assembly produced\n" unless( -e qq[$output_dir/contigs.fa]);  
  system("sed -i -e 's/|quiver\\\$//' $output_dir/contigs.fa"); # remove the |quiver from end of contig names
  system("mv $correction_output_dir/corrected.fastq.gz $self->{lane_path}/$lane_name.hgap_corrected.fastq.gz");
  
  # ~~~~~~ Circlator ~~~~~~~
  if(defined($self->{circularise}) && $self->{circularise} == 1) {
        # rename contigs here so that circlator logs have final contig names
	Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs->new(
                                        input_assembly => qq[$output_dir/contigs.fa],
                                        base_contig_name => qq[$contigs_base_name])->run();
  
       system(qq[$self->{circlator_exec} --assembler canu $output_dir/contigs.fa $corrected_reads $output_dir/circularised]);
       my \$circlator_final_file = qq[$output_dir/circularised/06.fixstart.fasta];
    
	# ~~~~~~ Quiver/Resequencing ~~~~~~~~~~
	if(-e \$circlator_final_file) {
		system("rm -f $resequence_output_dir");
		system("pacbio_smrtpipe -t $threads -r \$circlator_final_file -o $resequence_output_dir resequence $files");
		system("sed -i -e 's/|quiver\\\$//' $resequence_output_dir/contigs.fa");
		system(qq[mv $output_dir/contigs.fa $output_dir/hgap.contigs.fa]);
		system(qq[mv $resequence_output_dir/contigs.fa $output_dir/contigs.fa]);
		system("rm -rf $resequence_output_dir");
	}
  }#end if circularised
  
  # map corrected reads to assembly
  system("$self->{minimap2_exec} -ax map-pb -t $threads $output_dir/contigs.fa $lane_name.corrected.fasta.gz | $self->{samtools_htslib_exec} sort -o $output_dir/contigs.mapped.sorted.bam -");
  system("$self->{samtools_htslib_exec} index $output_dir/contigs.mapped.sorted.bam");
  system("$self->{samtools_htslib_exec} stats -r $output_dir/contigs.fa $output_dir/contigs.mapped.sorted.bam > $output_dir/contigs.mapped.sorted.bam.bc");

  system("assembly-stats $output_dir/contigs.fa  > $output_dir/contigs.fa.stats");
  system("rm -f $output_dir/contigs.mapped.sorted.bam"); # delete BAM and BAI files to save space
  system("rm -f $output_dir/contigs.mapped.sorted.bam.bai");
  
  if(not defined($self->{circularise}) || $self->{circularise} == 0) {
	# If circularisation has not been done, rename here (i.e. after bamcheck so that there are no problems with contig
	# names not matching the names in the BAM file generated by hgap)
  	  Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs->new(
  					input_assembly => qq[$output_dir/contigs.fa],
  					base_contig_name => qq[$contigs_base_name])->run();    		      
  }
  
  # touch done file and pipeline_version_8
  system("touch $output_dir/pipeline_version_$pipeline_version");
  system("touch $self->{prefix}hgap_assembly_done");
  exit;
  };
  close $scriptfh;
 
  my $job_name = $self->{prefix}.'hgap_assembly';
      
    $self->delete_bsub_files($lane_path, $job_name);
    VertRes::LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-n$threads -q $queue -M${memory_in_mb} -R 'select[mem>$memory_in_mb] rusage[mem=$memory_in_mb] span[hosts=1]'"}, qq{perl -w $script_name});

    return $self->{No};
}

sub modification_provides {
    my ($self) = @_;
    return [$self->{lane_path}."/".$self->{hgap_version_str}."_assembly/motifs.gff"];
}

sub modification_requires {
    my ($self) = @_;
    my @files = @{$self->get_subreads_bams()};
	push(@files, "$self->{prefix}hgap_assembly_done");
	
	my $assembly_file = $self->{lane_path}."/".$self->{hgap_version_str}."_assembly/contigs.fa";
	die 'assembly file empty' if( -z $assembly_file);
	push(@files, $assembly_file);
    
    return \@files;
}

# Find modified bases (Methylation)
sub modification {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $memory_in_mb = 10000;
    my $threads = 8;
	my $files = join(' ', @{$self->get_subreads_bams()});
    my $output_dir= $self->{lane_path}."/".$self->{hgap_version_str}."_assembly";
	my $modification_output_dir = $output_dir."/tmp_modification";
    my $queue = $self->{queue}|| "normal";
    my $umask    = $self->umask_str;
    my $lane_name = $self->{vrlane}->name;

    my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."modification.pl");

    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
  use strict;
  $umask

  system("rm -rf $modification_output_dir");
  system("pacbio_smrtpipe -t $threads -r $output_dir/contigs.fa -o $modification_output_dir modification $files");
  system(qq[mv $modification_output_dir/motifs.gff $output_dir/motifs.gff]);
  system("rm -rf $modification_output_dir");
 
  exit;
  };
  close $scriptfh;
 
  my $job_name = $self->{prefix}.'modification';
      
    $self->delete_bsub_files($lane_path, $job_name);
    VertRes::LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-n$threads -q $queue -M${memory_in_mb} -R 'select[mem>$memory_in_mb] rusage[mem=$memory_in_mb] span[hosts=1]'"}, qq{perl -w $script_name});

    return $self->{No};
}


sub generate_fastq_filename_from_subreads_filename
{
	my ($self, $filename) = @_;
	$filename =~ s!.subreads.bam!.fastq.gz!;
	return $filename;
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


sub get_subreads_bams
{
	my ($self) = @_;
	my $file_regex = $self->{lane_path}."/".'*.subreads.bam';
	my @files = glob( $file_regex);
	die 'no input BAM files' if(@files == 0);
	return \@files;
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
	my @all_assemblies = (@{$self->hgap_assembly_provides()},@{$self->canu_assembly_provides()},  @{$self->modification_provides()});
    return \@all_assemblies;
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
    foreach my $file (qw( hgap_assembly canu_assembly correct_reads modification)) 
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

