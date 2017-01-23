=head1 NAME

VertRes::Pipelines::Assembly - Assemble genomes

=head1 SYNOPSIS

echo '__Assembly__ assembly.conf' > pipeline.config
# where assembly.conf contains:
root    => '/lustre/scratch118/infgen/pathogen/pathpipe/prokaryotes/seq-pipelines',
module  => 'VertRes::Pipelines::Assembly',
prefix  => '_',
limit   => 30,

limits => {
 project => ['ABC']
},
vrtrack_processed_flags => {  stored => 1, import => 1, assembled => 0 },

db  => {
    database => 'pathogen_prok_track',
    host     => 'xxx',
    port     => yyyy,
    user     => 'pathpipe_rw',
    password => 'xxx',
},

data => {
    db  => {
        database => 'pathogen_prok_track',
        host     => 'xxx',
        port     => yyyy,
        user     => 'pathpipe_rw',
        password => 'xxx',
    },
    # rough estimate so we know how much RAM to request
    genome_size => 10000000,

    seq_pipeline_root    => '/lustre/scratch118/infgen/pathogen/pathpipe/prokaryotes/seq-pipelines',
    no_scaffolding => 0,
    improve_assembly => 1, # To run abacas (if reference given), scaffolding and gap filling
    annotation     => 1,
    error_correct  => 1, # Should the reads be put through an error correction stage?
    normalise	   => 1, # Should we do digital normalisation?
    remove_primers => 1,
    primers_file   => /path/to/primers/file, #Essential if primers are to be removed
    primer_removal_tool => 'quasr' # quasr is the only option for now.

    assembler => 'velvet',
    assembler_exec => '/software/pathogen/external/apps/usr/bin/velvet',
    optimiser_exec => '/software/pathogen/external/apps/usr/bin/VelvetOptimiser.pl',
    sga_exec       => '/software/pathogen/external/apps/usr/local/src/SGA/sga',
    khmer_exec		=> '/software/pathogen/external/apps/usr/local/khmer/scripts/normalize-by-median.py',
    QUASR_exec		=> '/software/pathogen/internal/pathdev/java/QUASR702_Parser/readsetProcessor.jar',
    trimmomatic_jar => '/software/pathogen/external/apps/usr/local/Trimmomatic-0.32/trimmomatic-0.32.jar',
    remove_adapters => 1, # should we remove adapters?
    adapters_file   => '/lustre/scratch118/infgen/pathogen/pathpipe/usr/share/solexa-adapters.fasta',
    adapter_removal_tool => 'trimmomatic' # trimmomatic is the only option for now
    max_threads => 1,
    single_cell => 1, # Put this in to assemble single cell data. For normal assemblies, leave it out.
    iva_qc => 1, # If set, run iva_qc. Default - do not run iva_qc
    kraken_db => 'path to kraken db', #for iva_qc
},


=head1 DESCRIPTION

This pipeline requires velvet, prokka, smalt, SSPACE, SGA, khmer, QUASR and GapFiller to be separately installed from the original authors software repositories.


=head1 AUTHOR

path-help@sanger.ac.uk

=cut
package VertRes::Pipelines::Assembly;

use strict;
use warnings;
use VertRes::Utils::Hierarchy;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VertRes::Parser::fastqcheck;
use VRTrack::Library;
use VRTrack::File;
use File::Basename;
use Time::Format;
use File::Copy;
use Cwd;
use VertRes::LSF;
use Data::Dumper;
use FileHandle;
use Utils;
use List::Util qw(min max);
use File::Path qw(make_path remove_tree);
use File::Touch;
use VertRes::Utils::Assembly;
use VertRes::Utils::Sam;
use VertRes::Wrapper::smalt;
use Bio::AssemblyImprovement::Scaffold::Descaffold;
use Bio::AssemblyImprovement::Scaffold::SSpace::PreprocessInputFiles;
use Bio::AssemblyImprovement::Scaffold::SSpace::Iterative;
use Bio::AssemblyImprovement::FillGaps::GapFiller::Iterative;
use Bio::AssemblyImprovement::Assemble::SGA::Main;
use Bio::AssemblyImprovement::DigitalNormalisation::Khmer::Main;
use Bio::AssemblyImprovement::PrimerRemoval::Main;
use Bio::AssemblyImprovement::Util::FastqTools;
use Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs;
use Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs;
use Bio::AssemblyImprovement::Util::FastaTools;
use Bio::AssemblyImprovement::Util::OrderContigsByLength;
use Bio::AssemblyImprovement::IvaQC::Main;


use base qw(VertRes::Pipeline);

our $actions = [ { name     => 'pool_fastqs',
                   action   => \&pool_fastqs,
                   requires => \&pool_fastqs_requires,
                   provides => \&pool_fastqs_provides },
                 { name     => 'optimise_parameters',
                   action   => \&optimise_parameters,
                   requires => \&optimise_parameters_requires,
                   provides => \&optimise_parameters_provides },
                 { name     => 'assembly_improvement',
                   action   => \&assembly_improvement,
                   requires => \&assembly_improvement_requires,
                   provides => \&assembly_improvement_provides },
                 { name     => 'iva_qc',
                   action   => \&iva_qc,
                   requires => \&iva_qc_requires,
                   provides => \&iva_qc_provides },
                 { name     => 'map_back',
                   action   => \&map_back,
                   requires => \&map_back_requires,
                   provides => \&map_back_provides },
                 { name     => 'update_db',
                   action   => \&update_db,
                   requires => \&update_db_requires,
                   provides => \&update_db_provides },
                 { name     => 'cleanup',
                   action   => \&cleanup,
                   requires => \&cleanup_requires,
                   provides => \&cleanup_provides }
                ];

our %options = (
                do_cleanup       => 1,
                scaffolder_exec  => '/software/pathogen/external/apps/usr/local/SSPACE-BASIC-2.0_linux-x86_64/SSPACE_Basic_v2.0.pl',
                gap_filler_exec  => '/software/pathogen/external/apps/usr/local/GapFiller_v1-11_linux-x86_64/GapFiller.pl',
                abacas_exec      => '/software/pathogen/internal/prod/bin/abacas.pl',
                sga_exec         => '/software/pathogen/external/apps/usr/bin/sga',
                khmer_exec		 => '/software/pathogen/external/apps/usr/local/khmer/scripts/normalize-by-median.py',
                QUASR_exec	     => '/software/pathogen/external/apps/usr/local/QUASR/readsetProcessor.jar',
                trimmomatic_jar  => '/software/pathogen/external/apps/usr/local/Trimmomatic-0.32/trimmomatic-0.32.jar',
                iva_qc_exec		 => '/software/pathogen/external/bin/iva_qc',
                iva_qc	=> 0,
                kraken_db => '/lustre/scratch118/infgen/pathogen/pathpipe/kraken/assemblyqc_fluhiv_20150728',
                adapters_file    => '/lustre/scratch118/infgen/pathogen/pathpipe/usr/share/solexa-adapters.fasta',
                primers_file     => '',
                remove_primers   => 0,
                remove_adapters  => 0,
                primer_removal_tool  => 'quasr',
                adapter_removal_tool => 'trimmomatic',
                no_scaffolding   => 0,
                annotation       => 0,
                single_cell      => 0,
                improve_assembly => 1,
                iva_seed_ext_min_cov   => 5,
                iva_seed_ext_min_ratio => 2,
                iva_ext_min_cov        => 5,
                iva_ext_min_ratio      => 2,
                iva_insert_size		   => 800,
                iva_strand_bias        => 0,

               );

sub new {
  my ($class, @args) = @_;

  my $self = $class->SUPER::new(%options, actions => $actions, @args);
  if(defined($self->{db}))
  {
    $self->{vrtrack} = VRTrack::VRTrack->new($self->{db}) or $self->throw("Could not connect to the database\n");
  }
  $self->{fsu} = VertRes::Utils::FileSystem->new;

  unless(defined($self->{max_threads}))
  {
    $self->{max_threads} = 1;
  }

  if(defined($self->{vrlane}))
  {
    # being run on individual lanes
    $self->{pools} =
       [{
            lanes => [$self->{vrlane}->name()],
            type => 'shortPaired',
       }];
  }

  my $assembly_util = VertRes::Utils::Assembly->new(assembler => $self->{assembler});
  my $assembler_class = $assembly_util->find_module();
  eval "require $assembler_class;";
  $self->{assembler_class} = $assembler_class;

  $self->{pipeline_version} = 2 unless(defined($self->{pipeline_version}) );

  return $self;
}

###########################
# Begin map_back
###########################

sub map_back_requires
{
  my ($self) = @_;
  if (exists $self->{iva_qc} and $self->{iva_qc} and defined(qq[$self->{kraken_db}])) {
    return $self->iva_qc_provides();
  }
  else {
    return $self->optimise_parameters_provides();
  }
}

sub map_back_provides
{
   my ($self) = @_;

   return [$self->{lane_path}."/".$self->{prefix}.$self->{assembler}.'_map_back_done'];

}

sub map_back
{
  my ($self, $build_path, $action_lock) = @_;
  $self->mapping_and_generate_stats($build_path, $action_lock, "optimised_directory", "map_back",$self->{reference});

  return $self->{No};
}


###########################
# End map_back
###########################



sub mapping_and_generate_stats
{
  my ($self, $build_path, $action_lock, $working_directory_method_name, $action_name_suffix, $reference) = @_;

  my $assembler_class = $self->{assembler_class};
  my $output_directory = $self->{lane_path};
  eval("use $assembler_class; ");
  my $assembler_util= $assembler_class->new( output_directory => qq[$output_directory]);
  my $base_path = $self->{seq_pipeline_root};
  my $umask    = $self->umask_str;
  my $job_name = $self->{prefix}.$self->{assembler}."_$action_name_suffix";
  my $script_name = $self->{fsu}->catfile($output_directory, $self->{prefix}.$self->{assembler}."_$action_name_suffix.pl");

  my $lane_names = $self->get_all_lane_names($self->{pools});
  my @lane_paths;
  for my $lane_name (@$lane_names)
  {
    push(@lane_paths,$base_path.'/'.$self->{vrtrack}->hierarchy_path_of_lane_name($lane_name).'/'.$lane_name);
  }
  my $lane_paths_str = '("'.join('","', @lane_paths).'")';
  my $reference_str = '';
  if($reference)
  {
    $reference_str = qq{,reference => qq[$reference]} ;
  }
  my $annotation = $self->{annotation};

  open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
  print $scriptfh qq{
  use strict;
  use $assembler_class;
  $umask
  my \@lane_paths = $lane_paths_str;

  my \$assembler_util= $assembler_class->new( output_directory => qq[$output_directory] $reference_str );
  my \$directory = \$assembler_util->${working_directory_method_name}();
  \$assembler_util->map_and_generate_stats(\$directory,\$directory, \\\@lane_paths );

  unlink("\$directory/contigs.mapped.sorted.bam.bai");
  unlink("\$directory/contigs.mapped.sorted.bam");
  unlink("\$directory/contigs.fa.small.sma");
  unlink("\$directory/contigs.fa.small.smi");

  system("touch \$directory/$self->{prefix}$self->{assembler}_${action_name_suffix}_done");

  system("touch $self->{prefix}$self->{assembler}_${action_name_suffix}_done");
  exit;
                };
  close $scriptfh;

  my $memory_required_mb = 1700;

  VertRes::LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => " -M${memory_required_mb} -R 'select[mem>$memory_required_mb] rusage[mem=$memory_required_mb]'"}, qq{perl -w $script_name});

}


###########################
# Begin optimise
###########################

sub optimise_parameters_provides
{
  my $self = shift;

  return  [$self->{lane_path}."/".$self->{prefix}."$self->{assembler}_optimise_parameters_done"];
}

sub optimise_parameters_requires
{
  my $self = shift;
  return $self->pool_fastqs_provides();
}

sub optimise_parameters
{
    my ($self, $build_path,$action_lock) = @_;

    my $lane_names = $self->get_all_lane_names($self->{pools});
    my $output_directory = $self->{lane_path};
    my $pool_directory = $self->{lane_path}."/".$self->{prefix}.$self->{assembler}."_pool_fastq_tmp_files";
    my $assembler_class = $self->{assembler_class};
    eval("use $assembler_class; ");
    my $assembler_util = $assembler_class->new();
    my $files_str = $assembler_util->generate_files_str($self->{pools}, $pool_directory);
    my $job_name = $self->{prefix}.$self->{assembler}.'_optimise_parameters';
    my $script_name = $self->{fsu}->catfile($self->{lane_path}, $self->{prefix}.$self->{assembler}."_optimise_parameters.pl");
    my $lane_paths_str = $self->get_lane_paths_str();
    my $umask_str = $self->umask_str;
    my ($memory_required_mb, $num_threads) = $self->get_memory_and_threads();
    my $insert_size = $self->get_insert_size();
    my $tmp_directory = $self->{tmp_directory}.'/'.$self->{prefix}.$self->{assembler}.'_'.$lane_names->[0] || getcwd();
    my $pipeline_version = join('/',($output_directory, $self->{assembler}.'_assembly','pipeline_version_'.$self->{pipeline_version}));
    my $contigs_base_name = $self->generate_contig_base_name();

    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");

    (my $script_string = qq{
        use strict;
        use $assembler_class;
        use VertRes::Pipelines::Assembly;
        use File::Copy;
        use Cwd;
        use File::Path qw(make_path remove_tree);
        use File::Touch;
        use Bio::AssemblyImprovement::Util::FastqTools;
        use Bio::AssemblyImprovement::IvaQC::Main;
        $umask_str

        my \$assembly_pipeline = VertRes::Pipelines::Assembly->new(
          assembler => "$self->{assembler}"
        );
        system("rm -rf $self->{assembler}_assembly_*");

        remove_tree(qq[$tmp_directory]) if(-d qq[$tmp_directory]);
        make_path(qq[$tmp_directory]) or die "Error mkdir $tmp_directory";
        chdir(qq[$tmp_directory]) or die "Error chdir $tmp_directory";

        # Calculate 66-90% of the median read length as min and max kmer values
        my \$fastq_tools  = Bio::AssemblyImprovement::Util::FastqTools->new(
            input_filename   => "$pool_directory/pool_1.fastq.gz",
            single_cell => $self->{single_cell},
        );

        my \%kmer;
        if ('$self->{assembler}' eq 'iva')
        {
          \%kmer = (min => 0, max => 0);
        }
        else
        {
          \%kmer = \$fastq_tools->calculate_kmer_sizes();
        }

        my \$assembler = $assembler_class->new(
          assembler => qq[$self->{assembler}],
          optimiser_exec => qq[$self->{optimiser_exec}],
          min_kmer => \$kmer{min},
          max_kmer => \$kmer{max},
          files_str => qq[$files_str],
          output_directory => qq[$tmp_directory],
          single_cell => $self->{single_cell},
          spades_kmer_opts => qq[$self->{spades_kmer_opts}],
          spades_opts => qq[$self->{spades_opts}],
          trimmomatic_jar => qq[$self->{trimmomatic_jar}],
          remove_primers => $self->{remove_primers},
          primer_removal_tool => qq[$self->{primer_removal_tool}],
          primers_file => qq[$self->{primers_file}],
          remove_adapters => $self->{remove_adapters},
          adapter_removal_tool => qq[$self->{adapter_removal_tool}],
          adapters_file => qq[$self->{adapters_file}],
          iva_seed_ext_min_cov   => qq[$self->{iva_seed_ext_min_cov}],
          iva_seed_ext_min_ratio => qq[$self->{iva_seed_ext_min_ratio}],
          iva_ext_min_cov        => qq[$self->{iva_ext_min_cov}],
          iva_ext_min_ratio      => qq[$self->{iva_ext_min_ratio}],
          iva_insert_size	     => qq[$self->{iva_insert_size}],
          iva_strand_bias		 => qq[$self->{iva_strand_bias}],
        );

        my \$ok = \$assembler->optimise_parameters($num_threads);

        # this is to stop the map back stage crashing because smalt needs all
        # contigs to be at least a kmer long
        my \$fasta_processor = Bio::AssemblyImprovement::Util::FastaTools->new(input_filename => \$assembler->optimised_assembly_file_path());
        \$fasta_processor->remove_small_contigs(50, 0);

        Bio::AssemblyImprovement::Util::OrderContigsByLength->new(
            input_filename => \$fasta_processor->output_filename,
            output_filename => \$assembler->optimised_assembly_file_path()
        )->run();

        copy(\$assembler->optimised_assembly_file_path(),\$assembler->optimised_directory().'/unscaffolded_contigs.fa');

        if('$self->{assembler}' eq 'velvet')
        {
          die "Missing assembly logfile - velvet optimiser didnt complete"  unless(-e qq[$tmp_directory].'/$self->{assembler}_assembly_logfile.txt');
          move(qq[$tmp_directory].'/$self->{assembler}_assembly_logfile.txt', qq[$output_directory].'/$self->{assembler}_assembly_logfile.txt');
        }

        die "Missing assembly directory - assembly didnt complete" unless(-d "$tmp_directory/$self->{assembler}_assembly");
        system("mv $tmp_directory/$self->{assembler}_assembly $output_directory") and die "Error mv $tmp_directory/$self->{assembler}_assembly $output_directory";
        chdir(qq[$output_directory]);
        remove_tree(qq[$tmp_directory]) or die "Error remove_tree $tmp_directory";
        #unlink(qq[$tmp_directory].'/contigs.fa.scaffolded.filtered');
        touch(qq[$pipeline_version]) or die "Error touch $pipeline_version";
        my \$done_file = '$output_directory/$self->{prefix}$self->{assembler}_optimise_parameters_done';
        touch(\$done_file) or die "Error touch \$done_file";
    }) =~ s/^ {8}//mg;

    print $scriptfh $script_string;
    close $scriptfh;

    my $total_memory_mb;
    if ($self->{assembler} eq 'iva') {
        $total_memory_mb = 2000;
    }
    else {
        $total_memory_mb = $num_threads*$memory_required_mb < 5000 ? 5000 : $num_threads*$memory_required_mb;
    }

    my $queue = $self->decide_appropriate_queue($memory_required_mb);

    VertRes::LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => "-n$num_threads -q $queue -M${total_memory_mb} -R 'select[mem>$total_memory_mb] rusage[mem=$total_memory_mb] span[hosts=1]'", dont_wait=>1}, qq{perl -w $script_name});

    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
    return $self->{No};
}


sub assembly_improvement_provides
{
  my $self = shift;

  return  [$self->{lane_path}."/".$self->{prefix}."$self->{assembler}_assembly_improvement_done"];
}


sub assembly_improvement_requires
{
  my $self = shift;
  return $self->optimise_parameters_provides();
}


sub assembly_improvement
{
    my ($self, $build_path,$action_lock) = @_;

    my $assembler_class = $self->{assembler_class};
    eval("use $assembler_class; ");
    my $insert_size = $self->get_insert_size();
    my $output_directory = $self->{lane_path};
    my $assembly_directory = $self->{fsu}->catfile($output_directory, "$self->{assembler}_assembly");
    my $fasta_for_improvement = $self->{fsu}->catfile($assembly_directory, "contigs.fa");
    my $lane_names = $self->get_all_lane_names($self->{pools});
    my $tmp_directory = $self->{tmp_directory}.'/'.$self->{prefix}.$self->{assembler}.'_'.$lane_names->[0] || getcwd();
    my $pool_directory = $self->{lane_path}."/".$self->{prefix}.$self->{assembler}."_pool_fastq_tmp_files";
    my ($temp_reads_dir, $reads_to_filter_fwd, $reads_to_filter_rev) = $self->get_split_fastq_dir_and_filenames();
    my $filtered_reads_fwd = $self->{fsu}->catfile($temp_reads_dir, 'forward.for_assembly_improvement.fastq');
    my $filtered_reads_rev = $self->{fsu}->catfile($temp_reads_dir, 'reverse.for_assembly_improvement.fastq');
    my ($memory_required_mb, $num_threads) = $self->get_memory_and_threads();
    my $lane_paths_str = $self->get_lane_paths_str();
    my $umask_str = $self->umask_str;
    my $contigs_base_name = $self->generate_contig_base_name();
    my $script_name = $self->{fsu}->catfile($self->{lane_path}, $self->{prefix}.$self->{assembler}."_assembly_improvement.pl");
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");

    (my $script_string = qq{
        use strict;
        use File::Touch;
        use File::Path qw(make_path remove_tree);
        use $assembler_class;
        use VertRes::Pipelines::Assembly;
        use File::Path qw(make_path);
        use Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs;
        use Pathogens::Utils::FastqPooler;
        $umask_str;

        remove_tree(qq[$tmp_directory]) if(-d qq[$tmp_directory]);
        make_path(qq[$tmp_directory]) or die "Error mkdir $tmp_directory";
        chdir(qq[$tmp_directory]) or die "Error chdir $tmp_directory";

        if ($self->{improve_assembly}) {
            unless (-e qq[$reads_to_filter_fwd] && -e qq[$reads_to_filter_rev]) {
              my \@lane_paths = $lane_paths_str;
              Pathogens::Utils::FastqPooler->new(
                  output_directory => qq[$temp_reads_dir],
                  lane_paths => \\\@lane_paths
              )->run();
            }
            my \$assembly_pipeline = VertRes::Pipelines::Assembly->new(
              assembler => "$self->{assembler}"
            );
            my \$ok = \$assembly_pipeline->map_and_filter_perfect_pairs(qq[$fasta_for_improvement], qq[$pool_directory]);

            \$ok = \$assembly_pipeline->improve_assembly(qq[$fasta_for_improvement],[qq[$filtered_reads_fwd],qq[$filtered_reads_rev]],$insert_size,$num_threads);
            unlink(qq[$filtered_reads_fwd]) or die "Error unlink $filtered_reads_fwd";
            unlink(qq[$filtered_reads_rev]) or die "Error unlink $filtered_reads_rev";
        }

        # Rename contigs
        Bio::AssemblyImprovement::PrepareForSubmission::RenameContigs->new(input_assembly => qw[$fasta_for_improvement],base_contig_name => qq[$contigs_base_name])->run();

        chdir(qq[$output_directory]);
        remove_tree(qq[$tmp_directory]) or die "Error remove_tree $tmp_directory";
        my \$done_file = qq[$self->{lane_path}/$self->{prefix}$self->{assembler}_assembly_improvement_done];
        touch(\$done_file) or die "Error touch \$done_file"; #The prefix in the config file is not always the name of assembler, so we append assembler name
    }) =~ s/^ {8}//mg;

    print $scriptfh $script_string;
    close $scriptfh;

    my $job_name = $self->{prefix}.$self->{assembler}.'_assembly_improvement';

    my $total_memory_mb = 3000;
    my $queue = 'normal';

    VertRes::LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => "-n$num_threads -q $queue -M${total_memory_mb} -R 'select[mem>$total_memory_mb] rusage[mem=$total_memory_mb] span[hosts=1]'", dont_wait=>1}, qq{perl -w $script_name});

    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
    return $self->{No};
}

sub iva_qc_provides
{
  my $self = shift;
  if (defined($self->{iva_qc}) and $self->{iva_qc} and defined(qq[$self->{kraken_db}])) {
    return  [$self->{lane_path}."/".$self->{prefix}."$self->{assembler}_iva_qc_done"];
  }
  else {
    return [];
  }
}


sub iva_qc_requires
{
  my $self = shift;
  return $self->assembly_improvement_provides();
}


sub iva_qc
{
    my ($self, $build_path,$action_lock) = @_;

    unless (defined($self->{iva_qc}) and $self->{iva_qc} and defined(qq[$self->{kraken_db}])) {
        return $self->{Yes};
    }

    my $output_directory = $self->{lane_path};
    my $lane_names = $self->get_all_lane_names($self->{pools});
    my $assembly_directory = $self->{fsu}->catfile($output_directory, "$self->{assembler}_assembly");
    my $assembly_fasta = $self->{fsu}->catfile($assembly_directory, "contigs.fa");
    my $tmp_directory = $self->{tmp_directory}.'/'.$self->{prefix}.$self->{assembler}.'_'.$lane_names->[0] || getcwd();
    my ($split_reads_dir, $split_reads_fwd, $split_reads_rev) = $self->get_split_fastq_dir_and_filenames();
    my $lane_paths_str = $self->get_lane_paths_str();
    my $umask_str = $self->umask_str;

    # Run iva_qc if needed
    if(defined($self->{iva_qc}) && $self->{iva_qc} && defined(qq[$self->{kraken_db}])) {
        my $script_name = $self->{fsu}->catfile($self->{lane_path}, $self->{prefix}.$self->{assembler}."_iva_qc.pl");
        open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");

        (my $script_string = qq{
            use strict;
            use File::Touch;
            use File::Path qw(make_path remove_tree);
            use File::Copy::Recursive qw(dircopy);
            use Pathogens::Utils::FastqPooler;
            use Bio::AssemblyImprovement::IvaQC::Main;
            $umask_str
            remove_tree(qq[$tmp_directory]) if(-d qq[$tmp_directory]);
            make_path(qq[$tmp_directory]) or die "Error mkdir $tmp_directory";
            chdir(qq[$tmp_directory]) or die "Error chdir $tmp_directory";

            # If improve assembly has already been run, the forward and reverse fastq files should already be here
            # If not, rerun the split reads method
            unless (-e qq[$split_reads_fwd] && -e qq[$split_reads_rev]){
                my \@lane_paths = $lane_paths_str;
                Pathogens::Utils::FastqPooler->new(
                    output_directory => qq[$split_reads_dir],
                    lane_paths => \\\@lane_paths
                )->run();
            }
            my \$iva_qc = Bio::AssemblyImprovement::IvaQC::Main->new(
                    'db'      			  => qq[$self->{kraken_db}],
                    'forward_reads'       => qq[$split_reads_fwd],
                    'reverse_reads'       => qq[$split_reads_rev],
                    'assembly'			  => qq[$assembly_fasta],
                    'iva_qc_exec'         => qq[$self->{iva_qc_exec}],
                    'working_directory'	  => qq[$tmp_directory],
            );
            \$iva_qc->run();

            chdir(qq[$output_directory]);

            if (-e qq[$tmp_directory/iva_qc/iva_qc.stats.txt]) {
                dircopy(qq[$tmp_directory/iva_qc/], qq[$assembly_directory/iva_qc]);
            }
            else {
                my \$fail_file = qq[$self->{lane_path}/$self->{prefix}$self->{assembler}_iva_qc_failed];
                touch(\$fail_file) or die "Error touch \$fail_file";
            }

            # whether or not iva qc actually worked, we want to write the done
            # fileanyway so pipeline can continue
            remove_tree(qq[$tmp_directory]) or die "Error remove_tree $tmp_directory";
            my \$done_file = qq[$self->{lane_path}/$self->{prefix}$self->{assembler}_iva_qc_done];
            touch(\$done_file) or die "Error touch \$done_file"; #The prefix in the config file is not always the name of assembler, so we append assembler name
        }) =~ s/^ {12}//mg;

        my $job_name = $self->{prefix}.$self->{assembler}.'_iva_qc';
        my $num_threads = 1;
        my $queue = 'normal';
        my $total_memory_mb = 2000;
        VertRes::LSF::run($action_lock, $output_directory, $job_name, {bsub_opts => "-n$num_threads -q $queue -M${total_memory_mb} -R 'select[mem>$total_memory_mb] rusage[mem=$total_memory_mb] span[hosts=1]'", dont_wait=>1}, qq{perl -w $script_name});

        print $scriptfh $script_string;
        close $scriptfh;
    }
}


sub decide_appropriate_queue
{
  my ($self, $memory_required_mb) = @_;
  my $queue = 'long';
  if($memory_required_mb > 200000)
  {
    $queue = 'hugemem';
  }
  elsif($memory_required_mb < 12000 and $self->{assembler} ne 'iva')
  {
    $queue = 'normal';
  }
  return $queue;
}


##Improve assembly step. Runs SSPACE, abacas (if a reference is provided) and GapFiller.
## descaffolds the resulting assembly if necessary and cleans up any small contigs.

sub map_and_filter_perfect_pairs
{
  my($self,$reference, $working_directory)= @_;

  my $mapper = VertRes::Wrapper::smalt->new();
  $mapper->setup_custom_reference_index($reference,'-k 13 -s 4','small');

  `smalt map -x -i 3000 -f samsoft -y 0.95 -o $working_directory/contigs.mapped.sam $reference.small $working_directory/forward.fastq $working_directory/reverse.fastq`;
  $self->throw("Sam file not created") unless(-e "$working_directory/contigs.mapped.sam");

  `samtools faidx $reference`;
  $self->throw("Reference index file not created") unless(-e "$reference.fai");

  # Filter reads which are flagged as perfect pairs
  `samtools view -F 2 -bt $reference.fai $working_directory/contigs.mapped.sam > $working_directory/contigs.mapped.bam`;
  $self->throw("Couldnt convert from sam to BAM") unless(-e "$working_directory/contigs.mapped.bam");
  unlink("$working_directory/contigs.mapped.sam");

  `samtools sort -m 500000000 $working_directory/contigs.mapped.bam $working_directory/contigs.mapped.sorted`;
  $self->throw("Couldnt sort the BAM") unless(-e "$working_directory/contigs.mapped.sorted.bam");

  `samtools index $working_directory/contigs.mapped.sorted.bam`;

  VertRes::Utils::Sam->new(verbose => 1, quiet => 0)->bam2fastq(qq[$working_directory/contigs.mapped.sorted.bam], qq[subset]);
  unlink("$working_directory/contigs.mapped.sorted.bam");
  unlink("$reference.small.sma");
  unlink("$reference.small.smi");
  unlink("$reference.fai");
  `mv $working_directory/subset_1.fastq $working_directory/forward.for_assembly_improvement.fastq`;
  `mv $working_directory/subset_2.fastq $working_directory/reverse.for_assembly_improvement.fastq`;
}


sub improve_assembly
{
  my ($self,$assembly_file, $input_files, $insert_size, $threads) = @_;

  my $preprocess_input_files = Bio::AssemblyImprovement::Scaffold::SSpace::PreprocessInputFiles->new(
      input_files    => $input_files,
      input_assembly => $assembly_file,
      reference      => $self->{reference},
  );
  my $process_input_files_tmp_dir_obj = $preprocess_input_files->_temp_directory_obj();

  # scaffold and extend contigs
  my $scaffolding_obj = Bio::AssemblyImprovement::Scaffold::SSpace::Iterative->new(
      input_files     => $preprocess_input_files->processed_input_files,
      input_assembly  => $preprocess_input_files->processed_input_assembly,
      insert_size     => $insert_size,
      threads		  => $threads,
      scaffolder_exec => $self->{scaffolder_exec}
  );
  $scaffolding_obj->run();

  my $scaffolding_output = $scaffolding_obj->final_output_filename;

  # order contigs on an assembly
  if(defined($self->{reference}) && (-e $self->{reference}))
  {
    $scaffolding_obj = Bio::AssemblyImprovement::Abacas::Iterative->new(
      reference      => $preprocess_input_files->processed_reference,
      input_assembly => $scaffolding_output,
      abacas_exec    => $self->{abacas_exec},
    );
    $scaffolding_obj->run();
  }

  # fill gaps
  my $fill_gaps_obj = Bio::AssemblyImprovement::FillGaps::GapFiller::Iterative->new(
      input_files     => $preprocess_input_files->processed_input_files,
      input_assembly  => $scaffolding_obj->final_output_filename,
      insert_size     => $insert_size,
      threads		  => $threads,
      gap_filler_exec => $self->{gap_filler_exec},
      _output_prefix  => 'gapfilled'
  )->run();
  move($fill_gaps_obj->final_output_filename,$assembly_file);

  # descaffold if needed
  if(defined($self->{no_scaffolding}) && $self->{no_scaffolding} == 1)
  {
    my $descaffold_obj = Bio::AssemblyImprovement::Scaffold::Descaffold->new(input_assembly => $assembly_file);
    $descaffold_obj->run();
    move($descaffold_obj->output_filename,$assembly_file);
  }

  # remove all contigs shorter than a read length
  my $fastq_processor = Bio::AssemblyImprovement::Util::FastqTools->new(input_filename => $preprocess_input_files->processed_input_files->[0]);
  my $read_length = $fastq_processor->first_read_length();
  my $fasta_processor = Bio::AssemblyImprovement::Util::FastaTools->new(input_filename => $assembly_file);
  $fasta_processor->remove_small_contigs($read_length, 0);
  move($fasta_processor->output_filename,$assembly_file);

  # filter short contigs if needed
  if(defined($self->{post_contig_filtering}) && $self->{post_contig_filtering} > 0)
  {
    $fasta_processor = Bio::AssemblyImprovement::Util::FastaTools->new(input_filename => $assembly_file);
    $fasta_processor->remove_small_contigs($self->{post_contig_filtering}, 95);
    move($fasta_processor->output_filename,$assembly_file);
  }

  # order contigs by length
  my $order_contigs = Bio::AssemblyImprovement::Util::OrderContigsByLength->new( input_filename  => $assembly_file );
  $order_contigs->run();
  move($order_contigs->output_filename,$assembly_file);




}



sub generate_contig_base_name
{
  my ($self) = @_;
  my $lane_names = $self->get_all_lane_names($self->{pools});

  for my $lane_name (@{$lane_names})
  {
    my $vrlane  = VRTrack::Lane->new_by_name($self->{vrtrack}, $lane_name) or $self->throw("No such lane in the DB: [".$lane_name."]");

    if(defined($vrlane->acc()))
    {
      # use the first one available which has an accession number, normally there will only be 1
      return join('.',($vrlane->acc(),$lane_name));
    }
  }
  return join('.',('',$lane_names->[0]));
}

# Get the requested insert size of the first lane. Not suitable for mixed insert sizes, should be run with standalone scripts in that case.
sub get_insert_size
{
  my ($self) = @_;
  my $lane_names = $self->get_all_lane_names($self->{pools});
  my $insert_size = 350;

  for my $lane_name (@{$lane_names})
  {
    my $vrlane  = VRTrack::Lane->new_by_name($self->{vrtrack}, $lane_name) or $self->throw("No such lane in the DB: [".$lane_name."]");
    eval
    {
      $insert_size = VRTrack::Library->new($self->{vrtrack}, $vrlane->library_id())->insert_size() ;
    };
  }
  return $insert_size;
}


sub lane_read_length
{
  my ($self) = @_;
  my $lane_names = $self->get_all_lane_names($self->{pools});
  my $vrlane  = VRTrack::Lane->new_by_name($self->{vrtrack}, @$lane_names[0]) or $self->throw("No such lane in the DB: [".@$lane_names[0]."]");

  my $read_length = $vrlane->read_len();

  if((!defined($read_length)) || $read_length <= 0)
  {
    for my $file_name (@{$vrlane->files})
    {
      my $fastqcheck_filename = "$file_name.fastqcheck";
      next unless(-e $fastqcheck_filename);
      my $pars = VertRes::Parser::fastqcheck->new(file => $fastqcheck_filename);
      $read_length = $pars->max_length();
      last if($read_length > 0);
    }
  }

  if((!defined($read_length)) || $read_length <= 0)
  {
    $read_length = 36;
  }

  return $read_length;
}

sub number_of_threads
{
  my ($self, $memory_required_mb) = @_;
  my $normal_queue_mem_limit = 100000;
  my $num_threads = 1;

  if($normal_queue_mem_limit/$memory_required_mb > 2)
  {
    $num_threads = int($normal_queue_mem_limit/$memory_required_mb);
  }
  if($num_threads  > $self->{max_threads})
  {
    $num_threads = $self->{max_threads};
  }

  return $num_threads;
}

###########################
# End optimise
###########################




###########################
# Begin pool fastqs
###########################


sub pool_fastqs {
    my ( $self, $build_path, $action_lock ) = @_;

    my $lane_names       = $self->get_all_lane_names( $self->{pools} );
    my $output_directory = $self->{lane_path} . "/" . $self->{prefix} . $self->{assembler} . "_pool_fastq_tmp_files";
    my $base_path        = $self->{seq_pipeline_root};
    my $assembler        = $self->{assembler};
	my $umask    = $self->umask_str;

    my $script_name = $self->{fsu}->catfile( $build_path, $self->{prefix} . "pool_fastqs.pl" );
    open( my $scriptfh, '>', $script_name ) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
use strict;
use VertRes::Pipelines::Assembly;
use Bio::AssemblyImprovement::Assemble::SGA::Main;
use Bio::AssemblyImprovement::DigitalNormalisation::Khmer::Main;
use Bio::AssemblyImprovement::PrimerRemoval::Main;
use Bio::AssemblyImprovement::AdapterRemoval::Trimmomatic::Main;
use File::Copy;
$umask
my \$assembly= VertRes::Pipelines::Assembly->new(assembler => qq[$assembler]);
my \@lane_names;

system("rm -rf $output_directory");
system("mkdir -p $output_directory");
};

    my $lane_name = $self->{lane};
        my $lane_path = $self->{vrtrack}->hierarchy_path_of_lane_name($lane_name);
        my $vlane = VRTrack::Lane->new_by_name( $self->{vrtrack}, $lane_name );

        my @file_names;
        my @file_names_with_path;

        if ( @{ $vlane->files } == 2 ) {

            for my $file_name ( @{ $vlane->files } ) {
                push( @file_names, $file_name->name );
            }
            @file_names = sort @file_names
              ; # Unfortunately normalisation code cannot handle reads interleaved in any other way besides /1, /2, /1 and so on. We rely here on the files being named with _1 and _2 so that the order is maintained.
        }

        my $forward_reads_filename = $self->{lane_path} . "/" . $file_names[0];
        my $reverse_reads_filename = $self->{lane_path} . "/" . $file_names[1];

        # Adapter removal
        my $using_trimmomatic = (defined( $self->{remove_adapters} )
            and $self->{remove_adapters} == 1
            and defined( $self->{adapters_file} )
            and defined( $self->{adapter_removal_tool} )
            and $self->{adapter_removal_tool} eq 'trimmomatic'
            and $self->{assembler} ne 'iva');

        my $adpater_trimmed_reads_1 = "$output_directory/adapter_trimmed.forward.fastq.gz";
        my $adpater_trimmed_reads_2 = "$output_directory/adapter_trimmed.reverse.fastq.gz";

        if ($using_trimmomatic) {
            print $scriptfh qq{
my \$adapter_remover = Bio::AssemblyImprovement::AdapterRemoval::Trimmomatic::Main->new(
  'reads_in_1'       => "$forward_reads_filename",
  'reads_in_2'       => "$reverse_reads_filename",
  'paired_out_1'     => "$adpater_trimmed_reads_1",
  'paired_out_2'     => "$adpater_trimmed_reads_2",
  'trimmomatic_exec' => "$self->{trimmomatic_jar}",
  'adapters_file'    => "$self->{adapters_file}",
)->run();
};
            $forward_reads_filename = $adpater_trimmed_reads_1;
            $reverse_reads_filename = $adpater_trimmed_reads_2;
        }

        #Primer removal
        my $using_quasr = (defined( $self->{remove_primers} )
          and $self->{remove_primers} == 1
          and defined( $self->{primers_file} )
          and defined( $self->{primer_removal_tool} )
          and $self->{primer_removal_tool} eq 'quasr'
          and $self->{assembler} ne 'iva');

        if ($using_quasr) {
            #Replace code here with an alternative way of removing primers that can accept a shuffled file
            print $scriptfh qq{
my \$primer_remover = Bio::AssemblyImprovement::PrimerRemoval::Main->new(
  forward_file     => "$forward_reads_filename",
  reverse_file     => "$reverse_reads_filename",
  primers_file     => "$self->{primers_file}",
  output_directory => "$output_directory",
  QUASR_exec       => "$self->{QUASR_exec}",
)->run();
};

            $forward_reads_filename = $output_directory . "/" . 'primer_removed.forward.fastq.gz';
            $reverse_reads_filename = $output_directory . "/" . 'primer_removed.reverse.fastq.gz';
        }

        # Create a shuffled sequence. This shuffled file will be the input for any processing steps below (i.e. normalisation, error correction etc)
        my $shuffled_filename = $output_directory . '/' . $lane_name . '.fastq.gz';
        my $output_filename   = $lane_name
          . '.fastq.gz'
          ; # Each step below (normalisation and error correction), should produce an output file with this name (which is the same as the shuffled filename)

        print $scriptfh qq{
\$assembly->shuffle_sequences_fastq_gz("$forward_reads_filename","$reverse_reads_filename", "$shuffled_filename");
};

#Clean up primer removed files. This is quite messy. Will be better when primer removal can accept a shuffled file so it fits in like normalisation and error_correction does
        if ($using_quasr) {
            print $scriptfh qq{
unlink("$output_directory/primer_removed.forward.fastq.gz");
unlink("$output_directory/primer_removed.reverse.fastq.gz");
};
        }

        if ($using_trimmomatic) {
            print $scriptfh qq{
unlink("$adpater_trimmed_reads_1");
unlink("$adpater_trimmed_reads_2");
};
        }

        # Digital normalisation
        if ( defined( $self->{normalise} ) and $self->{normalise} == 1 ) {
            print $scriptfh qq{
my \$diginorm = Bio::AssemblyImprovement::DigitalNormalisation::Khmer::Main->new(
input_file       => "$shuffled_filename",
khmer_exec       => "$self->{khmer_exec}",
output_filename  => "$output_filename",
output_directory => "$output_directory",
)->run();
};
        }

        # Error correction
        if ( defined( $self->{error_correct} ) and $self->{error_correct} == 1 ) {
            print $scriptfh qq{
my \$sga = Bio::AssemblyImprovement::Assemble::SGA::Main->new(
input_files     => ["$shuffled_filename"],
output_filename => "$output_filename",
output_directory => "$output_directory",
pe_mode		    => 2,
sga_exec        => "$self->{sga_exec}",
)->run();
};
        }

    print $scriptfh qq{
system("mv $shuffled_filename $output_directory/pool_1.fastq.gz");
unlink("$output_directory/$lane_name.fastq.gz");
system("touch $self->{lane_path}/$self->{prefix}$self->{assembler}_pool_fastqs_done");
exit;
      };
    close $scriptfh;
    my $job_name = $self->{prefix} . 'pool_fastqs';

    my $memory_in_mb = 500;
    my $queue        = 'normal';
    if ( defined( $self->{error_correct} ) and $self->{error_correct} == 1 ) {
        $memory_in_mb = 4000;
        $queue        = 'long';
    }

    VertRes::LSF::run(
        $action_lock, $self->{lane_path}, $job_name,
        { bsub_opts => "-q $queue -M${memory_in_mb} -R 'select[mem>$memory_in_mb] rusage[mem=$memory_in_mb]'" },
        qq{perl -w $script_name}
    );

    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
    return $self->{No};
}


sub pool_fastqs_requires
{
  [];
}

sub pool_fastqs_provides
{
   my $self = shift;
   my @provided_files ;
   push(@provided_files, $self->{lane_path}.'/'.$self->{prefix}.$self->{assembler}.'_pool_fastqs_done');

   return \@provided_files;
}

sub get_all_lane_names
{
  my ($self, $pooled_lanes) = @_;
  my @all_lane_names ;

  for my $lane_pool (@$pooled_lanes)
  {
    for my $lane_name (@{$lane_pool->{lanes}})
    {
      push(@all_lane_names, $lane_name);
    }
  }
  return \@all_lane_names;
}


# adapted from https://github.com/dzerbino/velvet/blob/master/shuffleSequences_fastq.pl
sub shuffle_sequences_fastq_gz
{
  my ($self, $input_file_1, $input_file_2, $output_file_name) = @_;


  open( my $FILEA, "-|",'gunzip -c '.$input_file_1);
  open( my $FILEB, "-|",'gunzip -c '.$input_file_2);
  open( my $OUTFILE, "|-", "gzip -c  > $output_file_name");

  # FIXME: if FileB contains more reads than fileA they will get missed! This is also a bug in Velvets shuffleSequences_fastq.pl
  while(<$FILEA>) {
    print $OUTFILE $_;
    $_ = <$FILEA>;
    print $OUTFILE $_;
    $_ = <$FILEA>;
    print $OUTFILE $_;
    $_ = <$FILEA>;
    print $OUTFILE $_;

    my $file_b_line = <$FILEB>;
    next unless(defined($file_b_line));
    print $OUTFILE $file_b_line;
    $_ = <$FILEB>;
    print $OUTFILE $_;
    $_ = <$FILEB>;
    print $OUTFILE $_;
    $_ = <$FILEB>;
    print $OUTFILE $_;
  }
  close($FILEA);
  close($FILEB);
  close($OUTFILE);
}

###########################
# End pool fastqs
###########################



sub total_number_of_reads
{
  my ($self) = @_;
  my $lane_names = $self->get_all_lane_names($self->{pools});
  my $total_reads = 0;

  for my $lane_name (@{$lane_names})
  {
    my $vrlane  = VRTrack::Lane->new_by_name($self->{vrtrack}, $lane_name) or $self->throw("No such lane in the DB: [".$lane_name."]");
    $total_reads += $vrlane->raw_reads();
  }
  return $total_reads;
}


#Gives the answer in kb.
sub estimate_memory_required
{
  my ($self, $output_directory, $kmer_size) = @_;

  my %memory_params ;
  $memory_params{total_number_of_reads} = $self->total_number_of_reads();
  $memory_params{genome_size}           = $self->{genome_size};
  $memory_params{read_length}           = $self->lane_read_length();
  $memory_params{kmer_size}             = $kmer_size;
  $memory_params{error_correct}         = defined ($self->{error_correct}) ? $self->{error_correct}:0;

  my $assembler_class = $self->{assembler_class};
  eval("use $assembler_class; ");
  my $assembler_util= $assembler_class->new(output_directory => $output_directory);
  my $memory_required_in_kb = $assembler_util->estimate_memory_required(\%memory_params);
  return $memory_required_in_kb;
}


=head2 cleanup_requires

 Title   : cleanup_requires
 Usage   : my $required_files = $obj->cleanup_requires('/path/to/lane');
 Function: Find out what files the cleanup action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_requires {
  my ($self) = @_;
  return [ $self->{prefix}."assembly_update_db_done"];
}

=head2 cleanup_provides

 Title   : cleanup_provides
 Usage   : my $provided_files = $obj->cleanup_provides('/path/to/lane');
 Function: Find out what files the cleanup action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_provides {
    my ($self) = @_;
    return [$self->{prefix}."assembly_cleanup_done"];
}

=head2 cleanup

 Title   : cleanup
 Usage   : $obj->cleanup('/path/to/lane', 'lock_filename');
 Function: Unlink all the pipeline-related files (_*) in a lane, as well
           as the split directory.
 Returns : $VertRes::Pipeline::Yes
 Args    : lane path, name of lock file to use

=cut

sub cleanup {
  my ($self, $lane_path, $action_lock) = @_;
#  return $self->{Yes} unless $self->{do_cleanup};

  my $all_files_prefix = $self->{fsu}->catfile($self->{lane_path}, $self->{prefix}.$self->{assembler});

  # remove job files
  for my $action (@$actions) {
    my $action_prefix;
    if ($action->{name} eq 'pool_fastqs') {
      $action_prefix = $self->{fsu}->catfile($self->{lane_path}, $self->{prefix}. "pool_fastqs");
    }
    else {
      $action_prefix = $all_files_prefix . '_' . $action->{name};
    }

    my $iva_qc_fail_file = $all_files_prefix . "_iva_qc_failed";
    if ($action->{name} eq 'iva_qc' and -e $iva_qc_fail_file) {
      next;
    }

    for my $suffix (qw/o e pl/) {
      my $filename = "$action_prefix.$suffix";
      if (-e $filename) {
        unlink($filename) or $self->throw("Error unlink $filename");
      }
    }
  }

  # remove the tmp directory that was storing the pooled reads etc
  my $pool_directory = $self->{lane_path}."/".$self->{prefix}.$self->{assembler}."_pool_fastq_tmp_files";
  if (-d $pool_directory) {
      remove_tree($pool_directory) or $self->throw("Error remove_tree $pool_directory");
  }

  # remove the tmp directory from lustre
  my $lane_names = $self->get_all_lane_names($self->{pools});
  my $tmp_directory = $self->{tmp_directory}.'/'.$self->{prefix}.$self->{assembler}.'_'.$lane_names->[0];

  if (-e $tmp_directory) {
      remove_tree($tmp_directory) or $self->throw("Error remove_tree $tmp_directory");
  }

  # remove all the other unwanted files
  my @unwanted_files = qw/
    before_rr.fasta
    contigs.paths
    input_dataset.yaml
    split_input
    tmp
    contigs.fa.scaffolded.filtered
    .RData
    contigs.fa.png.Rout
    scaffolded.summaryfile.txt
    forward.fastq
    reverse.fastq
  /;

  my $assembly_dir = $self->{fsu}->catfile($self->{lane_path}, $self->{assembler}.'_assembly');

  foreach my $file (@unwanted_files) {
    my $to_remove = $self->{fsu}->catfile($assembly_dir, $file);

    if (-e $to_remove) {
      if (-d $to_remove) {
        remove_tree($to_remove) or $self->throw("Error remove_tree $to_remove");
      }
      else {
        unlink($to_remove) or $self->throw("Error unlink $to_remove");
      }
    }
  }

  Utils::CMD("touch ".$self->{fsu}->catfile($lane_path,"$self->{prefix}assembly_cleanup_done")   );
  $self->update_file_permissions($lane_path);
  return $self->{Yes};
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
    return $self->map_back_provides();
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
    return [ $self->{lane_path}."/".$self->{prefix}."assembly_update_db_done"];
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

    unless($vrlane->is_processed('assembled')){
      $vrtrack->transaction_start();
      $vrlane->is_processed('assembled',1);
      $vrlane->update() || $self->throw("Unable to set assembled status on lane $lane_path");
      $vrtrack->transaction_commit();
    }


    my $job_status =  File::Spec->catfile($lane_path, $self->{prefix} . 'job_status');
    Utils::CMD("rm $job_status") if (-e $job_status);
    Utils::CMD("touch ".$self->{fsu}->catfile($lane_path,"$self->{prefix}assembly_update_db_done")   ) unless( -e "$self->{prefix}assembly_update_db_done");

    return $$self{'Yes'};
}



=head2 get_split_fastq_dir_and_filenames

 Title   : get_split_fastq_dir_and_filenames
 Usage   : $obj->get_split_fastq_dir_and_filenames();
 Function: Determines directory name where forward.fastq and reverse.fastq split reads files live
 Returns : directory name, forward reads filename, reverse reads filename
 Args    : None

=cut
sub get_split_fastq_dir_and_filenames {
    my ($self) = @_;
    my $dir = $self->{lane_path}."/".$self->{prefix}.$self->{assembler}."_pool_fastq_tmp_files";
    unless (-d qq[$dir]) {
        make_path(qq[$dir]) or $self->throw("Error make_path $dir: $!");
    }
    return ($dir, "$dir/forward.fastq", "$dir/reverse.fastq");
}



=head2 get_lane_paths_str

 Title   : get_lane_paths_str
 Usage   : $obj->get_lane_paths_str();
 Function: Determines lane paths from path to lane
 Returns : lane paths string
 Args    : None

=cut
sub get_lane_paths_str {
    my ($self) = @_;
    my $lane_names = $self->get_all_lane_names($self->{pools});
    my @lane_paths;
    my $base_path = $self->{seq_pipeline_root};
    for my $lane_name (@$lane_names) {
      push(@lane_paths, $base_path.'/'.$self->{vrtrack}->hierarchy_path_of_lane_name($lane_name).'/'.$lane_name);
    }
    my $lane_paths_str = '("' . join('","', @lane_paths) . '")';
    return \@lane_paths, $lane_paths_str;
}




=head2 get_memory_and_threads

 Title   : get_memory_and_threads
 Usage   : $obj->get_memory_and_threads();
 Function: Determines assembly memory and threads
 Returns : memory, threads
 Args    : None

=cut
sub get_memory_and_threads {
    my ($self) = @_;
    my $read_length = $self->lane_read_length();
    my $kmer_for_memory_calculation =  int($read_length*0.33);
    my $memory_required_mb = int($self->estimate_memory_required($self->{lane_path}, $kmer_for_memory_calculation)/1000);
    my $num_threads = $self->number_of_threads($memory_required_mb);
    return ($memory_required_mb, $num_threads);
}
