=head1 NAME

VertRes::Utils::Mapping - mapping utility functions

=head1 SYNOPSIS

use VertRes::Utils::Mapping;

# get the mapping utility object appropriate for the lane's technology
my $chooser = VertRes::Utils::Mapping->new();
my $class = $chooser->lane_to_module('/path/to/SLX/lane');
# $class is VertRes::Utils::Mappers::SLX
require $class;
my $mapping_util = $class->new();

# use any of the utility functions described here, eg.
$mapping_util->split_fastq(read1 => 'reads_1.fastq',
                           read2 => 'reads_2.fastq',
                           split_dir => '/path/to/desired/split_dir',
                           chunk_size => 1000000);

=head1 DESCRIPTION

Lets you do mapping-related things without worring about what technology your
lane was done with.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Utils::Mapping;

use strict;
use warnings;
use VertRes::Utils::FastQ;
use VertRes::IO;
use VertRes::Parser::flagstat;
use SamTools;
use File::Basename;

use base qw(VertRes::Base);

our %tech_to_module = (454 => 'VertRes::Utils::Mappers::ssaha',
                       SLX => 'VertRes::Utils::Mappers::bwa');

our %do_mapping_args = (insert_size => 1);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Mapping->new();
 Function: Create a new VertRes::Utils::Mapping object.
 Returns : VertRes::Utils::Mapping object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 lane_to_module

 Title   : lane_to_module
 Usage   : my $module = $obj->lane_to_module('/path/to/lane');
 Function: Find out what mapping utility module to use on your lane.
 Returns : class string (call new on it)
 Args    : path string

=cut

sub lane_to_module {
    my ($self, $arg) = @_;
    
    if ($arg =~ /\/SLX\//i) {
        return $tech_to_module{SLX};
    }
    elsif ($arg =~ /\/454\//) {
        return $tech_to_module{454};
    }
    else {
        $self->throw("Encountered an argument that doesn't correspond to a technology: $arg");
    }
}

=head2 split_fastq

 Title   : split_fastq
 Usage   : $obj->split_fastq(read1 => 'reads_1.fastq',
                             read2 => 'reads_2.fastq',
                             split_dir => '/path/to/desired/split_dir',
                             chunk_size => 1000000);
 Function: Split the fastq(s) into multiple smaller files. This is just a
           convienience alias to VertRes::Utils::FastQ::split, with syntax
           more similar to do_mapping().
 Returns : int (the number of splits created)
 Args    : read1 => 'reads_1.fastq', read2 => 'reads_2.fastq'
           -or-
           read0 => 'reads.fastq'

           split_dir => '/path/to/desired/split_dir'
           chunk_size => int (max number of bases per chunk, no default - throws)

=cut

sub split_fastq {
    my ($self, %args) = @_;
    my %read_args = $self->_do_read_args(%args);
    my $chunk_size = $args{chunk_size} || $self->throw("chunk_size is required");
    my $split_dir = $args{split_dir} || $self->throw("split_dir must be supplied");
    
    my @fastq_files = values %read_args;
    
    my $fastq_util = VertRes::Utils::FastQ->new(verbose => $self->verbose);
    return $fastq_util->split(\@fastq_files, split_dir => $split_dir, chunk_size => $chunk_size);
}

=head2 wrapper

 Title   : wrapper
 Usage   : do not call here; this is supposed to be overriden
 Function: Get a wrapper to actually do some mapping with.
 Returns : VertRes::Wrapper::WrapperI-based object (call do_mapping() on it)
 Args    : n/a

=cut

sub wrapper {
    my $self = shift;
    $self->throw("This is supposed to be overriden");
}

=head2 do_mapping

 Title   : do_mapping
 Usage   : $obj->do_mapping(ref => 'ref.fa',
                            read1 => 'reads_1.fastq',
                            read2 => 'reads_2.fastq',
                            output => 'output.sam',
                            insert_size => 2000);
 Function: A convienience method that calls do_mapping() on the return value of
           wrapper(), translating generic options to those suitable for the
           wrapper. Also converts the output to a sam file if that isn't the
           default format of the mapper.
 Returns : boolean (true on success)
 Args    : required options:
           ref => 'ref.fa'
           output => 'output.sam'

           read1 => 'reads_1.fastq', read2 => 'reads_2.fastq'
           -or-
           read0 => 'reads.fastq'

           and optional generic options:
           insert_size => int (default 2000)

=cut

sub do_mapping {
    my $self = shift;
    $self->throw("This is supposed to be overriden");
}

# Mappers should override this as appropriate
sub _bsub_opts {
    my ($self, $lane_path, $action) = @_;
    
    # by default, no extra bsub options will be supplied to LSF::Run
    return {bsub_opts => ''};
}

=head2 _do_mapping_args

 Title   : _do_mapping_args
 Usage   : my %args_for_my_wrapper_do_mapping = $obj->_do_mapping_args(
               \%conversion_hash, @user_args);
 Function: Internal method of Mapper module authors; convert the generic
           do_mapping args of this utility module to the wrapper-specific
           args of your Mapper.
 Returns : hash of args
 Args    : a hash ref with keys as the generic optional args understood by
           VertRes::Utils::Mapping::do_mapping(), and values as the args
           understood by your wrapper's do_mapping(). And the user's args as
           an associative array.

=cut

sub _do_mapping_args {
    my ($self, $converter, @args) = @_;
    my %args = (insert_size => 2000, @args);
    
    my %out_hash = $self->_do_read_args(@args);
    $self->throw("ref is required") unless defined $args{ref};
    $self->throw("ref file ($args{ref}) must exist") unless -s $args{ref};
    $self->throw("output is required") unless defined $args{output};
    
    $out_hash{ref} = $args{ref};
    $out_hash{output} = $args{output};
    
    foreach my $arg (keys %do_mapping_args) {
        if (defined $args{$arg}) {
            $out_hash{$converter->{$arg}} = $args{$arg};
        }
    }
    
    return %out_hash;
}

sub _do_read_args {
    my ($self, %args) = @_;
    
    my %out_hash;
    
    if (defined $args{read0}) {
        $self->throw("read0 and read1/2 are mutually exclusive") if (defined $args{read1} || defined $args{read2});
        $self->throw("read0 file ($args{read0}) must exist") unless -s $args{read0};
        $out_hash{read0} = $args{read0};
    }
    elsif (defined $args{read1}) {
        $self->throw("read2 must be supplied with read1") unless defined $args{read2};
        $self->throw("read1 file ($args{read1}) must exist") unless -s $args{read1};
        $self->throw("read2 file ($args{read2}) must exist") unless -s $args{read2};
        $out_hash{read1} = $args{read1};
        $out_hash{read2} = $args{read2};
    }
    elsif (defined $args{read2}) {
        $self->throw("read1 must be supplied with read2");
    }
    
    return %out_hash;
}

=head2 mapping_hierarchy_report

 Title   : mapping_hierarchy_report
 Usage   : $obj->mapping_hierarchy_report('output.csv', @paths_to_lanes);
 Function: Create a summary report for all lanes, describing how the mapping
           went.
 Returns : n/a
 Args    : output file, list of absolute paths to lanes

=cut

sub mapping_hierarchy_report {
    my ($self, $output_csv, @lane_paths) = @_;
    
    # write a header for our output csv
    open(my $csvfh, '>', $output_csv) or $self->throw("Cannot create output file: $output_csv");
    print $csvfh "Status,Center,Project,Individual,Tech,Library,Lane,Fastq0,ReadLength,Fastq1,ReadLength,Fastq2,ReadLength,#Reads,#Bases,#RawReadsMapped,#RawBasesMapped,#RawReadsPaired,#RmdupReadsMapped,#RmdupBasesMapped,ErrorRate\n";
    
    # print mapping stats for each lane
    foreach my $lane_path (@lane_paths) {
        if (! -d $lane_path) {
            $self->warn("Cant find lane directory: $lane_path");
            next;
        }
        
        # get the fastq info
        my @fastq_info = @{HierarchyUtilities::getFastqInfo($lane_path)};
        my $num_bases = 0;
        my $num_reads = 0;
        for my $i (0..$#fastq_info) {
            $num_reads += $fastq_info[$i][2];
            $num_bases += $fastq_info[$i][3];
        }
        my $fastq_csv = "$fastq_info[0][0],$fastq_info[0][1],$fastq_info[1][0],$fastq_info[1][1],$fastq_info[2][0],$fastq_info[2][1]";
        
        # get the lane info
        my %lane_info = %{HierarchyUtilities::lane_info($lane_path)};
        my $laneinfo_csv = "$lane_info{project},$lane_info{sample},$lane_info{technology},$lane_info{library},$lane_info{lane},$fastq_csv,$num_reads,$num_bases";
        
        my $io = VertRes::IO->new();
        my $raw_bam = $io->catfile($lane_path, 'raw.sorted.bam');
        my $rmdup_bam = $io->catfile($lane_path, 'rmdup.bam');
        
        if (-s $raw_bam && -s $rmdup_bam) {
            my @stats = $self->get_mapping_stats($raw_bam, $rmdup_bam);
            
            print $csvfh "MAPPED,,$laneinfo_csv";
            foreach (@stats) {
                print $csvfh ",$_";
            }
            print $csvfh "\n";
        }
        else {
            print $csvfh "UNMAPPED,,$laneinfo_csv\n";
        }
    }
    close($csvfh);
    
    return;
}

=head2 get_mapping_stats

 Title   : get_mapping_stats
 Usage   : my @stats = $obj->get_mapping_stats('raw.bam', 'rmdup.bam');
 Function: Get certain mapping stats for a bam file and its rmdup. Assumes
           VertRes::Utils::Sam->stats() has already been run on the two bam
           files.
 Returns : list of these stats: number of mapped reads, number of mapped bases,
           number of paired reads, number of mapped reads after removing
           duplicates, number of mapped bases after removing duplicates,
           error rate.
           NB: error rate is always 0, since this is a hold-over from maq, when
           maq mapstat was used on .map files, which told the error rate; we
           have no way of determining this with bam files.
 Args    : bam file, rmdup'd version of the bam file (both must have been
           supplied to stats() already)

=cut

sub get_mapping_stats {
    my ($self, $raw_bam, $rmdup_bam) = @_;
    
    my $fsparser = VertRes::Parser::flagstat->new(file => $raw_bam.'.flagstat');
    my $mapped_reads = $fsparser->mapped_reads();
    my $paired_reads = $fsparser->mapped_proper_paired_reads(); #*** or is this supposed to be mapped_paired_reads() ?
    $fsparser->file($rmdup_bam.'.flagstat');
    my $rmdup_mapped_reads = $fsparser->mapped_reads();
    
    # correct number of mapped bases isn't reported in either flagstat or
    # bamstat files; we have to calculate it
    my $stats = SamTools::collect_detailed_bam_stats($raw_bam, undef, {do_chrm => 0, do_rmdup => 0});
    my $mapped_bases = $stats->{bases_mapped};
    $stats = SamTools::collect_detailed_bam_stats($rmdup_bam, undef, {do_chrm => 0, do_rmdup => 0});
    my $rmdup_mapped_bases = $stats->{bases_mapped};
    
    return ($mapped_reads, $mapped_bases, $paired_reads, $rmdup_mapped_reads, $rmdup_mapped_bases, 0);
}

1;
