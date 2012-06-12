package Pathogens::Reports::Mapping::Row ;

use Moose;
use VRTrack::VRTrack; # Includes Lane, Mapstats, Etc. 

has 'vrtrack'             => ( is => 'ro', isa => 'VRTrack::VRTrack',  required => 1); # database
has 'lane'                => ( is => 'ro', isa => 'VRTrack::Lane',     required => 1); # lane
has 'mapstats'            => ( is => 'ro', isa => 'VRTrack::Mapstats', required => 1); # mapstats

# Checks
has 'is_qc_mapstats'      => ( is => 'ro', isa => 'Bool',              lazy_build => 1); # qc or mapping mapstats.
has 'is_mapping_complete' => ( is => 'ro', isa => 'Maybe[Bool]',       lazy_build => 1); # Mapping completed

# Internals
has '_vrtrack_project'    => ( is => 'ro', isa => 'VRTrack::Project',  lazy_build => 1); # Assembly - from mapststs 
has '_vrtrack_sample'     => ( is => 'ro', isa => 'VRTrack::Sample',   lazy_build => 1); # Mapper - from mapstats
has '_vrtrack_assembly'   => ( is => 'ro', isa => 'VRTrack::Assembly', lazy_build => 1); # Assembly - from mapststs 
has '_vrtrack_mapper'     => ( is => 'ro', isa => 'VRTrack::Mapper',   lazy_build => 1); # Mapper - from mapstats

# Cells
has 'study_id'            => ( is => 'ro', isa => 'Int',        lazy_build => 1); # study ssid
has 'sample'              => ( is => 'ro', isa => 'Str',        lazy_build => 1); # sample name
has 'lanename'            => ( is => 'ro', isa => 'Str',        lazy_build => 1); # lane name
has 'cycles'              => ( is => 'ro', isa => 'Int',        lazy_build => 1); # cycles/readlength
has 'reads'               => ( is => 'ro', isa => 'Int',        lazy_build => 1); # lane yield (reads)
has 'bases'               => ( is => 'ro', isa => 'Int',        lazy_build => 1); # lane yield (bases)
has 'map_type'            => ( is => 'ro', isa => 'Str',        lazy_build => 1); # qc or mapping reported value
has 'reference'           => ( is => 'ro', isa => 'Str',        lazy_build => 1); # reference name
has 'reference_size'      => ( is => 'ro', isa => 'Int',        lazy_build => 1); # reference size
has 'mapper'              => ( is => 'ro', isa => 'Str',        lazy_build => 1); # mapper name
has 'adapter_perc'        => ( is => 'rw', isa => 'Maybe[Num]', lazy_build => 1); # percent adaptor reads
has 'transposon_perc'     => ( is => 'rw', isa => 'Maybe[Num]', lazy_build => 1); # percent transposon reads
has 'mapped_perc'         => ( is => 'ro', isa => 'Num',        lazy_build => 1); # percent mapped reads
has 'paired_perc'         => ( is => 'ro', isa => 'Num',        lazy_build => 1); # percent paired reads
has 'mean_insert_size'    => ( is => 'ro', isa => 'Maybe[Num]', lazy_build => 1); # mean insert size
has 'genome_covered'      => ( is => 'ro', isa => 'Maybe[Num]', lazy_build => 1); # genome covered QC
has 'genome_covered_1x'   => ( is => 'ro', isa => 'Maybe[Num]', lazy_build => 1); # genome covered Mapping
has 'genome_covered_5x'   => ( is => 'ro', isa => 'Maybe[Num]', lazy_build => 1); # ditto
has 'genome_covered_10x'  => ( is => 'ro', isa => 'Maybe[Num]', lazy_build => 1); # ditto
has 'genome_covered_50x'  => ( is => 'ro', isa => 'Maybe[Num]', lazy_build => 1); # ditto
has 'genome_covered_100x' => ( is => 'ro', isa => 'Maybe[Num]', lazy_build => 1); # ditto
has 'depth_of_coverage'   => ( is => 'ro', isa => 'Maybe[Num]', lazy_build => 1); # mean depth of coverage
has 'depth_of_coverage_sd'=> ( is => 'ro', isa => 'Maybe[Num]', lazy_build => 1); # mean depth of coverage sd
has 'duplication_rate'    => ( is => 'ro', isa => 'Maybe[Num]', lazy_build => 1); # duplication rate qc
has 'error_rate'          => ( is => 'ro', isa => 'Maybe[Num]', lazy_build => 1); # error rate qc
has 'npg_qc'              => ( is => 'ro', isa => 'Maybe[Str]', lazy_build => 1); # npg qc
has 'manual_qc'           => ( is => 'ro', isa => 'Maybe[Str]', lazy_build => 1); # manual qc

# Is mapstats entry from QC or Mapping
sub _build_is_qc_mapstats
{
    my ($self) = @_;
    return $self->{mapstats}->is_qc();
}

sub _build_is_mapping_complete
{
    my ($self) = @_;
    return undef unless defined $self->mapstats->bases_mapped; # Mapping not complete or failed;
    return $self->mapstats->bases_mapped ? 1:0;                # Mapping complete and bases found;
}

# Build Internals
sub _build__vrtrack_project
{
    my($self) = @_;
    my $project;
    eval
    {
	$project = VRTrack::Project->new($self->vrtrack,$self->_vrtrack_sample->project_id);
    };
    return $project;
}

sub _build__vrtrack_sample
{
    my($self) = @_;
    my $sample;
    eval
    {
	my $library = VRTrack::Library->new($self->vrtrack,$self->lane->library_id);
	$sample = VRTrack::Sample->new($self->vrtrack,$library->sample_id) if defined $library;
    };
    return $sample;
}

sub _build__vrtrack_assembly
{
    my($self) = @_;
    my $assembly = $self->mapstats->assembly() or die "debug: build assembly failed\n";
    return $assembly; 
}

sub _build__vrtrack_mapper
{
    my($self) = @_;
    my $mapper = $self->mapstats->mapper() or die "debug: build mapper failed\n";
    return $mapper; 
}

# Build Cells
sub _build_study_id
{
    my($self) = @_;
    return $self->_vrtrack_project->ssid();
}

sub _build_sample
{
    my($self) = @_;
    return $self->_vrtrack_sample->name();
}

sub _build_lanename
{
    my($self) = @_;
    return $self->lane->name();
}

sub _build_cycles
{
    my($self) = @_;
    return $self->lane->read_len();
}

sub _build_reads
{
    my($self) = @_;
    return $self->lane->raw_reads();
}

sub _build_bases
{
    my($self) = @_;
    return $self->lane->raw_bases();
}

sub _build_map_type
{
    my($self) = @_;
    return $self->is_qc_mapstats() ? 'QC':'Mapping';
}

sub _build_reference
{
    my($self) = @_;
    return $self->_vrtrack_assembly->name();
}

sub _build_reference_size
{
my($self) = @_;
return $self->_vrtrack_assembly->reference_size();
}

sub _build_mapper
{
    my($self) = @_;
    return $self->mapstats->mapper->name();
}

sub _build_adapter_perc
{
    my($self) = @_;
    return undef unless $self->is_qc_mapstats; # not defined for mapping.

    my $adapter_perc;
    if( defined $self->mapstats->adapter_reads && $self->mapstats->raw_reads)
    {
	my $adapter_reads = $self->mapstats->adapter_reads;
	$adapter_perc  = sprintf("%.1f",($adapter_reads/$self->mapstats->raw_reads)*100);
    }
    return $adapter_perc;
}

sub _build_transposon_perc
{
    my($self) = @_;
    return undef unless $self->is_qc_mapstats; # not defined for mapping.

    my $transposon_perc;
    if(defined $self->mapstats->percentage_reads_with_transposon)
    {
	$transposon_perc = sprintf("%.1f", $self->mapstats->percentage_reads_with_transposon);
    }
    return $transposon_perc;
}

sub _build_mapped_perc
{
    my($self) = @_;
    my $reads_mapped_perc = '0.0';
    if( $self->is_mapping_complete )
    {
	my $reads_mapped = $self->mapstats->reads_mapped;
	my $raw_reads    = $self->mapstats->raw_reads;
	$reads_mapped_perc = sprintf("%.1f", ($reads_mapped/$raw_reads)*100);
    }
    return $reads_mapped_perc;
}

sub _build_paired_perc
{
    my($self) = @_;
    my $reads_paired_perc = '0.0';
    if( $self->is_mapping_complete )
    {
	my $reads_paired = $self->mapstats->reads_paired;
	my $raw_reads    = $self->mapstats->raw_reads;
	$reads_paired_perc = sprintf("%.1f", ($reads_paired/$raw_reads)*100);
    }
    return $reads_paired_perc;
}

sub _build_mean_insert_size
{
    my($self) = @_;
    return $self->mapstats->mean_insert;
}

sub _build_genome_covered
{
    my($self) = @_;
    return undef unless $self->is_qc_mapstats; # Not calculated for mapping

    my $genome_cover_perc; 
    if($self->is_mapping_complete)
    {
	my $target_bases_mapped = $self->mapstats->target_bases_mapped;
	my $genome_covered = $self->reference_size ? ($target_bases_mapped/$self->reference_size)*100 : undef;
	$genome_cover_perc = sprintf("%5.2f", $genome_covered) if defined $genome_covered;
    }
    return $genome_cover_perc;
}

sub _build_genome_covered_1x
{
    my($self) = @_;
    return undef if $self->is_qc_mapstats;
    return $self->_target_bases_X_perc(1);
}

sub _build_genome_covered_5x
{
    my($self) = @_;
    return undef if $self->is_qc_mapstats;
    return $self->_target_bases_X_perc(5);
}

sub _build_genome_covered_10x
{
    my($self) = @_;
    return undef if $self->is_qc_mapstats;
    return $self->_target_bases_X_perc(10);
}

sub _build_genome_covered_50x
{
    my($self) = @_;
    return undef if $self->is_qc_mapstats;
    return $self->_target_bases_X_perc(50);
}

sub _build_genome_covered_100x
{
    my($self) = @_;
    return undef if $self->is_qc_mapstats;
    return $self->_target_bases_X_perc(100);
}

sub _build_depth_of_coverage
{
    my($self) = @_;

    # Get value from mapstats
    my $depth = $self->mapstats->mean_target_coverage;
   
    # QC
    if($self->is_qc_mapstats && $self->is_mapping_complete)
    {
        my $genome_size        = $self->reference_size;
        my $rmdup_bases_mapped = $self->mapstats->rmdup_bases_mapped;
        my $qc_bases           = $self->mapstats->raw_bases;
        my $bases              = $self->bases;

        # if no mapstats value then calculate from mapped bases / genome size.
        $depth = ($genome_size ? $rmdup_bases_mapped/$genome_size : undef) unless defined $depth;

        # scale by lane bases / sample bases
        $depth = $depth * $bases / $qc_bases if defined $depth;
    }

    # Format and return
    $depth = sprintf("%.2f",$depth) if defined $depth;
    return $depth;
}

sub _build_depth_of_coverage_sd
{
    my($self) = @_;

    # Get value from mapstats
    my $depth_sd = $self->mapstats->target_coverage_sd;

    # QC
    if($self->is_qc_mapstats && $self->is_mapping_complete)
    {
        my $qc_bases           = $self->mapstats->raw_bases;
        my $bases              = $self->bases;

        # scale by lane bases / sample bases
        $depth_sd = $depth_sd * $bases / $qc_bases if defined $depth_sd;
    }

    # Format and return
    $depth_sd = sprintf("%.2f",$depth_sd) if defined $depth_sd;
    return $depth_sd;
}

sub _build_duplication_rate
{
    my($self) = @_;
    return undef unless $self->is_qc_mapstats;

    my $dupe_rate;
    if($self->is_mapping_complete)
    {
	my $rmdup_reads_mapped = $self->mapstats->rmdup_reads_mapped;
	my $reads_mapped       = $self->mapstats->reads_mapped;
	$dupe_rate = sprintf("%.4f", (1-$rmdup_reads_mapped/$reads_mapped));
    }
    return $dupe_rate;
}

sub _build_error_rate
{
    my($self) = @_;
    return undef unless $self->is_qc_mapstats;
    return $self->is_mapping_complete ? sprintf("%.3f",$self->mapstats->error_rate) : undef;
}

sub _build_npg_qc
{
    my($self) = @_;
    return $self->lane->npg_qc_status();
}

sub _build_manual_qc
{
    my($self) = @_;
    return $self->lane->qc_status();
}

# Return genome coverage percent greater then 1X, 5X etc.
sub _target_bases_X_perc
{
    my($self,$coverdepth) = @_;

    # Check to constrain $coverdepth
    my %allowed = (1 => 1, 2 => 1, 5 => 1, 10 => 1, 20 => 1, 50 => 1, 100 => 1);
    return undef unless exists $allowed{$coverdepth};

    # return percentage or undef
    my $coverdepth_field = 'target_bases_'.$coverdepth.'X';
    my $cover_perc = $self->mapstats->$coverdepth_field;
    $cover_perc = sprintf("%.1f",$cover_perc) if defined $cover_perc;
    return $cover_perc;
}


# Checks Project, Sample, Assembly and Mapper set.
sub is_all_tables_set
{
    my ($self) = @_;
    my @vrtrack_table_obj = qw(_vrtrack_project _vrtrack_sample _vrtrack_assembly _vrtrack_mapper);
    for my $table_obj (@vrtrack_table_obj)
    {
	return 0 unless defined $self->$table_obj;
    }
    return 1;
}

# Add selected QC mapstat values to Mapping mapstat values.
sub transfer_qc_values
{
    my($self, $qc_row) = @_;
    
    # error check
    return 0 unless defined $qc_row && $qc_row->is_qc_mapstats;
    return 0 unless defined $self->is_mapping_complete; # Do not update failed/running mapping.

    # update cells
    my @list_cells = qw(adapter_perc transposon_perc);

    for my $cell (@list_cells)
    {
	$self->$cell($qc_row->$cell) unless defined $self->$cell;
    }

    return 1;
}

1;
