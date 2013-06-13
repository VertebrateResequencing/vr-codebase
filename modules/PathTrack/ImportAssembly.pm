# To run this pipeline step:
# run-pipeline -c pipeline-test.conf -o -v -v

package PathTrack::ImportAssembly;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use VRTrack::VRTrack;
use Utils;
use VertRes::LSF;
use VertRes::Utils::FileSystem;

our @actions =
(
    {
        'name'     => 'import_assembly',
        'action'   => \&import_assembly,
        'requires' => \&import_assembly_requires, 
        'provides' => \&import_assembly_provides,
    },
    {
        'name'     => 'create_symlink',
        'action'   => \&create_symlink,
        'requires' => \&create_symlink_requires, 
        'provides' => \&create_symlink_provides,
    },
    {
        'name'     => 'process_assembly',
        'action'   => \&process_assembly,
        'requires' => \&process_assembly_requires, 
        'provides' => \&process_assembly_provides,
    },
);

our $options = 
{
    # Executables
    samtools_exec => 'samtools',
    refstats_exec => 'ref-stats',
    bwa_exec => 'bwa',
    stats_exec => '/nfs/users/nfs_a/ap12/genlibpy/genepy/pathtrack/stats.py',

    'bsub_opts'       => "-q normal",
};


### ---------------------------------------------------------------------------
### new
### ---------------------------------------------------------------------------
sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,@args);
    
    # check required options are provided
    $self->throw("Missing root option in config.\n") unless $self->{root};
    $self->throw("Missing db option in config") unless $self->{db};
    $self->throw("Missing vrlane object.\n") unless $self->{vrlane};
    $self->throw("Missing assembly_id option in config.\n") unless $self->{assembly_id};
    $self->throw("Missing contigs option in config.\n") unless $self->{contigs};
    $self->throw("Missing scaffolds option in config.\n") unless $self->{scaffolds};
    $self->throw("Missing samtools_exec option in config.\n") unless $self->{samtools_exec};

    # database connection
    my $vrtrack = VRTrack::VRTrack->new($self->{db}) or $self->throw("Could not connect to the database\n");

    # set assembly path
    my $assembly_path = $vrtrack->hierarchy_path_of_lane($self->{vrlane}, "genus:species-subspecies:ASSEMBLY:assemblyid");
    $assembly_path =~ s/assemblyid/$self->{assembly_id}/g;
    $self->{assembly_path} = $self->{root} . '/' . $assembly_path;

    # set contigs and scaffolds default filenames
    $self->{contigs_filename} = "Contigs.fna"; 
    $self->{scaffolds_filename} = "Scaffolds.fna";

    # set fastqs symlink
    $self->{fastqs_symlink} = $self->{assembly_path} . '/fastqs/' . $self->{lane};

    # temp debug print statements
    #print "lane          : ", $self->{lane}, "\n";
    #print "lane_path     : ", $self->{lane_path}, "\n";
    #print "fastqs_symlink: ", $self->{fastqs_symlink}, "\n";
    #print "assembly_id   : ", $self->{assembly_id}, "\n";
    #print "contigs       : ", $self->{contigs}, "\n";
    #print "scaffolds     : ", $self->{scaffolds}, "\n";
    #print "root          : ", $self->{root}, "\n";
    #print "assembly path : ", $self->{assembly_path}, "\n";

    return $self;
}


### ---------------------------------------------------------------------------
### import_assembly
### ---------------------------------------------------------------------------
sub import_assembly_requires
{
    my ($self) = @_;
    my @requires;
    push(@requires, "$$self{contigs}");
    push(@requires, "$$self{scaffolds}");
    return \@requires;
}

sub import_assembly_provides
{
    my ($self) = @_;
    my @provides;
    push(@provides, "$$self{assembly_path}/$$self{contigs_filename}");
    push(@provides, "$$self{assembly_path}/$$self{scaffolds_filename}");
    return \@provides;
}

sub import_assembly
{
    my ($self, $lane_path, $action_lock) = @_;
    my $fsu = VertRes::Utils::FileSystem->new();

    # create assembly hierarchy on disk
    Utils::create_dir($self->{assembly_path});

    # import contigs and scaffolds files
    $fsu->copy($self->{contigs}, $$self{assembly_path} . "/" . $$self{contigs_filename});
    $fsu->copy($self->{scaffolds}, $$self{assembly_path} . "/" . $$self{scaffolds_filename});
    
    return $$self{'Yes'};
}

### ---------------------------------------------------------------------------
### create_symlink
### ---------------------------------------------------------------------------
sub create_symlink_requires
{
    my ($self) = @_;
    my @requires;
    push(@requires, "$$self{assembly_path}");
    return \@requires;
}

sub create_symlink_provides
{
    my ($self) = @_;
    my @provides;
    push(@provides, "$$self{fastqs_symlink}");
    return \@provides;
}

sub create_symlink
{
    my ($self, $lane_path, $action_lock) = @_;

    # create fastqs dir
    my $fastqs_dir = $self->{assembly_path} . '/fastqs/';
    Utils::create_dir($fastqs_dir);

    # create symlink to fastqs files
    Utils::relative_symlink($$self{lane_path}, $$self{fastqs_symlink}); 
    
    return $$self{'Yes'};
}

### ---------------------------------------------------------------------------
### process_assembly
### ---------------------------------------------------------------------------
sub process_assembly_requires
{
    my ($self) = @_;
    my @requires;
    push(@requires, "$$self{assembly_path}/$$self{contigs_filename}");
    push(@requires, "$$self{assembly_path}/$$self{scaffolds_filename}");
    return \@requires;
}

sub process_assembly_provides
{
    my ($self) = @_;
    my @provides;
    push(@provides, "$$self{assembly_path}/_stats.pl");
    return \@provides;
}

sub process_assembly
{
    my ($self, $lane_path, $lock_file) = @_;

    # dynamic script to be run by LSF
    open(my $fh,'>', "$$self{assembly_path}/_stats.pl") or Utils::error("$$self{assembly_path}/_stats.pl: $!");
    print $fh
qq[
use strict;
use warnings;
use Utils;

# run samtools faidx
Utils::CMD("$$self{samtools_exec} faidx $$self{assembly_path}/$$self{scaffolds_filename}");
if ( ! -s "$$self{assembly_path}/$$self{scaffolds_filename}.fai" ) { 
    Utils::error("The command ended with an error:\\n\\t$$self{samtools_exec} faidx $$self{assembly_path}/$$self{scaffolds_filename}\\n");
}

# run ref-stats
Utils::CMD("$$self{refstats_exec} -r $$self{assembly_path}/$$self{scaffolds_filename} > $$self{assembly_path}/$$self{scaffolds_filename}.refstats");
if ( ! -s "$$self{assembly_path}/$$self{scaffolds_filename}.refstats" ) { 
    Utils::error("The command ended with an error:\\n\\t$$self{refstats_exec} -r $$self{assembly_path}/$$self{scaffolds_filename} > $$self{assembly_path}/$$self{scaffolds_filename}.refstats\\n");
}

# run bwa
Utils::CMD("$$self{bwa_exec} index $$self{assembly_path}/$$self{scaffolds_filename}");
if ( ! -s "$$self{assembly_path}/$$self{scaffolds_filename}.bwt" ) { 
    Utils::error("The command ended with an error:\\n\\t$$self{bwa_exec} index $$self{assembly_path}/$$self{scaffolds_filename}\\n");
}

# run stats
Utils::CMD("python $$self{stats_exec} -f $$self{assembly_path}/$$self{contigs_filename}");
if ( ! -s "$$self{assembly_path}/$$self{contigs_filename}.stats" ) { 
    Utils::error("The command ended with an error:\\n\\tpython $$self{stats_exec} -f $$self{assembly_path}/$$self{contigs_filename}\\n");
}
Utils::CMD("python $$self{stats_exec} -f $$self{assembly_path}/$$self{scaffolds_filename}");
if ( ! -s "$$self{assembly_path}/$$self{scaffolds_filename}.stats" ) { 
    Utils::error("The command ended with an error:\\n\\tpython $$self{stats_exec} -f $$self{assembly_path}/$$self{scaffolds_filename}\\n");
}
];
    close($fh);
    VertRes::LSF::run($lock_file, $$self{assembly_path}, "_$$self{assembly_id}_$$self{lane}", $self, qq[perl -w _stats.pl]);

    return $$self{'No'};
}

1;

