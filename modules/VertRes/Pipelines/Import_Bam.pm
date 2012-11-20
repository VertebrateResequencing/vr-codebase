package VertRes::Pipelines::Import_Bam;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use LSF;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VertRes::Utils::FileSystem;
use VertRes::Parser::bamcheck;
use VertRes::Wrapper::samtools;
use VertRes::Pipelines::Import_iRODS;

our @actions =
(
    # Creates the hierarchy path, downloads, gzips and checkbam the bam files.
    {
        'name'     => 'get_bams',
        'action'   => \&get_bams,
        'requires' => \&get_bams_requires, 
        'provides' => \&get_bams_provides,
    },

    # If all files were downloaded OK, update the VRTrack database.
    {
        'name'     => 'update_db',
        'action'   => \&VertRes::Pipelines::Import_iRODS::update_db,
        'requires' => \&update_db_requires, 
        'provides' => \&update_db_provides,
    },
);

our $options = 
{
    'bamcheck'        => 'bamcheck -q 20',
    'bsub_opts'       => "-q normal -R 'select[type==X86_64] rusage[thouio=1]'",
    'local_bam_dir'   => '',
};


# --------- OO stuff --------------

=head2 new

        Example    : my $qc = VertRes::Pipelines::TrackDummy->new( files=>[2451_1.bam] );
        Options    : See Pipeline.pm for general options.

                    bamcheck        .. 
                    files           .. Array reference to the list of files to be imported.

=cut

sub new 
{
    my ($class, %args) = @_;
    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,%args);
    $self->write_logs(1);
    if ( !$$self{files} ) { $self->throw("Missing the option files.\n"); }
    return $self;
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

#---------- get_bams ---------------------

# Requires nothing
sub get_bams_requires
{
    my ($self) = @_;
    return [];
}

sub get_bams_provides
{
    my ($self, $lane_path) = @_;
    my @provides;
    foreach my $file (@{$$self{files}}) {
        push @provides, File::Basename::basename($file);
    }
    @provides || $self->throw("Something went wrong; we don't seem to provide any bams!");
    return \@provides;
}

sub get_bams
{
    my ($self,$lane_path,$lock_file) = @_;

    my $opts = $self->dump_opts(qw(files bamcheck local_bam_dir db lane));

    my $prefix   = $$self{prefix};
    my $work_dir = $lane_path;

    # Create a script to be run on LSF.
    open(my $fh,'>', "$work_dir/${prefix}import_bams.pl") or $self->throw("$work_dir/${prefix}import_bams.pl: $!");
    print $fh qq[
use strict;
use warnings;
use VertRes::Pipelines::Import_Bam;

my $opts

my \$import = VertRes::Pipelines::Import_Bam->new(%\$opts);
\$import->get_files();

];

    close($fh);
    LSF::run($lock_file,$work_dir,"${prefix}import_bams",$self,qq[perl -w ${prefix}import_bams.pl]);

    return $$self{No};
}

sub get_files
{
    my ($self) = @_;
    
    if ( !$$self{db} ) { $self->throw("Expected the db key.\n"); }
    
    my $vrtrack = VRTrack::VRTrack->new($$self{db}) or $self->throw("Could not connect to the database\n");
    my $vrlane  = VRTrack::Lane->new_by_name($vrtrack,$$self{lane}) or $self->throw("No such lane in the DB: [$$self{lane}]\n");
    
    my $fsu = VertRes::Utils::FileSystem->new();
    my $samtools = VertRes::Wrapper::samtools->new();
    
    # Get files and run bamcheck on them
    for my $file (@{$$self{files}})
    {
        my $file_path = $file;
        if ($$self{local_bam_dir}) {
            $file_path = $fsu->catfile($$self{local_bam_dir}, $file);
        }
        unless ( -s $file_path ) { $self->throw("No such file on disk? [$file_path]\n"); next; }
        
        if ( !($file=~m{([^/]+)$}) ) { $self->throw("FIXME: [$file]"); }
        my $outfile = $1;
        if ( -s $outfile ) { next; }
        
        # Copy the BAM file over.
        $fsu->copy($file_path, "$outfile.tmp");
        $self->throw("FIXME: could not copy file $file_path to $outfile.tmp") unless (-s "$outfile.tmp");
        chmod 0664,"$outfile.tmp";
        
        # Get the md5sum and check
        my $vrfile = $vrlane->get_file_by_name($file);
        if ( !$vrfile ) { $self->throw("FIXME: the file not in the DB? [$file]"); }
        my $md5 = $vrfile->md5();
        if ( !$md5 ) { $self->throw("FIXME: the md5 not in the DB? [$file]"); }
        $self->throw("Could not verify md5 for $outfile.tmp") unless $fsu->verify_md5("$outfile.tmp", $md5);
        
        # sort the bam by coordinates and regenerate the md5
        $samtools->sort("$outfile.tmp", "sorted", m => 50000000);
        rename("$outfile.tmp.sorted.bam","$outfile.tmp") or $self->throw("rename $outfile.tmp.sorted.bam $outfile.tmp: $!");
        $md5 = $fsu->calculate_md5("$outfile.tmp");
        
        # Create the checksum file
        open(my $fh,'>',"$outfile.md5") or $self->throw("$outfile.md5: $!");
        print $fh "$md5  $outfile\n";
        close($fh);
        
        # Index the BAM
        $samtools->run_method('system');
        $samtools->index("$outfile.tmp", "$outfile.tmp.bai");
        $samtools->run_status >= 1 || $self->throw("Failed to create $outfile.tmp.bai");
        rename("$outfile.tmp.bai","$outfile.bai") or $self->throw("rename $outfile.tmp.bai $outfile.bai: $!");
        
        # Run bamcheck
        Utils::CMD(qq[$$self{bamcheck} $outfile.tmp > $outfile.tmp.bc]);
        rename("$outfile.tmp.bc","$outfile.bc") or $self->throw("rename $outfile.tmp.bc $outfile.bc: $!");
        rename("$outfile.tmp",$outfile) or $self->throw("rename $outfile.tmp $outfile: $!");
    }
}


#---------- update_db ---------------------

sub update_db_requires
{
    my ($self,$lane_path) = @_;
    my @requires = @{$self->get_bams_provides($lane_path)};
    @requires || $self->throw("Something went wrong; we don't seem to require any bams!");
    return \@requires;
}

# This subroutine will check existence of the key 'db'. If present, it is assumed
#   that Import should write the stats and status into the VRTrack database. In this
#   case, 0 is returned, meaning that the task must be run. The task will change the
#   QC status from NULL to pending, therefore we will not be called again.
#
#   If the key 'db' is absent, the empty list is returned and the database will not
#   be written.
#
sub update_db_provides
{
    my ($self) = @_;
    if ( exists($$self{db}) ) { return 0; }
    my @provides = ();
    return \@provides;
}

#---------- Debugging and error reporting -----------------

sub format_msg
{
    my ($self,@msg) = @_;
    return '['. scalar gmtime() ."]\t". join('',@msg);
}

sub warn
{
    my ($self,@msg) = @_;
    my $msg = $self->format_msg(@msg);
    if ($self->verbose > 0) 
    {
        print STDERR $msg;
    }
    $self->log($msg);
}

sub debug
{
    # The granularity of verbose messaging does not make much sense
    #   now, because verbose cannot be bigger than 1 (made Base.pm
    #   throw on warn's).
    my ($self,@msg) = @_;
    if ($self->verbose > 0) 
    {
        my $msg = $self->format_msg(@msg);
        print STDERR $msg;
        $self->log($msg);
    }
}

sub throw
{
    my ($self,@msg) = @_;
    my $msg = $self->format_msg(@msg);
    Utils::error($msg);
}

sub log
{
    my ($self,@msg) = @_;

    my $msg = $self->format_msg(@msg);
    my $status  = open(my $fh,'>>',$self->log_file);
    if ( !$status ) 
    {
        print STDERR $msg;
    }
    else 
    { 
        print $fh $msg; 
    }
    if ( $fh ) { close($fh); }
}


1;
