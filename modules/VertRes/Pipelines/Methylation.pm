=head1 NAME

VertRes::Pipelines::Methylation - Pipeline for methylation analysis on pacbio assembly

=head1 SYNOPSIS

# make the config files, which specifies the details for connecting to the
# VRTrack g1k-meta database and the data roots:
# #echo '__VRTrack_Methylation__ methylation.conf' > pipeline.config
# where methlyation.conf contains:
root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::Methylation',
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
    reference_fasta => '/foo/bar/assembly.fasta'  # if not set or is empty string, will use pipeline assembly fasta
    min_ipdratio => 0,
}


# __VRTrack_PacbioAssembly__

# run the pipeline:
run-pipeline -c pipeline.config

=head1 DESCRIPTION

Pipeline for assembling pacbio data

=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Pipelines::Methylation;

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
use Bio::PacbioMethylation::RSModificationRunner;

our $actions = [
    {
        name     => 'pacbio_methylation',
        action   => \&pacbio_methylation,
        requires => \&pacbio_methylation_requires,
        provides => \&pacbio_methylation_provides
    },
    {
        name     => 'update_db',
        action   => \&update_db,
        requires => \&update_db_requires,
        provides => \&update_db_provides
    }
];

our %options = (
    bsub_opts => '',
);

sub new {
    my ( $class, @args ) = @_;

    my $self = $class->SUPER::new( %options, actions => $actions, @args );
    if ( defined( $self->{db} ) ) {
        $self->{vrtrack} = VRTrack::VRTrack->new( $self->{db} ) or $self->throw("Could not connect to the database\n");
    }
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    unless (defined $self->{reference_fasta} and $self->{reference_fasta} ne '') {
        $self->{reference_fasta} = File::Spec->catfile($self->{lane_path}, 'pacbio_assembly', 'contigs.fa');
    }

    return $self;
}

sub pacbio_methylation_provides {
    my ($self) = @_;
    return [$self->{lane_path}."/".$self->{prefix}."methylation_done"];
}


sub pacbio_methylation_requires {
    my ($self) = @_;
    my $file_regex = $self->{lane_path}."/".'*.bax.h5';
    my @files = glob( $file_regex);
    die 'no .bax.h5 files found' if(@files == 0);
    return \@files;
}


=head2 pacbio_methylation

Title   : pacbio_methylation
Usage   : $obj->pacbio('/path/to/lane', 'lock_filename');
Function: Run Pacbio methylation pipeline
Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
Args    : lane path, name of lock file to use

=cut

sub pacbio_methylation {
    my ($self, $lane_path, $action_lock) = @_;
    my $memory_in_mb = $self->{memory} || 15000;
    my $threads = $self->{threads} || 8;
    my $output_dir= $self->{lane_path}."/pacbio_methylation";
    my $bax_h5_files = join(' ', @{$self->pacbio_methylation_requires()});
    my $queue = $self->{queue}|| "normal";
    my $pipeline_version = $self->{pipeline_version} || '7.0';
    my $umask    = $self->umask_str;
    my $lane_name = $self->{vrlane}->name;
    my $min_ipdratio = $self->{min_ipdratio} || 0;
    my $done_file = File::Spec->catfile($self->{lane_path}, $self->{prefix} . "methylation_done");
    my $tmp_directory = File::Spec->catfile($self->{tmp_directory}, $self->{prefix}. 'methylation_'.$self->{vrlane}->name()) || $output_dir;

    my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."methylation.pl");
    open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
    print $scriptfh qq{
use strict;
use File::Touch;
use Bio::PacbioMethylation::RSModificationRunner;
$umask

if (-e qw[$tmp_directory]) {
    system("rm -rf $tmp_directory") and die \$!;
}

my \@files = qw~ $bax_h5_files ~;

my \$obj = Bio::PacbioMethylation::RSModificationRunner->new(
   bax_h5_files => \\\@files,
   reference_fasta => qw[$self->{reference_fasta}],
   outdir => qw[$tmp_directory],
   threads => $threads,
   clean => 1,
   min_ipdratio => $min_ipdratio,
);
\$obj->run();
};

    if ($tmp_directory ne $output_dir) {
        print $scriptfh qq{system("rm -fr $output_dir") and die \$!;
mkdir qw[$output_dir] or die \$!;
system("rsync -a $tmp_directory/ $output_dir") and die \$!;
system("rm -fr $tmp_directory") and die \$!;
};
    }

    print $scriptfh qq{
my \$done_file = qw[$done_file];
touch(\$done_file) or die "Error touch \$done_file";
};

    close $scriptfh;

    my $job_name = $self->{prefix}.'methylation';

    $self->delete_bsub_files($lane_path, $job_name);
    VertRes::LSF::run($action_lock, $lane_path, $job_name, {bsub_opts => "-n$threads -q $queue -M${memory_in_mb} -R 'select[mem>$memory_in_mb] rusage[mem=$memory_in_mb] span[hosts=1]'"}, qq{perl -w $script_name});

    return $self->{No};
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
    return $self->pacbio_methylation_provides();
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
    return [ $self->{lane_path}."/".$self->{prefix}."methylation_update_db_done"];
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

    return $$self{'Yes'} if $vrlane->is_processed('methylation');

    unless($vrlane->is_processed('methylation')){
        $vrtrack->transaction_start();
        $vrlane->is_processed('methylation',1);
        $vrlane->update() || $self->throw("Unable to set methylation status on lane $lane_path");
        $vrtrack->transaction_commit();
    }


    my $job_status =  File::Spec->catfile($lane_path, $self->{prefix} . 'job_status');
    Utils::CMD("rm $job_status") if (-e $job_status);


    my $prefix = $self->{prefix};
    # remove job files
    foreach my $file (qw(pacbio_methylation ))
    {
        foreach my $suffix (qw(o e pl))
        {
            unlink($self->{fsu}->catfile($lane_path, $prefix.$file.'.'.$suffix));
        }
    }

    Utils::CMD("touch ".$self->{fsu}->catfile($lane_path,$self->{prefix}."methylation_update_db_done")   );
    $self->update_file_permissions($lane_path);
    return $$self{'Yes'};
}


1;
