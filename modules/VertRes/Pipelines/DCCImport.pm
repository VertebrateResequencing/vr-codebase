=head1 NAME

VertRes::Pipelines::DCCImport - pipeline for importing fastq files fromt the DCC

=head1 SYNOPSIS

# make the config files, which specifies the details for connecting to the
# VRTrack g1k-meta database and the import/mapping roots:
echo '__VRTrack_Import__ dcc_import.conf' > pipeline.config
# where dcc_import.conf contains:
root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::DCCImport',
prefix  => '_',

db  => {
    database => 'g1k_meta',
    host     => 'mcs4a',
    port     => 3306,
    user     => 'vreseq_ro',
},

data => {
    db  => 
    {
        database => 'g1k_meta',
        host     => 'mcs4a',
        port     => 3306,
        user     => 'vreseq_rw',
        password => 'xxxxxxx',
    },

    mapping_root => '/abs/path/to/root/mapping/dir',
},

# run the pipeline:
run-pipeline -c pipeline.config -s 30

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for importing new fastqs from the DCC: first it downloads the fastq
files, then it adds them to our data&mapping hierarchies for G1K. It depends on
si_to_vr.pl script to be run first, to update our database with the latest
sequence.index information, and afterwards the mapping pipeline should be run,
which will pick up on newly imported lanes and map them.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::DCCImport;

use strict;
use warnings;
use File::Basename;
use VertRes::IO;
use LSF;

use base qw(VertRes::Pipeline);

our $actions = [ { name     => 'import_fastqs',
                   action   => \&import_fastqs,
                   requires => \&import_fastqs_requires, 
                   provides => \&import_fastqs_provides },
                 { name     => 'update_db',
                   action   => \&update_db,
                   requires => \&update_db_requires, 
                   provides => \&update_db_provides } ];

our %options = (do_cleanup => 0);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::DCCImport->new();
 Function: Create a new VertRes::Pipelines::DCCImport object;
 Returns : VertRes::Pipelines::Mapping object
 Args    : lane => '/path/to/lane'
           do_cleanup => boolean (default false: don't do the cleanup action)
           files => \@files_for_import
           root => path to dir to hold the fastq files (data dir)
           mapping_root => path to dir in which mapping will be done
           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    foreach (qw(files mapping_root)) {
        $self->throw("$_ option is required") unless defined $self->{$_};
    }
    
    $self->{io} = VertRes::IO->new;
    
    return $self;
}

=head2 import_fastqs_requires

 Title   : import_fastqs_requires
 Usage   : my $required_files = $obj->import_fastqs_requires('/path/to/lane');
 Function: Find out what files the import_fastqs action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub import_fastqs_requires {
    my $self = shift;
    return [];
}

=head2 import_fastqs_provides

 Title   : import_fastqs_provides
 Usage   : my $provided_files = $obj->import_fastqs_provides('/path/to/lane');
 Function: Find out what files the import_fastqs action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub import_fastqs_provides {
    my ($self, $lane_path) = @_;
    return $self->_fastq_files($lane_path);
}

sub _fastq_files {
    my ($self, $lane_path) = @_;
    
    my @fastqs;
    for my $file (@{$self->{files}}) {
        push(@fastqs, basename($file));
        unless ($self->{hierarchy_path}) {
            $self->throw("no hierarchy path");
        }
        print "mapping_root[$self->{mapping_root}], hierarchy_path[$self->{hierarchy_path}], basename[".basename($file)."]\n";;
        push(@fastqs, $self->{io}->catfile($self->{mapping_root}, $self->{hierarchy_path}, basename($file)));
    }
    
    return \@fastqs;
}

=head2 import_fastqs

 Title   : get_fastqs
 Usage   : $obj->import_fastqs('/path/to/lane', 'lock_filename');
 Function: Downloads fastq files from the DCC into the data hierarchy. Checks
           they're OK. Creates a fastqcheck file. Creates symlinks to these in
           the mapping hierarchy.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub import_fastqs {
    my ($self, $lane_path, $action_lock) = @_;
    
    print "import_fastqs got lane_path $lane_path\n";
    for my $file (@{$self->{files}}) {
        print "  got a file $file\n"
    }
    
    exit;
    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
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
    return $self->_fastq_files($lane_path);
}

=head2 update_db_provides

 Title   : update_db_provides
 Usage   : my $provided_files = $obj->update_db_provides('/path/to/lane');
 Function: Find out what files the update_db action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub update_db_provides {
    my ($self, $lane_path) = @_;
    # run this action if we have the database details to write to the db
    return 0 if exists($self->{db});
    # otherwise, the empty list returned means that the action thinks it is
    # complete
    return [];
}

=head2 update_db

 Title   : update_db
 Usage   : $obj->update_db('/path/to/lane', 'lock_filename');
 Function: Update the g1k-meta database to indicate which fastqs have now been
           successfully downloaded and imported, and are ready to be mapped.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub update_db {
    my ($self, $lane_path, $action_lock) = @_;
    
    
    
    return $self->{Yes};
}

1;
