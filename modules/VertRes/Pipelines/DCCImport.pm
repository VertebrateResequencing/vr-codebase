=head1 NAME

VertRes::Pipelines::DCCImport - pipeline for importing fastq files from the DCC

=head1 SYNOPSIS

# make the config files, which specifies the details for connecting to the
# VRTrack g1k-meta database and the data roots:
echo '__VRTrack_Import__ dcc_import.conf' > pipeline.config
# where dcc_import.conf contains:
root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::DCCImport',
prefix  => '_',

ftp_limit => 100,

db  => {
    database => 'g1k_meta',
    host     => 'mcs4a',
    port     => 3306,
    user     => 'vreseq_rw',
    password => 'xxxxxxx',
},

# (the ftp_limit option is the maximum number of lanes worth of fastqs to
#  download at once)

# run the pipeline:
run-pipeline -c pipeline.config -s 1

# (and make sure it keeps running by adding that last to a regular cron job)

=head1 DESCRIPTION

A module for importing new fastqs from the DCC: first it downloads the fastq
files, then it adds them to our data hierarchy for G1K. It depends on
si_to_vr.pl script to be run first, to update our database with the latest
sequence.index information, and afterwards the mapping pipeline should be run,
which will pick up on newly imported lanes and map them.

It will also work to import non-G1K fastqs in a local directory into a G1K-style
hierarchy, assuming a fake sequence.index has been used with si_to_vr.pl on
a custom database. In this case, the import.conf file must contain the
fastq_base option:

root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::DCCImport',
prefix  => '_',

db  => {
    database => 'special_meta',
    host     => 'mcs4a',
    port     => 3306,
    user     => 'vreseq_rw',
    password => 'xxxxxxx',
},

data => {
    fastq_base => '/abs/path/to/flattened/fastq/dir'
},

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::DCCImport;

use strict;
use warnings;
use File::Basename;
use VertRes::IO;
use VertRes::Parser::fastqcheck;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::Library;
use VRTrack::Sample;
use LSF;

use base qw(VertRes::Pipeline);

our $actions = [ { name     => 'import_fastqs',
                   action   => \&import_fastqs,
                   requires => \&import_fastqs_requires, 
                   provides => \&import_fastqs_provides },
                 { name     => 'update_db',
                   action   => \&update_db,
                   requires => \&update_db_requires, 
                   provides => \&update_db_provides },
                 { name     => 'cleanup',
                   action   => \&cleanup,
                   requires => \&cleanup_requires, 
                   provides => \&cleanup_provides } ];

our %options = (do_cleanup => 0,
                fastq_base => 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/',
                bsub_opts => '-q normal');

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::DCCImport->new();
 Function: Create a new VertRes::Pipelines::DCCImport object;
 Returns : VertRes::Pipelines::Mapping object
 Args    : lane => '/path/to/lane'
           do_cleanup => boolean (default false: don't do the cleanup action)
           files => \@files_for_import
           root => path to dir where hierarchy will be built to contain the
                   downloaded files (data dir)
           fastq_base => ftp url containing the data dir (default exists);
                         can also be a local path to the dir containing the
                         fastq files
           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    $self->throw("files option is required") unless defined $self->{files};
    
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
    my $self = shift;
    
    my @fastqs;
    for my $file (@{$self->{files}}) {
        my $fastq = basename($file);
        my $fqc = $fastq.'.fastqcheck';
        foreach my $name ($fastq, $fqc) {
            push(@fastqs, $name);
        }
    }
    
    return \@fastqs;
}

=head2 import_fastqs

 Title   : get_fastqs
 Usage   : $obj->import_fastqs('/path/to/lane', 'lock_filename');
 Function: Downloads fastq files from the DCC into the data hierarchy. Checks
           they're OK. Creates a fastqcheck file. Notes the lane's gender by
           adding a text file called 'MALE' or 'FEMALE'.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub import_fastqs {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $lane_obj = $self->{vrlane};
    
    # download the files that haven't been downloaded already
    my @files = @{$lane_obj->files || []};
    
    my $did_one = 0;
    for my $file_obj (@files) {
        my $file = $file_obj->name;
        my $ftp_path = $self->{fastq_base}.$file;
        my $fastq = $file_obj->hierarchy_name;
        my $data_path = $self->{io}->catfile($lane_path, $fastq);
        
        my $precheck_file = $data_path;
        $precheck_file =~ s/\.fastq\.gz$/.precheck.fastq.gz/;
        my $fqc_file = $data_path.'.fastqcheck';
        
        next if (-s $data_path && -s $fqc_file);
        $did_one++;
        
        my $md5 = $file_obj->md5;
        my $total_reads = $file_obj->raw_reads;
        my $total_bases = $file_obj->raw_bases;
        
        my $script_name = $self->{io}->catfile($lane_path, $self->{prefix}."import_$fastq.pl");
        
        open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
        print $scriptfh qq{
use strict;
use VertRes::IO;
use File::Spec;
use File::Copy;
use VertRes::Wrapper::fastqcheck;
use VertRes::Parser::fastqcheck;

unless (-s '$data_path') {
    unless (-s '$precheck_file') {
        my \$io = VertRes::IO->new;
        
        if ('$ftp_path' =~ qr{^(?:ftp|http)://}) {
            # download to a temp-named file
            my \$got = \$io->get_remote_file('$ftp_path',
                                             save => '$precheck_file',
                                             md5  => '$md5');
            
            unless (\$got eq '$precheck_file') {
                die "Failed to download '$ftp_path' to '$precheck_file'\n";
            }
        }
        elsif (-s '$ftp_path') {
            copy('$ftp_path', '$precheck_file') || die "Failed to move '$ftp_path' to '$precheck_file'\n";
            my \$ok = \$io->verify_md5('$precheck_file', '$md5');
            unless (\$ok) {
                unlink('$precheck_file');
                die "md5 check on '$ftp_path' after moving it to hierarchy failed! The hierarchy copy was deleted.\n";
            }
        }
        else {
            die "The given path to the fastq file to import was invalid ('$ftp_path')\n";
        }
    }
    
    if (-s '$precheck_file') {
        # run fastqcheck on it and confirm it matches expectations
        unless (-s '$fqc_file') {
            # make the fqc file
            my \$wrapper = VertRes::Wrapper::fastqcheck->new();
            \$wrapper->run('$precheck_file', '$fqc_file.temp');
            die "fastqcheck on $precheck_file failed - try again?" unless \$wrapper->run_status >= 1;
            die "fastqcheck failed to make the file $fqc_file.temp" unless -s '$fqc_file.temp';
            
            # check it is parsable
            my \$parser = VertRes::Parser::fastqcheck->new(file => '$fqc_file.temp');
            my \$num_sequences = \$parser->num_sequences();
            if (\$num_sequences) {
                move('$fqc_file.temp', '$fqc_file') || die "failed to rename '$fqc_file.temp' to '$fqc_file'";
            }
            else {
                die "fastqcheck file '$fqc_file.temp' was created, but doesn't seem valid\n";
            }
        }
        
        if (-s '$fqc_file') {
            # check against expectation
            my \$parser = VertRes::Parser::fastqcheck->new(file => '$fqc_file');
            my \$num_sequences = \$parser->num_sequences();
            my \$num_bases = \$parser->total_length();
            
            my \$ok = 0;
            if (\$num_sequences == $total_reads) {
                \$ok++;
            }
            else {
                warn "fastqcheck reports number of reads as \$num_sequences, not the expected $total_reads\n";
            }
            if (\$num_bases == $total_bases) {
                \$ok++;
            }
            else {
                warn "fastqcheck reports number of bases as \$num_bases, not the expected $total_bases\n";
            }
            
            if (\$ok == 2) {
                move('$precheck_file', '$data_path') || die "Could not rename '$precheck_file' to '$data_path'\n";
            }
            else {
                die "fastqcheck of fastq file doesn't match expectations\n";
            }
        }
    }
}

exit;
        };
        close $scriptfh;
        
        LSF::run($action_lock, $lane_path, $self->{prefix}.'import_'.$fastq, $self, qq{perl -w $script_name});
    }
    
    my $male_file = $self->{io}->catfile($lane_path, 'MALE');
    my $female_file = $self->{io}->catfile($lane_path, 'FEMALE');
    my $unknown_file = $self->{io}->catfile($lane_path, 'UNKNOWN');
    if (! -s $male_file && ! -s $female_file && ! -s $unknown_file) {
        # find the gender for this lane
        my $gender;
        my $vrtrack = $lane_obj->vrtrack;
        my $lib_id = $lane_obj->library_id;
        my $lib = VRTrack::Library->new($vrtrack, $lib_id);
        my $samp_id = $lib->sample_id;
        my $samp = VRTrack::Sample->new($vrtrack, $samp_id);
        my $individual = $samp->individual();
        my $sex = $individual->sex();
        
        # create an empty file in the lane indicating the gender
        my $fh;
        if ($sex eq 'F') {
            open($fh, '>', $female_file);
        }
        elsif ($sex eq 'M') {
            open($fh, '>', $male_file);
        }
        else {
            open($fh, '>', $unknown_file);
        }
        close($fh);
    }
    
    return $did_one ? $self->{No} : $self->{Yes};
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
    return $self->import_fastqs_provides($lane_path);
}

=head2 update_db_provides

 Title   : update_db_provides
 Usage   : my $provided_files = $obj->update_db_provides('/path/to/lane');
 Function: Find out what files the update_db action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub update_db_provides {
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
    
    my $vrlane  = $self->{vrlane};
    my $vrtrack = $vrlane->vrtrack;
    my @files = @{$vrlane->files || []};
    
    $vrtrack->transaction_start();
    
    my $imported_files = 0;
    my $total_files = 0;
    my $attempted_updates = 0;
    my $successful_updates = 0;
    my $largest_read_len = 0;
    for my $file_obj (@files) {
        my $file = $file_obj->name;
        my $fastq = $file_obj->hierarchy_name;
        my $data_path = $self->{io}->catfile($lane_path, $fastq);
        
        $total_files++;
        next unless -s $data_path;
        $imported_files++;
        
        $file_obj->is_processed('import', 1);
        
        # set mean_q and read_len, based on the fastqcheck info
        my $fqc = VertRes::Parser::fastqcheck->new(file => "$data_path.fastqcheck");
        my $read_len = int($fqc->avg_length);
        $largest_read_len = $read_len if $read_len > $largest_read_len;
        $file_obj->read_len($read_len);
        $file_obj->mean_q($fqc->avg_qual);
        
        my $ok = $file_obj->update();
        $attempted_updates++;
        $successful_updates++ if $ok;
    }
    
    # update import status of whole lane if we imported all the files
    if ($imported_files && $imported_files == $total_files) {
        $vrlane->is_processed('import', 1);
        
        if ($largest_read_len) {
            $vrlane->read_len($largest_read_len);
        }
        
        my $ok = $vrlane->update();
        $attempted_updates++;
        $successful_updates++ if $ok;
    }
    
    $vrtrack->transaction_commit();
    
    if ($attempted_updates == $successful_updates) {
        return $self->{Yes};
    }
    else {
        $self->throw("Some database updates for $lane_path failed");
    }
}

=head2 cleanup_requires

 Title   : cleanup_requires
 Usage   : my $required_files = $obj->cleanup_requires('/path/to/lane');
 Function: Find out what files the cleanup action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_requires {
    return [];
}

=head2 cleanup_provides

 Title   : cleanup_provides
 Usage   : my $provided_files = $obj->cleanup_provides('/path/to/lane');
 Function: Find out what files the cleanup action generates on success.
 Returns : n/a - doesn't provide anything, and always runs
 Args    : lane path

=cut

sub cleanup_provides {
    return [];
}

=head2 cleanup

 Title   : cleanup
 Usage   : $obj->cleanup('/path/to/lane', 'lock_filename');
 Function: Unlink all the pipeline-related files (_*) in a lane.
 Returns : $VertRes::Pipeline::Yes
 Args    : lane path, name of lock file to use

=cut

sub cleanup {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $prefix = $self->{prefix};
    
    my $fastqs = $self->import_fastqs_provides($lane_path);
    my @files = ('log', 'job_status');
    foreach my $fastq (@{$fastqs}) {
        next if $fastq =~ /fastqcheck/;
        $fastq = basename($fastq);
        foreach my $suffix ('o', 'e', 'pl') {
            push(@files, 'import_'.$fastq.'.'.$suffix);
        }
    }
    
    foreach my $file (@files) {
        unlink($self->{io}->catfile($lane_path, $prefix.$file));
    }
    
    $self->{io}->rmtree($self->{io}->catfile($lane_path, 'split_se'));
    $self->{io}->rmtree($self->{io}->catfile($lane_path, 'split_pe'));
    
    return $self->{Yes};
}

sub is_finished {
    my ($self, $lane_path, $action) = @_;
    
    if ($action->{name} eq 'update_db') {
        $self->{vrlane}->is_processed('import') ? return $self->{Yes} : return $self->{No};
    }
    elsif ($action->{name} eq 'cleanup') {
        return $self->{No};
    }
    
    return $self->SUPER::is_finished($lane_path, $action);
}

1;
