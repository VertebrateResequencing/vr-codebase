=head1 NAME

VertRes::Pipelines::HierarchyUpdate - pipeline for importing fastq files and
                                      doing lane swaps etc.

=head1 SYNOPSIS

# make the config files, which specifies the details for connecting to the
# VRTrack g1k-meta database and the data roots:
echo '__VRTrack_HierarchyUpdate__ HierarchyUpdate.conf' > pipeline.config
# where HierarchyUpdate.conf contains:
root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::HierarchyUpdate',
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

A module for importing new or altered fastqs from some external source like the
DCC: first it downloads the fastq files, then it adds them to a data hierarchy.
It depends on update_vrmeta.pl script to be run first, to update the database
with the latest sequence.index information, and afterwards the mapping pipeline
should be run, which will pick up on newly imported/changed lanes and map them.

If the sequence.index had swaps (lanes moved to a different library or sample
etc.), and you allowed those changes to be made to the database when running
update_vrmeta.pl, the swaps will be made on disc, so that the disc hierarchy
reflects the database. If bams had already been made, their headers will be
corrected to reflect the new sample/library etc. Likewase .bas files will be
fixed.


It will also work to import non-DCC fastqs in a local directory into a G1K-style
hierarchy, assuming a fake sequence.index has been used with update_vrmeta.pl on
a custom database. In this case, the HierarchyUpdate.conf file must contain the
fastq_base option. If you don't have md5, number of reads and number of bases
information on the fastqs, you'll also need to supply the
calculate_missing_information option:

root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::HierarchyUpdate',
prefix  => '_',

db  => {
    database => 'special_meta',
    host     => 'mcs4a',
    port     => 3306,
    user     => 'vreseq_rw',
    password => 'xxxxxxx',
},

data => {
    fastq_base => '/abs/path/to/flattened/fastq/dir',
    calculate_missing_information => 1,
},

=head1 NOTES

It is the run-pipeline script that creates the lane directories for us. Hence
no make_path() or similar usage here.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::HierarchyUpdate;

use strict;
use warnings;
use File::Basename;
use File::Spec;
use File::Copy;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use VertRes::Parser::fastqcheck;
use VertRes::Utils::Hierarchy;
use VertRes::Utils::Sam;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::Library;
use VRTrack::Sample;
use LSF;

use base qw(VertRes::Pipeline);

our $actions = [ { name     => 'fix_swaps',
                   action   => \&fix_swaps,
                   requires => \&action_requires, 
                   provides => \&action_provides },
                 { name     => 'deletion',
                   action   => \&deletion,
                   requires => \&action_requires, 
                   provides => \&action_provides },
                 { name     => 'import_fastqs',
                   action   => \&import_fastqs,
                   requires => \&action_requires, 
                   provides => \&import_fastqs_provides },
                 { name     => 'update_db',
                   action   => \&update_db,
                   requires => \&update_db_requires, 
                   provides => \&action_provides },
                 { name     => 'cleanup',
                   action   => \&cleanup,
                   requires => \&action_requires, 
                   provides => \&action_provides } ];

our %options = (fastq_base => 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp',
                bsub_opts => '-q normal -R \'rusage[tmp=500]\'',
                calculate_missing_information => 0);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::HierarchyUpdate->new();
 Function: Create a new VertRes::Pipelines::HierarchyUpdate object.
 Returns : VertRes::Pipelines::HierarchyUpdate object
 Args    : lane => 'readgroup_id'
           files => \@files_for_import
           root => path to dir where hierarchy will be built to contain the
                   downloaded files (data dir)
           fastq_base => ftp url containing the data dir (DCC default exists);
                         can also be a local path to the dir containing the
                         fastq files
           calculate_missing_information => boolean (if this is true and the
                         database does not contain md5, reads or bases info on
                         a fastq, then these will be calculated and stored in
                         the db; otherwise it will fail with an error)
           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    $self->throw("files option is required") unless defined $self->{files};
    
    # if we've been supplied a list of lane paths to work with, instead of
    # getting the lanes from the db, we won't have a vrlane object; make one
    if (! $self->{vrlane}) {
        $self->throw("db option was not supplied in config") unless $self->{db};
        my $vrtrack = VRTrack::VRTrack->new($self->{db}) or $self->throw("Could not connect to the database\n");
        my $lane = $self->{lane} || $self->throw("lane readgroup not supplied, can't continue");
        my $vrlane  = VRTrack::Lane->new_by_name($vrtrack, $lane) or $self->throw("No such lane in the DB: [$lane]");
        $self->{vrlane} = $vrlane;
    }
    $self->{vrlane} || $self->throw("vrlane object missing");
    
    $self->{io} = VertRes::IO->new;
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    
    return $self;
}

=head2 action_requires

 Title   : action_requires
 Usage   : my $required_files = $obj->action_requires('/path/to/lane');
 Function: Find out what files the various actions need before they will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub action_requires {
    my $self = shift;
    return [];
}

=head2 action_provides

 Title   : action_provides
 Usage   : my $provided_files = $obj->action_provides('/path/to/lane');
 Function: Find out what files the various actions generate on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub action_provides {
    return [];
}

=head2 fix_swaps

 Title   : fix_swaps
 Usage   : $obj->fix_swaps('/path/to/lane', 'lock_filename');
 Function: Move a lane directory to the correct part of the hierarchy if there
           has been a swap. Corrects any bams and bas files that may have been
           created.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub fix_swaps {
    my ($self, $lane_path, $action_lock) = @_;
    
    # work out the pre-swap hierarchy path of the lane; our new lane location
    # already exists and we'll need to move files and correct headers
    my $vrlane = $self->{vrlane};
    my %info = VertRes::Utils::Hierarchy->new->lane_info($vrlane, pre_swap => 1);
    
    my $old_path = $info{hierarchy_path};
    if ($old_path ne $lane_path && -e $old_path) {
        # move the contents of lane_path (just now created by run-pipline,
        # containing pipline-related files like the lock file) to old_path, then
        # move old_path to the new location
        $self->{fsu}->move($lane_path, $old_path) || $self->throw("Failed to move contents of '$lane_path' to '$old_path'");
        move($old_path, $lane_path) || $self->throw("Failed to move '$old_path' to '$lane_path'");
        
        $self->_prune_path($old_path);
    }
    
    if ($vrlane->is_processed('mapped')) {
        # look for bam and bas files to correct
        opendir(my $dfh, $lane_path) || $self->throw("Could not open dir '$lane_path'");
        my (@bams, @bass);
        foreach (readdir($dfh)) {
            if (/\.bam$/) {
                push(@bams, $self->{fsu}->catfile($lane_path, $_));
            }
            elsif (/\.bam\.bas$/) {
                push(@bass, $self->{fsu}->catfile($lane_path, $_));
            }
        }
        
        unless (@bams && @bams == @bass) {
            $self->throw("$lane_path was supposed to be mapped, but either there were no bam files, or the number of bas files didn't match the number of bams");
        }
        
        my $su = VertRes::Utils::Sam->new();
        
        my %new_meta = ($info{lane} => { sample_name => $info{sample},
                                         library => $info{library_true},
                                         platform => $info{technology},
                                         centre => $info{centre},
                                         insert_size => $info{insert_size},
                                         project => $info{project} });
        
        foreach my $bam (@bams) {
            $su->rewrite_bam_header($bam, %new_meta) || $self->throw("Failed to correct the header of '$bam'");
        }
        
        foreach my $bas (@bass) {
            $su->rewrite_bas_meta($bas, %new_meta) || $self->throw("Failed to correct meta information in '$bam'");
        }
    }
    
    $vrlane->is_processed('swapped', 0);
    $vrlane->update || $self->throw("Having done the swap for $lane_path, failed to update the database");
    
    return $self->{Yes};
}

# if parent directory is now empty, delete the parent directory
sub _prune_path {
    my ($self, $path) = @_;
    
    my @dirs = File::Spec->splitdir($path);
    my $lane = pop(@dirs);
    $lane ||= pop(@dirs);
    $lane || $self->throw("Couldn't get lane directory off path '$path'");
    
    my $parent = File::Spec->catdir(@dirs);
    opendir(my $dfh, $parent) || $self->throw("Could not open dir '$parent'");
    my $empty = 1;
    foreach (readdir($dfh)) {
        next if /^\.{1,2}$/;
        $empty = 0;
        last;
    }
    
    if ($empty) {
        unlink($parent);
        $self->_prune_path($parent);
    }
}

=head2 deletion

 Title   : deletion
 Usage   : $obj->deletion('/path/to/lane', 'lock_filename');
 Function: If a lane was deleted entirely from the database, deletes it from
           disc. Liewise if fastqs were altered.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub deletion {
    my ($self, $lane_path, $action_lock) = @_;
    
    $self->throw("deletion not yet implemented");
    
    return $self->{Yes};
}

# copied over from update_vrmeta.pl; to become deletion() above...
#sub delete_lane {
#    my $lane = shift;
#    my $lane_name = $lane->hierarchy_name;
#    my $lane_path = $vrtrack->hierarchy_path_of_lane_hname($lane_name);
#    
#    if ($lane->is_processed('import')) {
#        # in case the lane has been moved and $lane_path is only a symlink, we
#        # get the absolute path
#        my $abs_lane_path = abs_path($lane_path);
#        my $deleted = rmtree($abs_lane_path);
#        rmtree($lane_path);
#        unlink($lane_path);
#        
#        unless ($deleted) {
#            issue_warning("Lane $lane_path was supposed to be removed from disc, but the attempt failed to delete anything; I'm still going to delete the lane from the database though");
#        }
#        else {
#            issue_warning("Lane $lane_path was removed from disc");
#        }
#    }
#    
#    my $deleted = $lane->delete;
#    if ($deleted) {
#        issue_warning("Lane $lane_name was deleted from the database");
#        return 1;
#    }
#    else {
#        issue_warning("Lane $lane_name was supposed to be deleted from the database, but the attempt failed!");
#        return 0;
#    }
#}

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
    unless (defined $self->{db}) {
        $self->{db} = $lane_obj->vrtrack->database_params;
    }
    
    # download the files that haven't been downloaded already
    my @files = @{$lane_obj->files || []};
    
    my $did_one = 0;
    for my $file_obj (@files) {
        my $file = $file_obj->name;
        my $ftp_path = $self->{fastq_base}.'/'.$file;
        my $fastq = $file_obj->hierarchy_name;
        my $data_path = $self->{fsu}->catfile($lane_path, $fastq);
        
        my $precheck_file = $data_path;
        $precheck_file =~ s/\.fastq\.gz$/.precheck.fastq.gz/;
        my $fqc_file = $data_path.'.fastqcheck';
        
        next if (-s $data_path && -s $fqc_file);
        $did_one++;
        
        my $md5 = $file_obj->md5 || '';
        my $total_reads = $file_obj->raw_reads || 0;
        my $total_bases = $file_obj->raw_bases || 0;
        
        my $file_id = $file_obj->id;
        unless ($self->{calculate_missing_information}) {
            $md5 || $self->throw("missing md5 for file $file");
            $total_reads || $self->throw("missing raw_reads for file $file");
            $total_bases || $self->throw("missing raw_bases for file $file");
        }
        
        my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."import_$fastq.pl");
        
        open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
        print $scriptfh qq{
use strict;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use File::Spec;
use File::Copy;
use VertRes::Wrapper::fastqcheck;
use VertRes::Parser::fastqcheck;

my \$data_path = '$data_path';
my \$precheck_file = '$precheck_file';
my \$ftp_path = '$ftp_path';
my \$md5 = '$md5';
my \$file_id = $file_id;
my \$fqc_file = '$fqc_file';
my \$total_reads = $total_reads;
my \$total_bases = $total_bases;

unless (-s \$data_path) {
    unless (-s \$precheck_file) {
        my \$io = VertRes::IO->new;
        
        if (\$ftp_path =~ qr{^(?:ftp|http)://}) {
            # download to a temp-named file
            my \$got = \$io->get_remote_file(\$ftp_path,
                                             save => \$precheck_file,
                                             \$md5 ? (md5  => \$md5) : ());
            
            unless (\$got && \$got eq \$precheck_file) {
                die "Failed to download \$ftp_path to \$precheck_file\\n";
            }
        }
        elsif (-s \$ftp_path) {
            copy(\$ftp_path, \$precheck_file) || die "Failed to move \$ftp_path to \$precheck_file\\n";
            
            if (\$md5) {
                my \$fsu = VertRes::Utils::FileSystem->new();
                my \$ok = \$fsu->verify_md5(\$precheck_file, \$md5);
                unless (\$ok) {
                    # we might have the md5 of the file in its uncompressed form
                    if (\$precheck_file =~ /\.gz\$/) {
                        system("gunzip -c \$precheck_file > \$precheck_file.uncompressed");
                        my \$uncomp_md5 = \$fsu->calculate_md5(\$precheck_file.'.uncompressed');
                        if (\$uncomp_md5 eq \$md5) {
                            \$ok = 1;
                            my \$opened = open(my \$md5fh, '>', \$data_path.'.md5');
                            unless (\$opened) {
                                unlink('\$precheck_file');
                                die "Could not write to \$data_path.md5\n";
                            }
                            print \$md5fh \$uncomp_md5, "\\n";
                            close(\$md5fh);
                            \$opened = open(\$md5fh, '<', \$data_path.'.md5');
                            unless (\$opened) {
                                unlink(\$precheck_file);
                                die "Could not open \$data_path.md5\n";
                            }
                            my \$confirm_md5 = <\$md5fh>;
                            chomp(\$confirm_md5);
                            unless (\$confirm_md5 eq \$uncomp_md5) {
                                unlink(\$precheck_file);
                                die "Tried to write the actual md5 to '\$data_path.md5' (\$uncomp_md5), but something went wrong\n";
                            }
                        }
                        unlink(\$precheck_file.'.uncompressed');
                    }
                    
                    unless (\$ok) {
                        unlink(\$precheck_file);
                        die "md5 check on '\$ftp_path' after moving it to hierarchy failed! The hierarchy copy was deleted.\n";
                    }
                }
            }
        }
        else {
            die "The given path to the fastq file to import was invalid ('$ftp_path')\n";
        }
    }
    
    if (-s \$precheck_file) {
        my \$vrfile;
        if ($self->{calculate_missing_information}) {
            require VRTrack::VRTrack;
            my \$vrtrack = VRTrack::VRTrack->new({ host => '$self->{db}->{host}',
                                                   port => $self->{db}->{port},
                                                   user => '$self->{db}->{user}',
                                                   password => '$self->{db}->{password}',
                                                   database => '$self->{db}->{database}' });
            require VRTrack::File;
            \$vrfile = VRTrack::File->new(\$vrtrack, \$file_id);
        }
        
        # calculate and store the md5 if we didn't already know it
        unless (\$md5) {
            my \$fsu = VertRes::Utils::FileSystem->new();
            \$md5 = \$fsu->calculate_md5(\$precheck_file);
            \$vrfile->md5(\$md5);
            \$vrfile->update;
        }
        
        # run fastqcheck on it and confirm it matches expectations
        unless (-s \$fqc_file) {
            # make the fqc file
            my \$wrapper = VertRes::Wrapper::fastqcheck->new();
            \$wrapper->run(\$precheck_file, \$fqc_file.'.temp');
            die "fastqcheck on \$precheck_file failed - try again?" unless \$wrapper->run_status >= 1;
            die "fastqcheck failed to make the file \$fqc_file.temp" unless -s \$fqc_file.'.temp';
            
            # check it is parsable
            my \$parser = VertRes::Parser::fastqcheck->new(file => \$fqc_file.'.temp');
            my \$num_sequences = \$parser->num_sequences();
            if (\$num_sequences) {
                move(\$fqc_file.'.temp', \$fqc_file) || die "failed to rename '\$fqc_file.temp' to '\$fqc_file'";
            }
            else {
                die "fastqcheck file '\$fqc_file.temp' was created, but doesn't seem valid\n";
            }
        }
        
        if (-s \$fqc_file) {
            # check against expectation
            my \$parser = VertRes::Parser::fastqcheck->new(file => \$fqc_file);
            my \$num_sequences = \$parser->num_sequences();
            my \$num_bases = \$parser->total_length();
            
            my \$ok = 0;
            if (\$total_reads) {
                if (\$num_sequences == \$total_reads) {
                    \$ok++;
                }
                else {
                    warn "fastqcheck reports number of reads as \$num_sequences, not the expected \$total_reads\n";
                }
            }
            else {
                \$vrfile->raw_reads(\$num_sequences);
                \$ok += \$vrfile->update || 0;
            }
            
            if (\$total_bases) {
                if (\$num_bases == \$total_bases) {
                    \$ok++;
                }
                else {
                    warn "fastqcheck reports number of bases as \$num_bases, not the expected \$total_bases\n";
                }
            }
            else {
                \$vrfile->raw_bases(\$num_sequences);
                \$ok += \$vrfile->update || 0;
            }
            
            if (\$ok == 2) {
                move(\$precheck_file, \$data_path) || die "Could not rename '\$precheck_file' to '\$data_path'\n";
            }
            else {
                die "fastqcheck of fastq file doesn't match expectations, or unable to store\n";
            }
        }
    }
}

exit;
        };
        close $scriptfh;
        
        my $job_name = $self->{prefix}.'import_'.$fastq;
        $self->archive_bsub_files($lane_path, $job_name);
        
        $script_name =~ s/ /\\ /g;
        LSF::run($action_lock, $lane_path, $job_name, $self, qq{perl -w $script_name});
    }
    
    my $male_file = $self->{fsu}->catfile($lane_path, 'MALE');
    my $female_file = $self->{fsu}->catfile($lane_path, 'FEMALE');
    my $unknown_file = $self->{fsu}->catfile($lane_path, 'UNKNOWN');
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

=head2 update_db

 Title   : update_db
 Usage   : $obj->update_db('/path/to/lane', 'lock_filename');
 Function: Update the  database to indicate which fastqs have now been
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
    my @md5_files;
    for my $file_obj (@files) {
        my $file = $file_obj->name;
        my $fastq = $file_obj->hierarchy_name;
        my $data_path = $self->{fsu}->catfile($lane_path, $fastq);
        
        $total_files++;
        next unless -s $data_path;
        $imported_files++;
        
        # check the md5; the existing md5 in the db may correspond to the md5
        # of the uncompressed fastq, but we want the compressed md5
        my $actual_md5_file = $data_path.'.md5';
        if (-e $actual_md5_file) {
            push(@md5_files, $actual_md5_file);
            open(my $md5fh, $actual_md5_file) || $self->throw("Could not open $actual_md5_file");
            my $actual_md5 = <$md5fh>;
            close($md5fh);
            chomp($actual_md5);
            $file_obj->md5($actual_md5);
        }
        
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
        foreach my $md5_file (@md5_files) {
            unlink($md5_file);
        }
        return $self->{Yes};
    }
    else {
        $self->throw("Some database updates for $lane_path failed");
    }
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
    
    my $fsu = VertRes::Utils::FileSystem->new();
    
    my $fastqs = $self->import_fastqs_provides($lane_path);
    my @files = ('log', 'job_status');
    foreach my $fastq (@{$fastqs}) {
        next if $fastq =~ /fastqcheck/;
        $fastq = basename($fastq);
        push(@files, 'import_'.$fastq.'.pl');
        
        $self->archive_bsub_files($lane_path, $prefix.'import_'.$fastq, 1);
    }
    
    foreach my $file (@files) {
        unlink($fsu->catfile($lane_path, $prefix.$file));
    }
    
    #$fsu->rmtree($fsu->catfile($lane_path, 'split_se'));
    #$fsu->rmtree($fsu->catfile($lane_path, 'split_pe'));
    
    return $self->{Yes};
}

sub is_finished {
    my ($self, $lane_path, $action) = @_;
    
    if ($action->{name} eq 'update_db') {
        $self->{vrlane}->is_processed('import') ? return $self->{Yes} : return $self->{No};
    }
    elsif ($action->{name} eq 'fix_swaps') {
        $self->{vrlane}->is_processed('swapped') ? return $self->{No} : return $self->{Yes};
    }
    elsif ($action->{name} eq 'deletion') {
        ($self->{vrlane}->is_processed('deleted') || $self->{vrlane}->is_processed('altered_fastq')) ? return $self->{No} : return $self->{Yes};
    }
    elsif ($action->{name} eq 'cleanup') {
        return $self->{No};
    }
    
    return $self->SUPER::is_finished($lane_path, $action);
}

1;
