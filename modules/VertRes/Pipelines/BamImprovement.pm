=head1 NAME

VertRes::Pipelines::BamImprovement - pipeline for improving bam files prior to
                                     calling

=head1 SYNOPSIS

# make the config files, which specifies the details for connecting to the
# VRTrack g1k-meta database and the data roots:
echo '__VRTrack_BamImprovement__ bamimprovement.conf' > bamimprovement.pipeline
# where bamimprovement.conf contains:
root    => '/abs/path/to/root/data/dir',
module  => 'VertRes::Pipelines::BamImprovement',
prefix  => '_',

db  => {
    database => 'g1k_meta',
    host     => 'mcs4a',
    port     => 3306,
    user     => 'vreseq_rw',
    password => 'xxxxxxx',
},

data => {
    slx_mapper => 'bwa',
    '454_mapper' => 'ssaha',
    reference => '/abs/path/to/ref.fa',
    assembly_name => 'NCBI37',
    dbsnp_rod => '/abs/path/to/dbsnp.rod',
    indel_sites => '/abs/path/to/indels.vcf',
    indel_intervals => '/abs/path/to/intervals.file',
    snp_sites => '/bas/path/to/snps.vcf'
},

# indel_intervals is the file produced by GATK RealignerTargetCreator on the
# same reference, dnsbp and indel_sites files as you're using here

# run the pipeline:
run-pipeline -c bamimprovement.pipeline -s 30

# (and make sure it keeps running by adding that last to a regular cron job)

# alternatively to determining bams to improve with the database, it can work
# on a file of directory paths by altering bamimprovement.pipeline:
# make a file of absolute paths to your desired lane directories, eg:
find $G1K/META -type d -name "*RR*" | grep -v "SOLID" > lanes.fod
echo '<lanes.fod bamimprovement.conf' > bamimprovement.pipeline
# in this case, the db => {} option must be moved in the the data => {} section
# of bamimprovement.conf

=head1 DESCRIPTION

A module for "improving" bams. It uses GATK to realign reads around known indels
and then recalibrates the quality scores. Optionally, it then does samtools
calmd -r. This should improve the quality of subsequent variant calling on the
bam.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Pipelines::BamImprovement;

use strict;
use warnings;
use VertRes::Utils::Hierarchy;
use VertRes::IO;
use VertRes::Utils::FileSystem;
use VertRes::Utils::Sam;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::File;
use File::Basename;
use File::Copy;
use Time::Format;
use LSF;

use base qw(VertRes::Pipeline);

our $actions = [ { name     => 'realign',
                   action   => \&realign,
                   requires => \&realign_requires, 
                   provides => \&realign_provides },
                 { name     => 'sort',
                   action   => \&sort,
                   requires => \&sort_requires, 
                   provides => \&sort_provides },
                 { name     => 'recalibrate',
                   action   => \&recalibrate,
                   requires => \&recalibrate_requires, 
                   provides => \&recalibrate_provides },
                 { name     => 'calmd',
                   action   => \&calmd,
                   requires => \&calmd_requires, 
                   provides => \&calmd_provides },
                 { name     => 'statistics',
                   action   => \&statistics,
                   requires => \&statistics_requires, 
                   provides => \&statistics_provides },
                 { name     => 'update_db',
                   action   => \&update_db,
                   requires => \&update_db_requires, 
                   provides => \&update_db_provides },
                 { name     => 'cleanup',
                   action   => \&cleanup,
                   requires => \&cleanup_requires, 
                   provides => \&cleanup_provides } ];

our %options = (slx_mapper => 'bwa',
                '454_mapper' => 'ssaha',
                dbsnp_rod => '/lustre/scratch102/projects/g1k/ref/broad_recal_data/dbsnp_129_b37.rod',
                indel_sites => '/lustre/scratch102/projects/g1k/ref/broad_recal_data/pilot_data/1kg.pilot_release.merged.indels.sites.hg19.vcf',
                indel_intervals => '/lustre/scratch102/projects/g1k/ref/broad_recal_data/pilot_data/indel.dbsnp_129_b37-vs-pilot.intervals',
                snp_sites => '/lustre/scratch102/projects/g1k/ref/broad_recal_data/pilot_data/pilot.snps.bed',
                calmdr => 0,
                calmde => 0,
                do_cleanup => 1);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Pipelines::BamImprovement->new(lane => '/path/to/lane');
 Function: Create a new VertRes::Pipelines::BamImprovement object;
 Returns : VertRes::Pipelines::BamImprovement object
 Args    : lane => 'readgroup_id'
           lane_path => '/path/to/lane'
           do_cleanup => boolean (default true: run the cleanup action)
           slx_mapper => 'bwa'|'maq' (default bwa; the mapper you used for
                                      mapping SLX lanes)
           '454_mapper' => 'ssaha', (default ssaha; the mapper you used for
                                     mapping 454 lanes)
           
           reference => '/path/to/ref.fa' (no default, either this or the
                        male_reference and female_reference pair of args must be
                        supplied)
           assembly_name => 'NCBI37' (no default, must be set to the name of
                                      the reference)
           dbsnp_rod => '/path/to/dbsnp.rod' (defaults to dnsnp129_b37)
           indel_sites => '/path/to/indels.vcf' (known indels VCF, defaults to
                                                 1000 genomes pilot dindel calls)
           indel_intervals => '/path/to/intervals' (file produced by GATK
                              RealignerTargetCreator using the above indel_sites
                              dbsnp_rod and reference; default for 1000 genomes
                              exists)
           snp_sites => '/path/to/snps.vcf' (known snps to mask out during
                                             quality score recalibration, in
                                             addition to those in the dbsnp_rod)
           calmdr => boolean (default false; when running calmd, setting this
                              to true will use -r option)
           calmde => boolean (default false; when running calmd, setting this
                              to true will use -e option)
           release_date => 'yyyymmdd' (the release date to be included in the
                                       .bas files made; not important - defaults
                                       to today's date)

           other optional args as per VertRes::Pipeline

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(%options, actions => $actions, @args);
    
    # we should have been supplied the option 'lane_path' which tells us which
    # lane we're in, which lets us choose which mapper module to use.
    my $lane = $self->{lane} || $self->throw("lane readgroup not supplied, can't continue");
    my $lane_path = $self->{lane_path} || $self->throw("lane path not supplied, can't continue");
    
    # if we've been supplied a list of lane paths to work with, instead of
    # getting the lanes from the db, we won't have a vrlane object; make one
    if (! $self->{vrlane}) {
        $self->throw("db option was not supplied in config") unless $self->{db};
        my $vrtrack = VRTrack::VRTrack->new($self->{db}) or $self->throw("Could not connect to the database\n");
        my $vrlane  = VRTrack::Lane->new_by_name($vrtrack, $lane) or $self->throw("No such lane in the DB: [$lane]");
        $self->{vrlane} = $vrlane;
    }
    $self->{vrlane} || $self->throw("vrlane object missing");
    
    # if we've been supplied both male and female references, we'll need to
    # pick one and set the reference option
    unless ($self->{reference} || ($self->{male_reference} && $self->{female_reference})) {
        $self->throw("reference or (male_reference and female_reference) options are required");
    }
    if ($self->{male_reference}) {
        my $vrtrack = $self->{vrlane}->vrtrack;
        my $lib_id = $self->{vrlane}->library_id;
        my $lib = VRTrack::Library->new($vrtrack, $lib_id);
        my $samp_id = $lib->sample_id;
        my $samp = VRTrack::Sample->new($vrtrack, $samp_id);
        my $individual = $samp->individual();
        my $gender = $individual->sex();
        
        if ($gender eq 'F') {
            $self->{reference} = $self->{female_reference};
        }
        elsif ($gender eq 'M') {
            $self->{reference} = $self->{male_reference};
        }
        else {
            $self->throw("gender could not be determined for lane $lane");
        }
    }
    $self->throw("assembly_name must be supplied in conf") unless $self->{assembly_name};
    
    # get a list of bams in this lane we want to improve
    my $hu = VertRes::Utils::Hierarchy->new(verbose => $self->verbose);
    $self->{assembly_name} || $self->throw("no assembly_name!");
    my @bams = $hu->lane_bams($lane_path, vrtrack => $self->{vrlane}->vrtrack,
                                          assembly_name => $self->{assembly_name},
                                          slx_mapper => $self->{slx_mapper},
                                          '454_mapper' => $self->{'454_mapper'});
    @bams || $self->throw("no bams to improve in lane $lane_path!");
    $self->{in_bams} = \@bams;
    $self->{mapper_class} = $hu->{mapper_class};
    $self->{mapper_obj} = $hu->{mapper_obj};
    $self->{mapstats_obj} = $hu->{mapstats_obj};
    
    unless (defined $self->{release_date}) {
        $self->{release_date} = "$time{'yyyymmdd'}"; 
    }
    
    $self->{io} = VertRes::IO->new;
    $self->{fsu} = VertRes::Utils::FileSystem->new;
    
    return $self;
}

=head2 realign_requires

 Title   : realign_requires
 Usage   : my $required_files = $obj->realign_requires('/path/to/lane');
 Function: Find out what files the realign action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub realign_requires {
    my $self = shift;
    return [@{$self->{in_bams}},
            $self->{indel_sites},
            $self->{dbsnp_rod},
            $self->{reference},
            $self->{indel_intervals}];
}

=head2 realign_provides

 Title   : realign_provides
 Usage   : my $provided_files = $obj->realign_provides('/path/to/lane');
 Function: Find out what files the realign action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub realign_provides {
    my ($self, $lane_path) = @_;
    
    my @provides;
    
    foreach my $in_bam (@{$self->{in_bams}}) {
        my $base = basename($in_bam);
        push(@provides, '.realign_complete_'.$base);
    }
    
    return \@provides;
}

sub _bam_name_conversion {
    my ($self, $in_bam) = @_;
    my $rel_bam = $in_bam;
    $rel_bam =~ s/.recal//;
    $rel_bam =~ s/.sorted//;
    $rel_bam =~ s/\.bam$/.realigned.bam/;
    my $sorted_bam = $rel_bam;
    $sorted_bam =~ s/\.bam$/.sorted.bam/;
    my $recal_bam = $sorted_bam;
    $recal_bam =~ s/\.bam$/.recal.bam/;
    my $calmd_bam = $recal_bam;
    $calmd_bam =~ s/\.bam$/.calmd.bam/;
    return ($rel_bam, $sorted_bam, $recal_bam, $calmd_bam);
}

=head2 realign

 Title   : realign
 Usage   : $obj->realign('/path/to/lane', 'lock_filename');
 Function: realign bam files around known indels.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub realign {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = '-q normal -M4800000 -R \'select[mem>4800] rusage[mem=4800]\'';
    my $verbose = $self->verbose;
    
    foreach my $in_bam (@{$self->{in_bams}}) {
        my $base = basename($in_bam);
        my ($rel_bam) = $self->_bam_name_conversion($in_bam);
        
        next if -s $rel_bam;
        
        my $working_bam = $rel_bam;
        $working_bam =~ s/\.bam$/.working.bam/;
        my $done_file = $self->{fsu}->catfile($lane_path, '.realign_complete_'.$base);
        
        # run realign in an LSF call to a temp script
        my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."realign_$base.pl");
        
        open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
        print $scriptfh qq{
use strict;
use VertRes::Wrapper::GATK;
use VertRes::Utils::Sam;
use File::Copy;

my \$in_bam = '$in_bam';
my \$rel_bam = '$rel_bam';
my \$working_bam = '$working_bam';
my \$done_file = '$done_file';
my \$intervals_file = '$self->{indel_intervals}';

my \$gatk = VertRes::Wrapper::GATK->new(verbose => $verbose,
                                        java_memory => 3800,
                                        dbsnp => '$self->{dbsnp_rod}',
                                        reference => '$self->{reference}',
                                        build => '$self->{assembly_name}');
\$gatk->set_b('indels,VCF,$self->{indel_sites}');

# do the realignment, generating an uncompressed, name-sorted bam
unless (-s \$rel_bam) {
    \$gatk->indel_realigner(\$in_bam, \$intervals_file, \$working_bam,
                            useOnlyKnownIndels => 1,
                            LODThresholdForCleaning => 0.4,
                            bam_compression => 0);
}

# check for truncation
if (-s \$working_bam) {
    my \$su = VertRes::Utils::Sam->new;
    my \$orig_records = \$su->num_bam_records(\$in_bam);
    my \$new_records = \$su->num_bam_records(\$working_bam);
    
    if (\$orig_records == \$new_records && \$new_records > 10) {
        move(\$working_bam, \$rel_bam) || die "Could not rename \$working_bam to \$rel_bam\n";
        
        # mark that we've completed
        open(my \$dfh, '>', \$done_file) || die "Could not write to \$done_file";
        print \$dfh "done\n";
        close(\$dfh);
    }
    else {
        move(\$working_bam, "\$rel_bam.bad");
        die "\$working_bam is bad (\$new_records records vs \$orig_records), renaming to \$rel_bam.bad";
    }
}

exit;
        };
        close $scriptfh;
        
        my $job_name = $self->{prefix}.'realign_'.$base;
        $self->archive_bsub_files($lane_path, $job_name);
        
        LSF::run($action_lock, $lane_path, $job_name, $self, qq{perl -w $script_name});
    }
    
    # we've only submitted to LSF, so it won't have finished; we always return
    # that we didn't complete
    return $self->{No};
}

=head2 sort_requires

 Title   : sort_requires
 Usage   : my $required_files = $obj->sort_requires('/path/to/lane');
 Function: Find out what files the sort action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub sort_requires {
    my ($self, $lane_path) = @_;
    
    my @requires;
    
    foreach my $in_bam (@{$self->{in_bams}}) {
        my ($rel_bam) = $self->_bam_name_conversion($in_bam);
        my $done_file = $self->{fsu}->catfile($lane_path, '.sort_complete_'.basename($in_bam));
        next if -s $done_file;
        push(@requires, $rel_bam);
    }
    
    return \@requires;
}

=head2 sort_provides

 Title   : sort_provides
 Usage   : my $provided_files = $obj->sort_provides('/path/to/lane');
 Function: Find out what files the sort action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub sort_provides {
    my ($self, $lane_path) = @_;
    
    my @provides;
    foreach my $in_bam (@{$self->{in_bams}}) {
        my $base = basename($in_bam);
        push(@provides, '.sort_complete_'.$base);
    }
    
    return \@provides;
}

=head2 sort

 Title   : sort
 Usage   : $obj->sort('/path/to/lane', 'lock_filename');
 Function: sorts and fixes mates in the realigned bam.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub sort {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = '-q normal -M6800000 -R \'select[mem>6800] rusage[mem=6800]\'';
    
    foreach my $in_bam (@{$self->{in_bams}}) {
        my $base = basename($in_bam);
        my ($rel_bam, $final_bam) = $self->_bam_name_conversion($in_bam);
        
        next if -s $final_bam;
        
        my $working_bam = $final_bam;
        $working_bam =~ s/\.bam$/.working.bam/;
        my $done_file = $self->{fsu}->catfile($lane_path, '.sort_complete_'.$base);
        
        # run sort in an LSF call to a temp script
        my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."sort_$base.pl");
        
        open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
        print $scriptfh qq{
use strict;
use VertRes::Wrapper::picard;
use VertRes::Utils::Sam;
use File::Copy;

my \$in_bam = '$rel_bam';
my \$final_bam = '$final_bam';
my \$working_bam = '$working_bam';
my \$done_file = '$done_file';

# sort and fix mates
unless (-s \$final_bam) {
    my \$picard = VertRes::Wrapper::picard->new();
    \$picard->FixMateInformation(\$in_bam, \$working_bam, COMPRESSION_LEVEL => 0);
}

# check for truncation
if (-s \$working_bam) {
    my \$su = VertRes::Utils::Sam->new;
    my \$orig_records = \$su->num_bam_records(\$in_bam);
    my \$new_records = \$su->num_bam_records(\$working_bam);
    
    if (\$orig_records == \$new_records && \$new_records > 10) {
        move(\$working_bam, \$final_bam) || die "Could not rename \$working_bam to \$final_bam\n";
        
        # mark that we've completed
        open(my \$dfh, '>', \$done_file) || die "Could not write to \$done_file";
        print \$dfh "done\n";
        close(\$dfh);
    }
    else {
        move(\$working_bam, "\$final_bam.bad");
        die "\$working_bam is bad (\$new_records records vs \$orig_records), renaming to \$final_bam.bad";
    }
}

exit;
        };
        close $scriptfh;
        
        my $job_name = $self->{prefix}.'sort_'.$base;
        $self->archive_bsub_files($lane_path, $job_name);
        
        LSF::run($action_lock, $lane_path, $job_name, $self, qq{perl -w $script_name});
    }
    
    $self->{bsub_opts} = $orig_bsub_opts;
    return $self->{No};
}

=head2 recalibrate_requires

 Title   : recalibrate_requires
 Usage   : my $required_files = $obj->recalibrate_requires('/path/to/lane');
 Function: Find out what files the recalibrate action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub recalibrate_requires {
    my ($self, $lane_path) = @_;
    
    # we need bams
    my @requires;
    foreach my $in_bam (@{$self->{in_bams}}) {
        my (undef, $sort_bam, $recal_bam) = $self->_bam_name_conversion($in_bam);
        next if -s $recal_bam;
        push(@requires, $sort_bam);
    }
    
    # and reference-related files
    my $ref = $self->{reference};
    my $dict = $ref;
    $dict =~ s/\.[^\.]+$/.dict/;
    push(@requires, $ref);
    push(@requires, $dict);
    push(@requires, $self->{dbsnp_rod});
    push(@requires, $self->{snp_sites});
    
    return \@requires;
}

=head2 recalibrate_provides

 Title   : recalibrate_provides
 Usage   : my $provided_files = $obj->recalibrate_provides('/path/to/lane');
 Function: Find out what files the recalibrate action generates on success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub recalibrate_provides {
    my ($self, $lane_path) = @_;
    
    my @provides;
    foreach my $in_bam (@{$self->{in_bams}}) {
        my $base = basename($in_bam);
        push(@provides, '.recalibrate_complete_'.$base);
    }
    
    return \@provides;
}

=head2 recalibrate

 Title   : recalibrate
 Usage   : $obj->recalibrate('/path/to/lane', 'lock_filename');
 Function: Recalibrate the quality values of bams. (Original quality values are
           retained in the file.)
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub recalibrate {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = '-q long -M6800000 -R \'select[mem>6800] rusage[mem=6800]\'';
    my $verbose = $self->verbose;
    
    foreach my $bam (@{$self->{in_bams}}) {
        my $base = basename($bam);
        my (undef, $in_bam, $out_bam) = $self->_bam_name_conversion($bam);
        
        next if -s $out_bam;
        my $working_bam = $out_bam;
        $working_bam =~ s/\.bam$/.working.bam/;
        my $done_file = $self->{fsu}->catfile($lane_path, '.recalibrate_complete_'.$base);
        
        # run recal in an LSF call to a temp script
        my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."recalibrate_$base.pl");
        
        open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
        print $scriptfh qq{
use strict;
use VertRes::Wrapper::GATK;

my \$in_bam = '$in_bam';
my \$recal_bam = '$out_bam';
my \$done_file = '$done_file';
my \$intervals_file = '$self->{indel_intervals}';

my \$gatk = VertRes::Wrapper::GATK->new(verbose => $verbose,
                                        java_memory => 6000,
                                        dbsnp => '$self->{dbsnp_rod}',
                                        reference => '$self->{reference}',
                                        build => '$self->{assembly_name}');
\$gatk->set_b('mask,BED,$self->{snp_sites}');

# do the recalibration
unless (-s \$recal_bam) {
    \$gatk->recalibrate(\$in_bam, \$recal_bam);
    
    die "recalibration failed for \$in_bam\n" unless -s \$recal_bam;
    
    # mark that we've completed
    open(my \$dfh, '>', \$done_file) || die "Could not write to \$done_file";
    print \$dfh "done\n";
    close(\$dfh);
}

exit;
        };
        close $scriptfh;
        
        my $job_name = $self->{prefix}.'recalibrate_'.$base;
        $self->archive_bsub_files($lane_path, $job_name);
        
        LSF::run($action_lock, $lane_path, $job_name, $self, qq{perl -w $script_name});
    }
    
    $self->{bsub_opts} = $orig_bsub_opts;
    return $self->{No};
}

=head2 calmd_requires

 Title   : statistics_requires
 Usage   : my $required_files = $obj->calmd_requires('/path/to/lane');
 Function: Find out what files the calmd action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub calmd_requires {
    my ($self, $lane_path) = @_;
    
    my @recal_bams;
    foreach my $bam (@{$self->{in_bams}}) {
        my (undef, undef, $recal_bam, $calmd_bam) = $self->_bam_name_conversion($bam);
        next if -s $calmd_bam;
        push(@recal_bams, $recal_bam);
    }
    
    return \@recal_bams;
}

=head2 calmd_provides

 Title   : calmd_provides
 Usage   : my $provided_files = $obj->calmd_provides('/path/to/lane');
 Function: Find out what files the calmd action generates on
           success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub calmd_provides {
    my ($self, $lane_path) = @_;
    
    my @calmd_bams;
    foreach my $bam (@{$self->{in_bams}}) {
        my (undef, undef, undef, $calmd_bam) = $self->_bam_name_conversion($bam);
        push(@calmd_bams, $calmd_bam);
    }
    
    return \@calmd_bams;
}

=head2 calmd

 Title   : calmd
 Usage   : $obj->calmd('/path/to/lane', 'lock_filename');
 Function: Recalculates the MD/NM tags, and optionally does read-independent
           local realignment, altering quality values. Also optionally replaces
           reference bases with '='.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub calmd {
    my ($self, $lane_path, $action_lock) = @_;
    
    my $orig_bsub_opts = $self->{bsub_opts};
    $self->{bsub_opts} = '-q normal -M3800000 -R \'select[mem>3800] rusage[mem=3800]\'';
    
    my $r = $self->{calmdr} ? 'r => 1' : 'r => 0';
    my $e = $self->{calmde} ? 'e => 1' : 'e => 0';
    
    foreach my $in_bam (@{$self->{in_bams}}) {
        my $base = basename($in_bam);
        my (undef, undef, $recal_bam, $final_bam) = $self->_bam_name_conversion($in_bam);
        
        next if -s $final_bam;
        
        # run calmd in an LSF call to a temp script
        my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."calmd_$base.pl");
        
        open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
        print $scriptfh qq{
use strict;
use VertRes::Wrapper::samtools;

my \$in_bam = '$recal_bam';
my \$final_bam = '$final_bam';

# run calmd
unless (-s \$final_bam) {
    my \$s = VertRes::Wrapper::samtools->new(verbose => 1);
    \$s->calmd_and_check(\$in_bam, '$self->{reference}', \$final_bam, $r, $e);
    \$s->run_status() >= 1 || die "calmd failed\n";
}

exit;
        };
        close $scriptfh;
        
        my $job_name = $self->{prefix}.'calmd_'.$base;
        $self->archive_bsub_files($lane_path, $job_name);
        
        LSF::run($action_lock, $lane_path, $job_name, $self, qq{perl -w $script_name});
    }
    
    $self->{bsub_opts} = $orig_bsub_opts;
    return $self->{No};
}

=head2 statistics_requires

 Title   : statistics_requires
 Usage   : my $required_files = $obj->statistics_requires('/path/to/lane');
 Function: Find out what files the statistics action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut

sub statistics_requires {
    my ($self, $lane_path) = @_;
    
    my @requires = @{$self->calmd_provides($lane_path)};
    
    @requires || $self->throw("Something went wrong; we don't seem to require any bams!");
    
    return \@requires;
}

=head2 statistics_provides

 Title   : statistics_provides
 Usage   : my $provided_files = $obj->statistics_provides('/path/to/lane');
 Function: Find out what files the statistics action generates on
           success.
 Returns : array ref of file names
 Args    : lane path

=cut

sub statistics_provides {
    my ($self, $lane_path) = @_;
    
    my @provides;
    foreach my $bam (@{$self->statistics_requires($lane_path)}) {
        foreach my $suffix ('bas', 'flagstat') {
            push(@provides, $bam.'.'.$suffix);
        }
    }
    
    @provides || $self->throw("Something went wrong; we don't seem to provide any stat files!");
    
    return \@provides;
}

=head2 statistics

 Title   : statistics
 Usage   : $obj->statistics('/path/to/lane', 'lock_filename');
 Function: Creates statistic files for the bam files.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub statistics {
    my ($self, $lane_path, $action_lock) = @_;
    
    # setup filenames etc. we'll use within our temp script
    my $mapper_class = $self->{mapper_class};
    my $verbose = $self->verbose;
    my $release_date = $self->{release_date};
    
    # we treat read 0 (single ended - se) and read1+2 (paired ended - pe)
    # independantly.
    foreach my $bam_file (@{$self->statistics_requires($lane_path)}) {
        my $basename = basename($bam_file);
        my $script_name = $self->{fsu}->catfile($lane_path, $self->{prefix}."statistics_$basename.pl");
        
        my @stat_files;
        foreach my $suffix ('bas', 'flagstat') {
            push(@stat_files, $bam_file.'.'.$suffix);
        }
        
        # run this in an LSF call to a temp script
        open(my $scriptfh, '>', $script_name) or $self->throw("Couldn't write to temp script $script_name: $!");
        print $scriptfh qq{
use strict;
use VertRes::Utils::Sam;

my \$sam_util = VertRes::Utils::Sam->new(verbose => $verbose);

my \$num_present = 0;
foreach my \$stat_file (qw(@stat_files)) {
    \$num_present++ if -s \$stat_file;
}
unless (\$num_present == ($#stat_files + 1)) {
    my \$ok = \$sam_util->stats('$release_date', '$bam_file');
    
    unless (\$ok) {
        foreach my \$stat_file (qw(@stat_files)) {
            unlink(\$stat_file);
        }
        die "Failed to get stats for the bam '$bam_file'!\n";
    }
}

exit;
        };
        close $scriptfh;
        
        my $job_name = $self->{prefix}.'statistics_'.$basename;
        $self->archive_bsub_files($lane_path, $job_name);
        
        LSF::run($action_lock, $lane_path, $job_name, $self->{mapper_obj}->_bsub_opts($lane_path, 'statistics'), qq{perl -w $script_name});
    }
    
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
    return $self->statistics_provides($lane_path);
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
 Function: Records in the database that the lane has been improved.
 Returns : $VertRes::Pipeline::Yes or No, depending on if the action completed.
 Args    : lane path, name of lock file to use

=cut

sub update_db {
    my ($self, $lane_path, $action_lock) = @_;
    
    # get the bas files that contain our mapping stats
    my $files = $self->update_db_requires($lane_path);
    my @bas_files;
    foreach my $file (@{$files}) {
        next unless $file =~ /\.bas$/;
        push(@bas_files, $file);
        -s $bas_files[-1] || $self->throw("Expected bas file $bas_files[-1] but it didn't exist!");
    }
    
    @bas_files || $self->throw("There were no bas files!");
    
    my $vrlane = $self->{vrlane};
    my $vrtrack = $vrlane->vrtrack;
    
    return $self->{Yes} if $vrlane->is_processed('improved');
    
    $vrtrack->transaction_start();
    
    # get the mapping stats from each bas file
    my %stats;
    foreach my $file (@bas_files) {
        my $bp = VertRes::Parser::bas->new(file => $file);
        my $rh = $bp->result_holder;
        my $ok = $bp->next_result; # we'll only ever have one line, since this
                                   # is only one read group
        $ok || $self->throw("bas file '$file' had no entries, suggesting the bam file is empty, which should never happen!");
        
        # We can have up to 2 bas files, one representing single ended mapping,
        # the other paired, so we add where appropriate and only set insert size
        # for the paired
        $stats{raw_reads} += $rh->[9];
        $stats{raw_bases} += $rh->[7];
        $stats{reads_mapped} += $rh->[10];
        $stats{reads_paired} += $rh->[12];
        $stats{bases_mapped} += $rh->[8];
        $stats{mean_insert} = $rh->[15] if $rh->[15];
        $stats{sd_insert} = $rh->[16] if $rh->[15];
    }
    
    $stats{raw_bases} || $self->throw("bas file(s) for $lane_path indicate no raw bases; bam file must be empty, which should never happen!");
    
    # add mapping details to db (overwriting any existing values, but checking
    # existing values so we can work out if we need to update())
    my $mapping = $self->{mapstats_obj};
    my $needs_update = 0;
    foreach my $method (qw(raw_reads raw_bases reads_mapped reads_paired bases_mapped mean_insert sd_insert)) {
        defined $stats{$method} || next;
        my $existing = $mapping->$method;
        if (! defined $existing || $existing != $stats{$method}) {
            $needs_update = 1;
            last;
        }
    }
    if ($needs_update) {
        $mapping->raw_reads($stats{raw_reads});
        $mapping->raw_bases($stats{raw_bases});
        $mapping->reads_mapped($stats{reads_mapped});
        $mapping->reads_paired($stats{reads_paired});
        $mapping->bases_mapped($stats{bases_mapped});
        $mapping->mean_insert($stats{mean_insert}) if $stats{mean_insert};
        $mapping->sd_insert($stats{sd_insert}) if $stats{sd_insert};
        $mapping->update || $self->throw("Unable to set mapping details on lane $lane_path");
    }
    
    # set improved status
    $vrlane->is_processed('improved', 1);
    $vrlane->update() || $self->throw("Unable to set improved status on lane $lane_path");
    
    $vrtrack->transaction_commit();
    
    return $self->{Yes};
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
 Returns : array ref of file names
 Args    : lane path

=cut

sub cleanup_provides {
    return [];
}

=head2 cleanup

 Title   : cleanup
 Usage   : $obj->cleanup('/path/to/lane', 'lock_filename');
 Function: Unlink all the pipeline-related files (_*) in a lane, as well as the
           no-longer needed intermediate bam files
 Returns : $VertRes::Pipeline::Yes
 Args    : lane path, name of lock file to use

=cut

sub cleanup {
    my ($self, $lane_path, $action_lock) = @_;
    return $self->{Yes} unless $self->{do_cleanup};
    
    my $prefix = $self->{prefix};
    
    foreach my $file (qw(log job_status)) {
        unlink($self->{fsu}->catfile($lane_path, $prefix.$file));
    }
    
    foreach my $in_bam (@{$self->{in_bams}}) {
        my $base = basename($in_bam);
        
        foreach my $action ('realign', 'sort', 'recalibrate', 'statistics') {
            foreach my $suffix (qw(o e pl)) {
                my $job_name = $prefix.$action.'_'.$base;
                if ($suffix eq 'o') {
                    $self->archive_bsub_files($lane_path, $job_name, 1);
                }
                unlink($self->{fsu}->catfile($lane_path, $job_name.'.'.$suffix));
            }
        }
    }
    
    unlink($self->{fsu}->catfile($lane_path, 'GATK_Error.log'));
    
    return $self->{Yes};
}

sub is_finished {
    my ($self, $lane_path, $action) = @_;
    
    # so that we can delete the temp bams as we go along, and not redo actions
    if ($action->{name} eq 'realign') {
        foreach my $in_bam (@{$self->{in_bams}}) {
            my ($realign_bam) = $self->_bam_name_conversion($in_bam);
            
            my $base = basename($in_bam);
            my $done_file = $self->{fsu}->catfile($lane_path, '.realign_complete_'.$base);
            
            # to avoid doubling disc space requirements we want to delete the
            # original bam (replacing it with an empty file to keep the first
            # action happy), but GATK realigner had bugs and might still have
            # bugs, so we want to be able to repeat - so don't delete for now.
            #if (-s $done_file && -s $in_bam && -s $realign_bam) {
                #unlink($in_bam);
                #unlink($in_bam.'.bai');
                #system("touch $in_bam");
            #}
        }
    }
    elsif ($action->{name} eq 'sort') {
        foreach my $in_bam (@{$self->{in_bams}}) {
            my ($realign_bam, $sorted_bam) = $self->_bam_name_conversion($in_bam);
            
            my $base = basename($in_bam);
            my $done_file = $self->{fsu}->catfile($lane_path, '.sort_complete_'.$base);
            
            if (-s $done_file && -s $realign_bam && -s $sorted_bam) {
                unlink($realign_bam);
                unlink($realign_bam.'.bai');
            }
        }
    }
    elsif ($action->{name} eq 'recalibrate') {
        foreach my $in_bam (@{$self->{in_bams}}) {
            my (undef, $sorted_bam, $recal_bam) = $self->_bam_name_conversion($in_bam);
            if (-s $recal_bam && -s $sorted_bam) {
                # recal processes internally checks for truncation, so if this
                # file exists it is OK
                
                # delete the temp sorted realigned bam now that we've finished
                # recalibrating it
                unlink($sorted_bam);
                unlink($sorted_bam.'.bai');
            }
        }
    }
    elsif ($action->{name} eq 'calmd') {
        foreach my $in_bam (@{$self->{in_bams}}) {
            my (undef, undef, $recal_bam, $calmd_bam) = $self->_bam_name_conversion($in_bam);
            if (-s $calmd_bam && -s $recal_bam) {
                # calmd processe internally checks for truncation, so if this
                # file exists it is OK
                
                # delete the temp sorted realigned bam now that we've finished
                # recalibrating it
                unlink($recal_bam);
                unlink($recal_bam.'.bai');
            }
        }
    }
    elsif ($action->{name} eq 'statistics' && $self->{bam_symlink}) {
        foreach my $file (@{$self->statistics_provides($lane_path)}) {
            # Check if each file exists, if so symlink
            my $fullpath = $self->{fsu}->catfile($lane_path, $file);
            my @partnames = split(/\.bam/, $file);
            my @mappers = split(/::/, $self->{mapper_class});
            my $newname = $self->{fsu}->catfile($lane_path, $self->{assembly_name}.'_'.$mappers[-1].".bam".$partnames[-1]);
            if (-e $fullpath) {
                if (!-e $newname) {
                    symlink $fullpath, $newname;
                }
            }
        }
    }
    elsif ($action->{name} eq 'cleanup' || $action->{name} eq 'update_db') {
        return $self->{No};
    }
    
    return $self->SUPER::is_finished($lane_path, $action);
}

1;
