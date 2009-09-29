package VertRes::Import;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use LSF;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::File;
use VertRes::Parser::fastqcheck;


our @actions =
(
    # Creates the hierarchy path, downloads, gzips and checkfastq the fastq files.
    {
        'name'     => 'get_fastqs',
        'action'   => \&get_fastqs,
        'requires' => \&get_fastqs_requires, 
        'provides' => \&get_fastqs_provides,
    },

    # If all files were downloaded OK, update the VRTrack database.
    {
        'name'     => 'update_db',
        'action'   => \&update_db,
        'requires' => \&update_db_requires, 
        'provides' => \&update_db_provides,
    },
);

our $options = 
{
    # Executables
    'fastqcheck'      => '/nfs/sf8/G1K/bin/fastqcheck',
    'mpsa'            => '/software/solexa/bin/mpsa_download',

    'bsub_opts'       => "-q normal -R 'select[type==X86_64]'",
};


# --------- OO stuff --------------

=head2 new

        Example    : my $qc = VertRes::TrackDummy->new( files=>[2451_6_1.fastq, 2451_6_2.fastq] );
        Options    : See Pipeline.pm for general options.

                    fastqcheck      .. The fastqcheck executable
                    files           .. Array reference to the list of files to be imported.
                    mapping_root    .. The root of the mapping hierarchy, where symlinks should be created
                    mpsa            .. The mpsa executable

=cut

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,@args);
    $self->write_logs(1);

    if ( !$$self{mpsa} ) { $self->throw("Missing the option mpsa.\n"); }
    if ( !$$self{fastqcheck} ) { $self->throw("Missing the option fastqcheck.\n"); }
    if ( !$$self{mapping_root} ) { $self->throw("Missing the option mapping_root.\n"); }
    if ( !$$self{files} ) { $self->throw("Missing the option files.\n"); }

    return $self;
}


#---------- get_fastqs ---------------------

# Requires nothing
sub get_fastqs_requires
{
    my ($self) = @_;
    my @requires = ();
    return \@requires;
}

# It may provide also _2.fastq.gz, but we assume that if
#   there is _1.fastq.gz but _2 is missing, it is OK.
#
sub get_fastqs_provides
{
    my ($self) = @_;
    my @provides = ("$$self{lane}_1.fastq.gz");
    return \@provides;
}

sub get_fastqs
{
    my ($self,$lane_path,$lock_file) = @_;

    my $must_be_run = 0;
    for my $file (@{$$self{files}})
    {
        # If all gzipped files are in place, everything has been done already.
        if ( !-e qq[$lane_path/$file.gz] ) { $must_be_run=1; last; }
        if ( !-e qq[$lane_path/$file.gz.fastqcheck] ) { $must_be_run=1; last; }
    }
    if ( !$must_be_run ) { return $$self{No}; }

    my $files    = 'q[' . join('],q[', @{$$self{files}}) . ']';
    my $prefix   = $$self{prefix};
    my $work_dir = $lane_path;

    # Create a script to be run on LSF.
    open(my $fh,'>', "$work_dir/${prefix}import_fastqs.pl") or $self->throw("$work_dir/${prefix}import_fastqs.pl: $!");
    print $fh
        qq[
use strict;
use warnings;
use VertRes::Import;

my \$opts = {
    fastqcheck   => q[$$self{fastqcheck}],
    mpsa         => q[$$self{mpsa}], 
    mapping_root => q[$$self{mapping_root}],
    files        => [ $files ],
};
my \$import = VertRes::Import->new(%\$opts);
\$import->get_files();

];

    close($fh);
    LSF::run($lock_file,$work_dir,"${prefix}import_fastqs",$self,qq[perl -w ${prefix}import_fastqs.pl]);

    return $$self{No};
}

sub get_files
{
    my ($self) = @_;

    my $prefix = $$self{prefix};
    my @gzip_files = ();    # to be gzipped .. only the new splitted files
    my @fastqcheck = ();    # to be fastqchecked .. only the new splitted files 
    for my $file (@{$$self{files}})
    {
        Utils::CMD(qq[$$self{mpsa} -m -f $file > $file.md5]) unless -e "$file.md5";
        Utils::CMD(qq[$$self{mpsa} -c -f $file > $file]) unless ( -e "$file.gz" || -e $file );
        if ( -e $file )
        {
            # This should always be true if everything goes alright. But if the subroutine is called
            #   again after an error, the file may be already gzipped. Will be executed only once, 
            #   when $file is not gzipped yet.
            Utils::CMD(qq[md5sum -c $file.md5]);
        }

        if ( $file=~/^(\d+)_s_(\d+)\./ )
        {
            my $run  = $1;
            my $lane = $2;
            my ($file1,$file2) = $self->split_single_fastq($file,$run,$lane);
            push @gzip_files, $file1,$file2;
            push @fastqcheck, $file1,$file2;
        }
        push @gzip_files,$file;
        push @fastqcheck,$file;
    }

    for my $file (@fastqcheck)
    {
        if ( -e $file && ! -e "$file.md5" ) { Utils::CMD(qq[md5sum $file > $file.md5]); }
        if ( -e "$file.gz.fastqcheck" ) { next; }

        if ( -e $file )
        {
            Utils::CMD(qq[cat $file | $$self{fastqcheck} > $file.gz.fastqcheck;]);
        }
        elsif ( -e "$file.gz" )
        {
            Utils::CMD(qq[zcat $file.gz | $$self{fastqcheck} > $file.gz.fastqcheck;]);
        }
    }

    for my $file (@gzip_files)
    {
        if ( -e "$file.gz" ) { next; }

        Utils::CMD(qq[mv $file $file.x; gzip $file.x;]);
        Utils::CMD(qq[mv $file.x.gz $file.gz]);
    }
}

# Splits the fastq file into two. Assuming that sequences and qualities have even length
#   and that the first half is the forward and the second the reverse read.
sub split_single_fastq
{
    my ($self,$file,$run,$lane) = @_;
    my $fname1 = qq[${run}_${lane}_1.fastq];
    my $fname2 = qq[${run}_${lane}_2.fastq];

    if ( -e "$fname1.gz" && -e "$fname2.gz" ) { return; }

    open(my $fh_in,'<',$file) or $self->throw("$file: $!");
    open(my $fh_out1,'>',$fname1) or $self->throw("$fname1: $!");
    open(my $fh_out2,'>',$fname2) or $self->throw("$fname2: $!");

    my $len;
    while (1)
    {
        my $id   = <$fh_in> || last;
        my $seq  = <$fh_in> || last;
        my $sep  = <$fh_in> || last;
        my $qual = <$fh_in> || last;

        # @IL7_692:2:1:878:794/2
        # TTTATTATTGCCATACTATGGGCAAAGGTACACTAA
        # +
        # @@@A?AAA@CABBBBBBBBBAAA@@@?:2=;>:;<:

        chomp($seq);
        chomp($qual);

        if ( !$len ) 
        { 
            $len = length $seq;
            if ( !$len ) { $self->throw("Zero length sequence? [$len] [$file] [$seq]\n"); }
            if ( $len%2 != 0 ) { $self->throw("The length of the sequence not even: [$len] [$file] [$seq]\n"); }
            $len = $len/2;
        }

        print $fh_out1 $id;
        print $fh_out2 $id;

        print $fh_out1 substr($seq,0,$len) . "\n";
        print $fh_out2 substr($seq,$len) . "\n";

        print $fh_out1 $sep;
        print $fh_out2 $sep;

        print $fh_out1 substr($qual,0,$len) . "\n";
        print $fh_out2 substr($qual,$len) . "\n";
    }
    close($fh_out1);
    close($fh_out2);
    close($fh_in);

    return ($fname1,$fname2);
}


sub get_hierarchy_path
{
    my ($self) = @_;
    return "$$self{project}/$$self{sample}/$$self{technology}/$$self{library}/$$self{lane}";
}


#---------- update_db ---------------------

# Requires the gzipped fastq files. How many? Find out how many .md5 files there are.
sub update_db_requires
{
    my ($self,$lane_path) = @_;
    my @requires = ();
    my $i = 1;
    while ( -e "$lane_path/$$self{lane}_$i.fastq.md5" )
    {
        push @requires, "$$self{lane}_$i.fastq.gz";
        $i++;
    }
    if ( !@requires ) { @requires = ("$$self{lane}_1.fastq.gz"); }
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

sub update_db
{
    my ($self,$lane_path,$lock_file) = @_;

    if ( !$$self{db} ) { $self->throw("Expected the db key.\n"); }

    my $vrtrack = VRTrack::VRTrack->new($$self{db}) or $self->throw("Could not connect to the database\n");
    my $vrlane  = VRTrack::Lane->new_by_name($vrtrack,$$self{lane}) or $self->throw("No such lane in the DB: [$$self{lane}]\n");

    my $mapping_dir = "$$self{mapping_root}/" . $self->get_hierarchy_path();
    Utils::CMD(qq[mkdir -p '$mapping_dir']);

    $vrtrack->transaction_start();

    my $i = 0;
    while (1)
    {
        # Check what fastq files actually exist in the hierarchy and update the processed flag
        $i++;
        my $name = "$$self{lane}_$i.fastq";

        if ( ! -e "$lane_path/$name.gz" ) { last; }

        if ( ! -e "$mapping_dir/$name.gz" )
        {
            Utils::relative_symlink("$lane_path/$name.gz","$mapping_dir/$name.gz");
        }
        if ( ! -e "$mapping_dir/$name.gz.fastqcheck" )
        {
            Utils::relative_symlink("$lane_path/$name.gz.fastqcheck","$mapping_dir/$name.gz.fastqcheck");
        }

        my $vrfile = $vrlane->get_file_by_name($name);
        if ( !$vrfile ) 
        { 
            $vrfile = $vrlane->add_file($name); 
            $vrfile->hierarchy_name($name);
        }
        $vrfile->md5(`awk '{printf "%s",\$1}' $lane_path/$name.md5`);
        my $fastq = VertRes::Parser::fastqcheck->new(file => "$lane_path/$name.gz.fastqcheck");
        $vrfile->read_len($fastq->avg_length());
        $vrfile->raw_bases($fastq->total_length());
        $vrfile->raw_reads($fastq->num_sequences());
        $vrfile->mean_q($fastq->avg_qual()); 
        $vrfile->is_processed('import',1);
        $vrfile->update();
    }

    # Change also the qc status of the single _s_ file, if there is any. To find out,
    #   loop through the files supplied by run-pipeline.
    for my $file (@{$$self{files}})
    {
        if ( !($file=~/^\d+_s_\d+\./) ) { next; }

        my $vrfile = $vrlane->get_file_by_name($file);
        $vrfile->is_processed('import',1);
        $vrfile->update();
    }

    # Change the qc status of the lane.
    $vrlane->is_processed('import',1);
    $vrlane->update();
    $vrtrack->transaction_commit();

    return $$self{Yes};
}



#---------- Debugging and error reporting -----------------

sub warn
{
    my ($self,@msg) = @_;
    my $msg = join('',@msg);
    if ($self->verbose > 0) 
    {
        print STDERR $msg;
    }
    $self->log($msg);
}

sub debug
{
    my ($self,@msg) = @_;
    if ($self->verbose > 1) 
    {
        my $msg = join('',@msg);
        print STDERR $msg;
        $self->log($msg);
    }
}

sub throw
{
    my ($self,@msg) = @_;
    Utils::error(@msg);
}

sub log
{
    my ($self,@msg) = @_;

    my $msg_str = join('',@msg);
    my $status  = open(my $fh,'>>',$self->log_file);
    if ( !$status ) 
    {
        print STDERR $msg_str;
    }
    else 
    { 
        print $fh $msg_str; 
    }
    if ( $fh ) { close($fh); }
}


1;

