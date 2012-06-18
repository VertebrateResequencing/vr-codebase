package VertRes::Pipelines::Import;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use LSF;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::File;
use VertRes::Parser::fastqcheck;
use VertRes::Pipelines::Import_iRODS;
use VertRes::Pipelines::Import_Bam;
use File::Spec;

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
    'fastqcheck'      => 'fastqcheck',
    'mpsa'            => '/software/solexa/bin/mpsa_download',
    'bsub_opts'       => "-q normal -R 'select[type==X86_64] rusage[thouio=1]'",
    'local_bam_dir'     => '',
};


# --------- OO stuff --------------

=head2 new

        Example    : my $qc = VertRes::Pipelines::TrackDummy->new( files=>[2451_6_1.fastq, 2451_6_2.fastq] );
        Options    : See Pipeline.pm for general options.

                    fastqcheck      .. The fastqcheck executable
                    files           .. Array reference to the list of files to be imported.
                    mpsa            .. The mpsa executable
                    paired          .. Is the lane from the paired end sequencing.

=cut

sub new 
{
    my ($class, %args) = @_;

    if ( exists($args{files}) )
    {
        # This should be always true
        my $bams_requested = 1;
        my $bams_local = 0;
        for my $file (@{$args{files}})
        {
            my $local_file = $file;
            if ( !($file=~/\.bam$/i) ) { $bams_requested=0; last; }
            if ( $args{local_bam_dir} ) { $local_file = File::Spec->catfile($args{local_bam_dir}, $file); }
            if ( -s $local_file ) { $bams_local=1; last; }
        }
        if ( $bams_requested )
        {
            # Not sure if this is going to work and most certainly is not clean.
            if ($bams_local)
            {
                return VertRes::Pipelines::Import_Bam->new(%args);
            }
            else
            {
                return VertRes::Pipelines::Import_iRODS->new(%args);
            }
        }
    }

    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,%args);
    $self->write_logs(1);

    if ( !$$self{mpsa} ) { $self->throw("Missing the option mpsa.\n"); }
    if ( !$$self{fastqcheck} ) { $self->throw("Missing the option fastqcheck.\n"); }
    if ( !$$self{files} ) { $self->throw("Missing the option files.\n"); }
    if ( !exists($$self{paired}) ) { $self->throw("Missing the option paired.\n"); }

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

# Check database for import filenames and set @provides to 
# list of gzipped filenames.
sub get_fastqs_provides
{
    my ($self) = @_;

    my @provides = ();
    if ( !$$self{db} ) { $self->throw("Expected the db key.\n"); }

    my $vrtrack = VRTrack::VRTrack->new($$self{db}) or $self->throw("Could not connect to the database\n");
    my $vrlane  = VRTrack::Lane->new_by_name($vrtrack,$$self{lane}) or $self->throw("No such lane in the DB: [$$self{lane}]\n");
    my $vfiles=$vrlane->files;
    for my $vfile(@$vfiles){
        my $gzstrippedfilename=$vfile->name;
        $gzstrippedfilename=~s/\.gz//g;
        push @provides,$gzstrippedfilename.".gz";
    }

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
        if ( !-e qq[$lane_path/$file.gz.fastqcheck] || !-s qq[$lane_path/$file.gz.fastqcheck] ) { $must_be_run=1; last; }
    }
    if ( !$must_be_run ) { return $$self{No}; }

    # The _1 should come as the last one, as it is the only one checked in _provides
    my @rev_files = sort { $b cmp $a } @{$$self{files}};

    my $files    = 'q[' . join('],q[', @rev_files) . ']';
    my $prefix   = $$self{prefix};
    my $work_dir = $lane_path;

    # Create a script to be run on LSF.
    open(my $fh,'>', "$work_dir/${prefix}import_fastqs.pl") or $self->throw("$work_dir/${prefix}import_fastqs.pl: $!");
    print $fh
        qq[
use strict;
use warnings;
use VertRes::Pipelines::Import;

my \$opts = {
    fastqcheck   => q[$$self{fastqcheck}],
    lane         => q[$$self{lane}],
    mpsa         => q[$$self{mpsa}], 
    paired       => $$self{paired},
    files        => [ $files ],
};
my \$import = VertRes::Pipelines::Import->new(%\$opts);
\$import->get_files();

];

    close($fh);
    LSF::run($lock_file,$work_dir,"${prefix}import_fastqs",$self,qq[perl -w ${prefix}import_fastqs.pl]);

    return $$self{No};
}

sub get_files
{
    my ($self) = @_;

    # Download the files and their md5 from mpsa.
    my $prefix = $$self{prefix};
    my @gzip_files = ();    # to be gzipped .. only the new splitted files
    my @fastqcheck = ();    # to be fastqchecked .. only the new splitted files 
    for my $file (@{$$self{files}})
    {
	# multiplexed fastq files
	if ($file=~/#(\d+)/){
	    my $tag=$1;
	    my $realfile=$file;
	    $realfile=~s/#\d+//;
	    unless ( (-e "$file.gz" && -s "$file.gz") || (-e $file && -s $file) ){
	    Utils::CMD(qq[$$self{mpsa} -c -f $realfile > $file]);
	    my $extractedfile = $self->filter_fastq_on_tag($file,$tag);
	    Utils::CMD(qq[mv $extractedfile $file]);
	    Utils::CMD(qq[md5sum $file > $file.md5]);
	    }
	}
	else{
        Utils::CMD(qq[$$self{mpsa} -m -f $file > $file.md5]) unless (-e "$file.md5" && -s "$file.md5");
        Utils::CMD(qq[$$self{mpsa} -c -f $file > $file]) unless ( (-e "$file.gz" && -s "$file.gz") || (-e $file && -s $file) );
	}

        if ( -e $file )
        {
            # This should always be true if everything goes alright. But if the subroutine is called
            #   again after an error, the file may be already gzipped. Will be executed only once, 
            #   when $file is not gzipped yet.
            Utils::CMD(qq[md5sum -c $file.md5]);
        }

        if ( $$self{paired} && $file=~/^(\d+)_s_(\d+)\./ )
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

    # Run fastqcheck for all fastqs.
    for my $file (@fastqcheck)
    {
        if ( -e $file && (! -e "$file.md5" || !-s "$file.md5") ) { Utils::CMD(qq[md5sum $file > $file.md5]); }
        if ( -e "$file.gz.fastqcheck" && -s "$file.gz.fastqcheck" ) { next; }

        if ( -e $file )
        {
            Utils::CMD(qq[cat $file | $$self{fastqcheck} > $file.gz.fastqcheck.part]);
            if ( -s "$file.gz.fastqcheck.part" ) 
            { 
                rename("$file.gz.fastqcheck.part","$file.gz.fastqcheck") or $self->throw("rename $file.gz.fastqcheck.part $file.gz.fastqcheck: $!"); 
            }
        }
        elsif ( -e "$file.gz" )
        {
            Utils::CMD(qq[zcat $file.gz | $$self{fastqcheck} > $file.gz.fastqcheck.part]);
            if ( -s "$file.gz.fastqcheck.part" ) 
            { 
                rename("$file.gz.fastqcheck.part","$file.gz.fastqcheck") or $self->throw("rename $file.gz.fastqcheck.part $file.gz.fastqcheck: $!"); 
            }
        }
    }

    # Gzip the fastqs
    for my $file (@gzip_files)
    {
        if ( -e "$file.gz" && -s "$file.gz" ) { next; }

        Utils::CMD(qq[mv $file $file.x; rm -f $file.x.gz; gzip $file.x;]);
        Utils::CMD(qq[mv $file.x.gz $file.gz]);
        Utils::CMD(qq[md5sum $file.gz > $file.gz.md5]);
    }

    # If there are any single files (e.g. *_s_*) which are not paired, create a symlink to it:
    #   First find out how many fastq files with the conventional naming there are (there 
    #   shouldn't be any).
    my $i = 1;
    while ( -e "$$self{lane}_$i.fastq.gz" ) { $i++; }
    for my $file (@{$$self{files}})
    {
        if ( $$self{paired} ) { next; }
        if ( $file=~/^$$self{lane}_(\d+).fastq/ ) { next; } # This one is named correctly

        my ($run,$lane);
        if ( $file=~/^(\d+)_(\d+)/ || $file=~/^(\d+)_s_(\d+)/ ) { $run=$1; $lane=$2; }
        else { $self->throw("Weird naming, could not parse '$file'"); }

        my $name = "${1}_${2}_$i";
        if ( ! -e "$name.fastq.gz" ) 
        { 
            Utils::relative_symlink("$file.gz","$name.fastq.gz"); 
        }
        if ( -e "$file.gz.fastqcheck" && ! -e "$name.fastq.gz.fastqcheck" )
        {
            Utils::relative_symlink("$file.gz.fastqcheck","$name.fastq.gz.fastqcheck");
        }
        if ( -e "$file.md5" && ! -e "$name.fastq.md5" )
        {
            Utils::relative_symlink("$file.md5","$name.fastq.md5");
        }

        $i++;
    }

}

# Splits the fastq file into two. Assuming that sequences and qualities have even length
#   and that the first half is the forward and the second the reverse read.
sub split_single_fastq
{
    my ($self,$file,$run,$lane) = @_;
    my $fname1 = qq[${run}_${lane}_1.fastq];
    my $fname2 = qq[${run}_${lane}_2.fastq];

    if ( -e "$fname1.gz" && -e "$fname2.gz" ) { return ($fname1,$fname2); }

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

#  Filters the multiplex fastqs and extracts the sequences relevant to a particular tag	 
sub filter_fastq_on_tag
{
    my ($self,$file,$tag) = @_;
    my $fname = qq[$file.temp];

    if ( -e "$fname.gz" ) { return $fname; }

    open(my $fh_in,'<',$file) or $self->throw("$file: $!");
    open(my $fh_out,'>',$fname) or $self->throw("$fname: $!");
    while (1)
    {
	#@IL5_4821:7:1:1154:8255#0/1
	#AAAGATTTGCCTGAGTCAGGNTGCGCAAGNNNNANNNNNNNNGNGNATNAANNTTCANAAAATAAAAANNNNNNNN	 
	#+
	#BB@@@@B@@A@AA0@@A=B=&9A=?16@=%$$$5$#########($6)%48$$*/7'%&$'(6$'$&/#$$$$$#$
	my $id   = <$fh_in> || last;
	# extract only if the tag matches 
	if($id=~/#(\d+)/ && $1==$tag){
	my $seq  = <$fh_in> || last; 
	my $sep  = <$fh_in> || last; 
	my $qual = <$fh_in> || last; 
	print $fh_out $id;
	print $fh_out $seq;
	print $fh_out $sep;
	print $fh_out $qual;
	}
	else{
	# skip the next 3 lines	 
	foreach (1..3){<$fh_in>;}	 
	}
    }
    close($fh_out);	 
    close($fh_in);	 
    return $fname;	 
}


sub get_hierarchy_path
{
    my ($self) = @_;
    return "$$self{project}/$$self{sample}/$$self{technology}/$$self{library}/$$self{lane}";
}


#---------- update_db ---------------------

# Requires the gzipped fastq files.
sub update_db_requires
{
    my ($self,$lane_path) = @_;

    if ( !$$self{db} ) { $self->throw("Expected the db key.\n"); }
    my @requires = ();
    my $vrtrack = VRTrack::VRTrack->new($$self{db}) or $self->throw("Could not connect to the database\n");
    my $vrlane  = VRTrack::Lane->new_by_name($vrtrack,$$self{lane}) or $self->throw("No such lane in the DB: [$$self{lane}]\n");
    my $vfiles=$vrlane->files;
    for my $vfile(@$vfiles){
        my $gzstrippedfilename=$vfile->name;
        $gzstrippedfilename=~s/\.gz//g;        
        push @requires,$gzstrippedfilename.".gz";
    }

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

    $vrtrack->transaction_start();

    # To determine file types, a simple heuristic is used: When there are two files 
    #   lane_1.fastq.gz and lane_2.fastq.gz, fwd (1) and rev (2) will be set for file type.
    #   When only one file is present, single-end (0) will be set.
    my $i=0;
    my %processed_files;
    while (1)
    {
        $i++;
        my $name = "$$self{lane}_$i.fastq";

        # Check what fastq files actually exist in the hierarchy
        if ( ! -e "$lane_path/$name.gz" ) { last; }

        $processed_files{$name} = $i;
    }

    # do the same for every _nonhuman.fastq files
    $i=0;
    while (1)
    {
	$i++;
	my $name = "$$self{lane}_$i"."_nonhuman.fastq";

	# Check what fastq files actually exist in the hierarchy
	if ( ! -e "$lane_path/$name.gz" ) { last; }

	$processed_files{$name} = $i;
    }

    my $nfiles = scalar keys %processed_files;
    while (my ($name,$type) = each %processed_files)
    {
        # The file may be absent from the database, if it was created by splitting the _s_ fastq.
        my $vrfile = $vrlane->get_file_by_name($name);
        if ( !$vrfile ) 
        { 
            $vrfile = $vrlane->add_file($name); 
            $vrfile->hierarchy_name($name);
        }
        $vrfile->md5(`awk '{printf "%s",\$1}' $lane_path/$name.md5`);

        # Hm, this must be evaled, otherwise it dies without rollback
        my ($avg_len,$tot_len,$num_seq,$avg_qual);
        eval {
            my $fastq = VertRes::Parser::fastqcheck->new(file => "$lane_path/$name.gz.fastqcheck");
            $avg_len  = $fastq->avg_length();
            $tot_len  = $fastq->total_length();
            $num_seq  = $fastq->num_sequences();
            $avg_qual = $fastq->avg_qual();
        };
        if ( $@ )
        {
            $vrtrack->transaction_rollback();
            $self->throw("Problem reading the fastqcheck file: $lane_path/$name.gz.fastqcheck\n");
        }
        $vrfile->read_len($avg_len);
        $vrfile->raw_bases($tot_len);
        $vrfile->raw_reads($num_seq);
        $vrfile->mean_q($avg_qual); 
        $vrfile->is_processed('import',1);

        # Only the gzipped variant will carry the latest flag
        $vrfile->is_latest(0);  

        # If there is only one file, it is single-end file
        if ( $nfiles == 1 ) { $type = 0; }
        $vrfile->type($type);
        $vrfile->update();

        # Now add also the fastq.gz file into the File table. (Requested for the Mapping pipeline.)
        $self->fix_lane_file_if_exists("$name.gz",$vrlane, $vrtrack);
        my $vrfile_gz = $vrlane->get_file_by_name("$name.gz");
        if ( !$vrfile_gz )
        {
            $vrfile_gz = $vrlane->add_file("$name.gz");
            vrtrack_copy_fields($vrfile,$vrfile_gz,[qw(file_id name hierarchy_name)]);
            $vrfile_gz->hierarchy_name("$name.gz");
            $vrfile_gz->md5(`awk '{printf "%s",\$1}' $lane_path/$name.gz.md5`);
            $vrfile_gz->is_processed('import',1);
            $vrfile_gz->update();
        }
    }

    # Update also the status of the _s_ file, if any. Loop over the
    #   files passed to the pipeline and filter out those which were not
    #   processed above.
    for my $file (@{$$self{files}})
    {
        if ( exists($processed_files{$file}) ) { next; }

        my $vrfile = $vrlane->get_file_by_name($file);
        if ( !$vrfile ) { $self->throw("FIXME: no such file [$file] for the lane [$lane_path]."); }
        $vrfile->is_processed('import',1);
        $vrfile->is_latest(0);
        $vrfile->update();
    }


    # Update the lane stats
    my $read_len=0;
    my $raw_reads=0;
    my $raw_bases=0;
    my $vfiles=$vrlane->files;
    for my $vfile(@$vfiles){
	$raw_reads+=$vfile->raw_reads;
	$raw_bases+=$vfile->raw_bases;
	$read_len=$vfile->read_len;
    }
    $vrlane->raw_reads($raw_reads);
    $vrlane->raw_bases($raw_bases);
    $vrlane->read_len($read_len);


    # Finally, change the import status of the lane, so that it will not be picked up again
    #   by the run-pipeline script.
    $vrlane->is_processed('import',1);
    $vrlane->update();
    $vrtrack->transaction_commit();

    return $$self{Yes};
}

=head2 fix_lane_file_if_exists

        Example : fix_lane_file_if_exists("my_file_name", $vlane);
        Args    : filename to check
                : lane object which will in future have the file object

  This will see if a file row already exists for a lane. If its attached to a different lane, then delete it so that it doesnt cause
  issues later.
=cut

sub fix_lane_file_if_exists
{
  my ($self,$file_name, $vrlane, $vrtrack) = @_;
  
  if(VRTrack::File->is_name_in_database($vrtrack,$file_name,$file_name) == 1)
  {
    my $direct_file_object = VRTrack::File->new_by_name($vrtrack, $file_name);
    my $lane_file_object = $vrlane->get_file_by_name($file_name);
    
    if(! defined($lane_file_object)  || ( defined($lane_file_object) &&  $direct_file_object->row_id() != $lane_file_object->row_id()  ) )
    {
      # file exists already but its not attached to our lane, delete it
      $direct_file_object->delete();
    }
  }
}


=head2 vrtrack_copy_fields

        Example : vrtrack_copy_fields($src_obj,$dst_obj,[qw(lane_id library_id submission_id)]);
        Args    : Source object
                : Target object
                : Array ref to the fields which should be omitted

=cut

sub vrtrack_copy_fields
{
    my ($src,$dst,$except) = @_;

    my %omit = $except ? map { $_=>1 } @$except : ();
    my $src_fields = $src->fields_dispatch();
    my $dst_fields = $dst->fields_dispatch();

    while (my ($key,$handler) = each %$src_fields)
    {
        if ( exists($omit{$key}) ) { next; }

        my $value = &{$handler}();
        &{$$dst_fields{$key}}($value);
    }
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

