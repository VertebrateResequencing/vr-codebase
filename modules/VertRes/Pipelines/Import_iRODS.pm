package VertRes::Pipelines::Import_iRODS;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use VertRes::LSF;
use VRTrack::VRTrack;
use VRTrack::Lane;
use VRTrack::File;
use VertRes::Wrapper::iRODS;
use VertRes::Parser::bamcheck;
use File::Basename;

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
        'action'   => \&update_db,
        'requires' => \&update_db_requires, 
        'provides' => \&update_db_provides,
    },
);

our $options = 
{
    'bamcheck'        => 'bamcheck -q 20',
    'bsub_opts'       => "-q normal -R 'select[type==X86_64] rusage[thouio=1]'",
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
    my ($self) = @_;
    return $$self{files};
}

sub get_bams
{
    my ($self,$lane_path,$lock_file) = @_;

    my $opts = $self->dump_opts(qw(files bamcheck));

    my $prefix   = $$self{prefix};
    my $work_dir = $lane_path;

    # Create a script to be run on LSF.
    open(my $fh,'>', "$work_dir/${prefix}import_bams.pl") or $self->throw("$work_dir/${prefix}import_bams.pl: $!");
    print $fh qq[
use strict;
use warnings;
use VertRes::Pipelines::Import_iRODS;

my $opts

my \$import = VertRes::Pipelines::Import_iRODS->new(%\$opts);
\$import->get_files();

];

    close($fh);
    VertRes::LSF::run($lock_file,$work_dir,"${prefix}import_bams",$self,qq[perl -w ${prefix}import_bams.pl]);

    return $$self{No};
}

sub get_files
{
    my ($self) = @_;

    my $irods = VertRes::Wrapper::iRODS->new();

    # Get files and run bamcheck on them
    for my $file (@{$$self{files}})
    {
        my $ifile = $irods->find_file_by_name($file);
        if ( !defined $ifile ) { $self->warn("No such file in iRODS? [$file]\n"); next; }

        if ( !($ifile=~m{([^/]+)$}) ) { $self->throw("FIXME: [$ifile]"); }
        my $outfile = $1;
        if ( -e $outfile ) { next; }

        # Check if the BAM contains any reads.
        if( $irods->get_total_reads($ifile) == 0)
        {
          $self->log("No reads in BAM file");
          next;
        }

        # Get the BAM index
        my $bai = $ifile;
        $bai =~ s/\.bam$/.bai/i;
        $irods->get_file($bai,"$outfile.bai");
        chmod 0664,"$outfile.bai";

        # Get the BAM file from iRODS. From some reasons, some of the files have wrong permissions
        $irods->get_file($ifile,"$outfile.tmp");
        chmod 0664,"$outfile.tmp";

        # Get the md5sum and check
        my $md5 = $irods->get_file_md5($ifile);
        open(my $fh,'>',"$outfile.md5") or $self->throw("$outfile.md5: $!");
        print $fh "$md5  $outfile.tmp\n";
        close($fh);

        # Check that everything went alright. Although iRODS is supposed to check, local IO failures went unnoticed the other day...
        Utils::CMD(qq[md5sum --status -c $outfile.md5]);

        # Recreate the checksum file to contain the correct file name
        open($fh,'>',"$outfile.md5") or $self->throw("$outfile.md5: $!");
        print $fh "$md5  $outfile\n";
        close($fh);

        Utils::CMD(qq[$$self{bamcheck} $outfile.tmp > $outfile.tmp.bc]);
        rename("$outfile.tmp.bc","$outfile.bc") or $self->throw("rename $outfile.tmp.bc $outfile.bc: $!");
        rename("$outfile.tmp",$outfile) or $self->throw("rename $outfile.tmp $outfile: $!");
    }
}


#---------- update_db ---------------------

# Requires the gzipped fastq files. How many? Find out how many .md5 files there are.
sub update_db_requires
{
    my ($self,$lane_path) = @_;
    return $$self{files};
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

    my ($tot_tot_len,$tot_num_seq);
    for my $file (@{$$self{files}})
    {
        my $base = File::Basename::basename($file);
        # Hm, this must be evaled, otherwise it dies without rollback
        my ($avg_len,$tot_len,$num_seq,$avg_qual,$nfirst,$nlast,$is_mapped,$is_sorted,$is_paired,$ok);
        eval {
            my $pars = VertRes::Parser::bamcheck->new(file => "$lane_path/$base.bc");
            $num_seq   = $pars->get('sequences');
            $tot_len   = $pars->get('total_length');
            $avg_len   = $pars->get('avg_length');
            $avg_qual  = $pars->get('avg_qual');
            $nfirst    = $pars->get('1st_fragments');
            $nlast     = $pars->get('last_fragments');
            $is_sorted = $pars->get('is_sorted');
            $is_paired = $pars->get('is_paired');

            # One sequence without BAM_FUNMAP flag (0x0004) is enough to say that the BAM file is mapped.
            my $unmapped = $pars->get('reads_unmapped');
            $is_mapped = ( $nfirst+$nlast-$unmapped>0 ) ? 1 : 0;

            $ok = 1;
        };
        if ( !$ok )
        {
            my $err_msg = $@;
            $vrtrack->transaction_rollback();
            $self->throw("Problem reading the bamcheck file: $lane_path/$base.bc\n$err_msg\n");
        }

        my ($md5) = Utils::CMD(qq[awk '{printf "%s",\$1}' $lane_path/$base.md5]);

        my $vrfile = $vrlane->get_file_by_name($file);
        if ( !$vrfile ) { $self->throw("FIXME: the file not in the DB? [$file]"); }

        $vrfile->read_len($avg_len);
        $vrfile->raw_bases($tot_len);
        $vrfile->raw_reads($num_seq);
        $vrfile->mean_q($avg_qual);
        $vrfile->md5($md5);
        $vrfile->type( (!$nfirst or !$nlast)?0:4 );
        $vrfile->is_processed('import',1);
        if ( $is_mapped ) 
        { 
            $vrfile->is_processed('mapped',1);
            $vrlane->is_processed('mapped',1);
        }
        $vrfile->update();

        # Create mapstats, fill in mapper and assembly and rename the files by mapstat ID
        my $assembly = 'Unknown';
        my $mapper   = 'Unknown';
        my $mapper_version = 'Unknown';

        my $parser = VertRes::Parser::bam->new(file=>"$lane_path/$base");
        my @info = $parser->sequence_info();
        if ( @info ) { $assembly = $info[1]{AS}; }

        my $prg = $parser->program();
        if ( $prg )
        {
            $mapper = $prg;
            my %prgs = $parser->program_info();
            $mapper_version = $prgs{$mapper}{VN};
        }

        my $vrmapping  = $vrlane->add_mapping();
        my $vrassembly = $vrmapping->assembly($assembly);
        if ( !$vrassembly ) { $vrassembly = $vrmapping->add_assembly($assembly); }

        my $vrmapper = $vrmapping->mapper($mapper,$mapper_version);
        if ( !$vrmapper ) { $vrmapper = $vrmapping->add_mapper($mapper,$mapper_version); }
        $vrmapping->update();

        # Now rename everything to '1234.pe.raw.sorted.*'
        my $new_fn = sprintf "%d.%s.raw%s.bam", 
            $vrmapping->id(), 
            $is_paired ? 'pe' : 'se',
            $is_sorted ? '.sorted' : '';

        rename("$lane_path/$base","$lane_path/$new_fn");
        rename("$lane_path/$base.bc","$lane_path/$new_fn.bc");
        rename("$lane_path/$base.bai","$lane_path/$new_fn.bai");
        rename("$lane_path/$base.md5","$lane_path/$new_fn.md5");
        Utils::CMD(qq[echo "$md5  $new_fn" > $lane_path/$new_fn.md5]);
        
        $vrfile->hierarchy_name($new_fn);
        $vrfile->update();

        $tot_num_seq += $num_seq;
        $tot_tot_len += $tot_len;
    }

    # Finally, change the import status of the lane, so that it will not be picked up again
    #   by the run-pipeline script.
    $vrlane->is_processed('import',1);
    $vrlane->is_withdrawn(0);
    $vrlane->raw_bases($tot_tot_len);
    $vrlane->raw_reads($tot_num_seq);
    $vrlane->update();
    $vrtrack->transaction_commit();

    return $$self{Yes};
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

