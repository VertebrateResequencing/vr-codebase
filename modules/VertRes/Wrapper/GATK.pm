=head1 NAME

VertRes::Wrapper::RecalQual - wrapper for Broad's bam recalibrator

=head1 SYNOPSIS

use VertRes::Wrapper::RecalQual;

my $wrapper = VertRes::Wrapper::RecalQual->new(recalibrate => 1);
$wrapper->run('input_bam', 'output_bam');

# check the status
my $status = $wrapper->run_status;
if ($status == -1) {
    # try and run it again?...
}

=head1 DESCRIPTION

Runs Broad's bam recalibrator RecalQual.py in a nice way.

For mouse, assumes you have the env variable MOUSE pointing to team145's mouse
directory.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::RecalQual;

use strict;
use warnings;
use Cwd qw (abs_path);
use File::Basename;
use File::Spec;

use base qw(VertRes::Wrapper::WrapperI);
use VertRes::Wrapper::samtools;

our $default_output_dir = './broad_evaluation_files/';
our $mouse_dbsnp_link = File::Spec->catfile($ENV{MOUSE}, 'ref', 'broad_recal_data', 'dbsnp.mouse_rod_link');
our $default_mouse_dnsnp = File::Spec->catfile($ENV{MOUSE}, 'ref', 'broad_recal_data', 'dbsnp.mouse.rod.out');


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::RecalQual->new();
 Function: Create a VertRes::Wrapper::RecalQual object.
 Returns : VertRes::Wrapper::RecalQual object
 Args    : quiet    => boolean
           recalibrate => boolean
           evaluate => boolean
           species => human|mouse
           dbsbp_file => /path/to/rod.out (currently only affects mouse)

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args,
                                  exe      => 'RecalQual.py',
                                  switches => [qw(recalibrate evaluate)]);
    
    # our bsub jobs will get killed if we don't select high-mem machines
    $self->bsub_options(M => 4000000, R => "'select[mem>4000] rusage[mem=4000]'");
    
    # reset mouse dbsnp to default
    unlink($mouse_dbsnp_link);
    system("ln -s $default_mouse_dnsnp $mouse_dbsnp_link");
    
    # allow user to change species and mouse dbsnp file
    $self->_set_from_args(\@args, methods => ['species', 'dbsnp_file']);
    
    # not user-configurable atm since this location is hard-coded within
    # RecalQual.py
    $self->{_output_dir} = $default_output_dir;
    
    return $self;
}

=head2 species

 Title   : species
 Usage   : $obj->species('mouse');
 Function: Get/set the species you are recalibrating.
 Returns : species as a string (default human)
 Args    : human|mouse

=cut

sub species {
    my ($self, $species) = @_;
    
    if ($species && ($species eq 'human' || $species eq 'mouse')) {
        $self->{_species} = $species;
    }
    
    return $self->{_species} || 'human';
}

=head2 dbsnp_file

 Title   : dbsnp_file
 Usage   : $obj->dbsnp_file('/path/to/dbsnp.rod.out');
 Function: Set the dbsnp file to use. Must already be in rod.out format. See
           http://www.broadinstitute.org/gsa/wiki/index.php/The_DBSNP_rod for
           details.
           NB: currently only affects mouse!
 Returns : n/a
 Args    : path to dbsnp.rod.out

=cut

sub dbsnp_file {
    my ($self, $dbsnp_file) = @_;
    $self->throw("dbsnp_file can only currently be set when working on mouse") unless $self->species eq 'mouse';
    
    if ($dbsnp_file && -s $dbsnp_file) {
        $dbsnp_file = abs_path($dbsnp_file);
        unlink($mouse_dbsnp_link);
        system("ln -s $dbsnp_file $mouse_dbsnp_link");
    }
}

=head2 run

 Title   : run
 Usage   : $wrapper->run('input.bam', 'output.bam');
 Function: Recalibrate and/or evaluate a bam file. Will create first
           input.bam.bai using samtools index if it doesn't already exist.
 Returns : n/a (in addition to the output.bam, you will find .png files in
           subdirs of broad_evaluation_files directory if you did evaluate())
 Args    : path to input .bam file, path to output .bam file

=cut

sub _pre_run {
    my ($self, $input_bam, $output_bam) = @_;
    $input_bam && $output_bam || $self->throw("input and output bam files are required");
    -s $input_bam || $self->throw("Your input bam file '$input_bam' is empty or non-existant!");
    ($self->recalibrate || $self->evaluate) || $self->throw("recalibrate() and/or evaluate() must have been set before running");
    
    $self->exe($self->species eq 'human' ? 'RecalQual.py' : 'RecalQual_mouse.py');
    
    # fix up the bam file so we know that RecalQual.pm will like it
    # *** doesn't actually help!
    #my $healed_bam = $input_bam;
    #$healed_bam =~ s/\.bam$/.healed.bam/;
    ## check that the healed is older than the original
    #if (-s $healed_bam) {
    #    unless ($self->_is_older($healed_bam, $input_bam)) {
    #        unlink($healed_bam);
    #    }
    #}
    #unless (-s $healed_bam) {
    #    # need to write a wrapper for this...
    #    system("java -jar $ENV{G1K}/bin/picard-tools-1.01/MergeSamFiles.jar I=$input_bam O=$healed_bam") && $self->throw("Failed to run MergeSamFiles: $! | $?");
    #}
    #$input_bam = $healed_bam;
    
    my $bai_file = $input_bam.'.bai';
    
    # check that the bai is older than the bam
    if (-e $bai_file && $self->_is_older($input_bam, $bai_file)) {
        unlink($bai_file);
    }
    
    # create a bai file if necessary
    unless (-e $bai_file) {
        # run via system instead of bsub until we have dependency stuff in place...
        my $sam_wrapper = VertRes::Wrapper::samtools->new(verbose => $self->verbose,
                                                          run_method => 'system');
        $sam_wrapper->index($input_bam, $bai_file);
        
        unless ($sam_wrapper->run_status >= 1) {
            if ($sam_wrapper->run_status == -1) {
                $self->warn("Failed to create index file '$bai_file', will try one more time and then give up...");
                sleep(5);
                $sam_wrapper->index($input_bam, $bai_file);
            }
            unless ($sam_wrapper->run_status >= 1) {
                $self->throw("Failed to create index file '$bai_file', giving up!");
            }
        }
    }
    
    $self->_set_params_string(double_dash => 1);
    mkdir($self->{_output_dir});
    
    $self->{_inbam} = $input_bam;
    $self->{_outbam} = $output_bam;
    
    return ($input_bam, $output_bam);
}

sub _post_run {
    my $self = shift;
    
    # check the output
    $self->_set_run_status(1);
    my $output_bam = $self->{_outbam};
    if (-e $output_bam) {
        # on true success, RecalQual creates an index of the output bam
        my $out_index = $output_bam;
        $out_index .= '.bai';
        unless (-s $out_index) {
            $self->_set_run_status(-1);
        }
        
        if ($self->evaluate()) {
            # .png file should have been created; we know its dir location,
            # but the file name is strange
            my $dir_name = basename($self->{_inbam});
            $dir_name =~ s/\.[^.]+$//;
            $dir_name = "$self->{_output_dir}/output.$dir_name/";
            my $png_file = `ls $dir_name/recalibrated.*.empirical_v_reported_quality.png`;
            unless ($png_file) {
                $self->_set_run_status(-1);
            }
        }
    }
    else {
        $self->_set_run_status(0);
    }
    
    # move .pngs and .cvs to somewhere nice and delete the output dir?
    # ./broad_evaluation_files/output.$input_bam_minus_suffix/ contains:
    # [initial|recalibrated].*.empirical_v_reported_quality.csv
    # ".*.png
    
    return @_;
}

1;
