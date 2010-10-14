package VertRes::Pipelines::TrackGType;
use base qw(VertRes::Pipeline);

use strict;
use warnings;
use LSF;
use HierarchyUtilities;
use VertRes::Utils::GTypeCheck;

our @actions =
(
    {
        'name'     => 'check_genotype',
        'action'   => \&check_genotype,
        'requires' => \&check_genotype_requires, 
        'provides' => \&check_genotype_provides,
    },
);


our $options = 
{
    # Executables
    'glf'            => 'glf',
    'samtools'       => 'samtools',

    'bsub_opts'      => "-q normal -M5000000 -R 'select[mem>5000] rusage[mem=5000,thouio=1]'",
    'min_glf_ratio'  => 5.0,
};


# --------- OO stuff --------------

=head2 new

        Description: This is an example of a simple Pipeline, runs a genotype check on a supplied .bam file.
        Example    : my $cg = VertRes::CheckGenotype->new('bam' => 'file.bam'); $cg->run_lane();
        Options    : See Pipeline.pm for general options.

                    # Executables
                    glf             .. glf executable
                    samtools        .. samtools executable

                    # Options specific to CheckGenotype
                    bam             .. 
                    bsub_opts       .. LSF bsub options for jobs
                    fa_ref          .. see HierarchyUtilities::lane_info
                    min_glf_ratio   .. the minimum distinctive glf likelihood ratio
                    snps            .. see HierarchyUtilities::lane_info

=cut

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$options,'actions'=>\@actions,@args);
    return $self;
}



#----------- check_genotype ---------------------

sub check_genotype_requires
{
    my ($self,$lane) = @_;

    my @requires;
    if ( exists($$self{'bam'}) ) 
    {
        my ($dir,$basename,$suffix) = Utils::basename($$self{'bam'});
        @requires = ("$basename$suffix");
    }
    else
    {
        my $lane_info  = HierarchyUtilities::lane_info($lane);
        @requires = ("$$lane_info{lane}.bam");
    }
    return \@requires;
}

sub check_genotype_provides
{
    my ($self,$lane) = @_;

    my @provides;
    if ( exists($$self{'bam'}) ) 
    {
        my ($dir,$basename,$suffix) = Utils::basename($$self{'bam'});
        @provides = ("$basename.gtype");
    }
    else
    {
        my $lane_info  = HierarchyUtilities::lane_info($lane);
        @provides = ("$$lane_info{lane}.gtype");
    }
    return \@provides;
}


# First tries to live without the hierarchy (that is, without calling
#   HierarchyUtilities::lane_info). For this to work, the keys 
#   'fa_ref','snps','genotype','bam' must be set.
#
sub check_genotype
{
    my ($self,$lane_path,$lock_file) = @_;

    my $needs_lane_info = 0;
    my $options = {};
    for my $key ('fa_ref','snps','genotype','bam')
    {
        if ( exists($$self{$key}) ) { $$options{$key} = $$self{$key} }
        else { $needs_lane_info = 1; print STDERR "no key $key\n"; }
    }

    if ( $needs_lane_info )
    {
        my $lane_info  = HierarchyUtilities::lane_info($lane_path);
        $$options{'fa_ref'}   = exists($$self{'fa_ref'}) ? $$self{'fa_ref'} : $$lane_info{'fa_ref'};
        $$options{'snps'}     = exists($$self{'snps'}) ? $$self{'snps'} : $$lane_info{'snps'};
        $$options{'genotype'} = exists($$lane_info{'genotype'}) ? $$lane_info{'genotype'} : '';
        $$options{'bam'}      = exists($$self{'bam'}) ? $$self{'bam'} : "$$lane_info{lane}.bam";
    }

    $$options{'bam'}           = "$lane_path/$$options{'bam'}";
    $$options{'bsub_opts'}     = $$self{'bsub_opts'};
    $$options{'glf'}           = $$self{'glf'};
    $$options{'samtools'}      = $$self{'samtools'};
    $$options{'min_glf_ratio'} = $$self{'min_glf_ratio'};
    $$options{'prefix'}        = $$self{'prefix'};
    $$options{'lock_file'}     = $lock_file;

    my $gtc = VertRes::Utils::GTypeCheck->new(%$options);
    $gtc->check_genotype();

    return $$self{'No'};
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

