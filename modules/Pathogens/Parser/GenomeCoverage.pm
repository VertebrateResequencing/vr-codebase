
=head1 NAME

GenomeCoverage.pm - Find genome coverage from bamcheck file and reference size.

=head1 SYNOPSIS

use GenomeCoverage;

my $gc = Pathogens::Parser::GenomeCoverage->new( bamcheck => 'file.bam.bc',
                                                 ref_size => 4646520);

my @bases_covered = $gc->coverage(1,2,5,10,20,50,100);
my ($genome_bases_covered, $depth_mean, $depth_sd) = $gc->coverage_depth();

=head1 DESCRIPTION

Find the amount of genome covered and coverage depth from a bamcheck file. 
The size of the reference genome is required for the mean coverage depth calculation.

=cut 


package Pathogens::Parser::GenomeCoverage;

use Moose;
use VertRes::Parser::bamcheck;
use VertRes::Utils::Math;
use base qw(VertRes::Base);

has 'bamcheck' => ( is => 'ro', isa => 'Str', required => 1 );
has 'ref_size' => ( is => 'rw', isa => 'Int', required => 0 );
has '_bc_coverage' => ( is => 'ro', isa => 'ArrayRef', lazy_build => 1 );

sub _build__bc_coverage
{
    my ($self) = @_;

    # get coverage from bamcheck file
    my $bc = VertRes::Parser::bamcheck->new( file => $self->bamcheck );
    my $bc_cov;
    eval{$bc_cov = $bc->get('coverage');};
    if($@){return [];}
    return $bc_cov;
}


=head2 coverage

 Title   : coverage
 Usage   : my @bases_covered = $obj->coverage(@coverage_depth);
 Function: Parses bamcheck file to find the number of bases covered to specified depths.
           If coverage depth is not supplied then function returns 1X coverage.
 Returns : array of int. 
 Args    : Optional list of coverage depths eg (1,2,5,10,20,50,100)

=cut

sub coverage
{
    my ($self,@bin) = @_;
    
    unless( -e $self->bamcheck ){ $self->throw("Input file not found: ".$self->bamcheck."\n"); }    
    unless(@bin){ @bin = (1); } # Default to 1X coverage.

    # Sort coverage depths into ascending order.
    # Zero base count.
    my @sorted_bin;
    my @coverage;
    for(my $i=0; $i < @bin; $i++)
    {
        $sorted_bin[$i][0] = $i;
        $sorted_bin[$i][1] = $bin[$i];
        $coverage[$i]      = 0;
    }
    @sorted_bin = sort {my @a = @{$a}; my @b = @{$b}; $a[1] <=> $b[1];} @sorted_bin;

    # Get coverage from bamcheck file
    my @bc_cover = @{$self->_bc_coverage};

    for(my $x=1; $x < @bc_cover; $x++)
    {
        for(my $i=0; $i < @bin; $i++)
        {
            if($sorted_bin[$i][1] <= $bc_cover[$x][1] ){ $coverage[$sorted_bin[$i][0]] += $bc_cover[$x][2]; }
            else{ last; }
        }
    }
    return @coverage;
}


=head2 coverage_depth

 Title   : coverage_depth
 Usage   : my($bases_covered, $depth, $depth_sd) = $obj->coverage_depth();
 Function: Parses bamcheck file to find the number of bases covered, depth of coverage and standard 
           deviation for depth of coverage. Reference size must be set for this function.
 Returns : array with number of bases covered, average depth of coverage and standard deviation for
           depth of coverage. 
 Args    : none

=cut

sub coverage_depth
{
    my($self) = @_;

    unless( -e $self->bamcheck ){ $self->throw("Input file not found: ".$self->bamcheck."\n"); }
    unless( defined $self->ref_size ){ $self->throw("Reference size must be set for mean coverage depth calculation.\n"); }
    unless( $self->ref_size =~ /^\d+$/ && $self->ref_size > 0 ){ $self->throw("Reference size must be non-zero integer\n"); }

    # Build depth histogram
    my $coverage = 0; # genome bases mapped
    my %depth_hist;   # Depth from bamcheck, key = coverage, value = total_bases.
    my @bc_cover = @{$self->_bc_coverage};

    for(my $x=1; $x < @bc_cover; $x++)
    {
        $coverage += $bc_cover[$x][2];
        $depth_hist{$bc_cover[$x][1]} += $bc_cover[$x][2] if $bc_cover[$x][2];
    }

    unless( $coverage <= $self->ref_size ){ $self->throw("Total bases found by bamcheck exceeds size of reference sequence.\n"); }

    # Add ummapped bases to histogram.
    $depth_hist{0} = $self->ref_size - $coverage;

    # Calculate Mean depth and SD
    my $math_util = VertRes::Utils::Math->new();
    my %stats = $math_util->histogram_stats(\%depth_hist);

    return($coverage, $stats{'mean'}, $stats{'standard_deviation'});
}

1;
