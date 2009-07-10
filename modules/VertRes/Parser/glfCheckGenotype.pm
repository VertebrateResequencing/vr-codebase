=head1 NAME

VertRes::Parser::glfCheckGenotype - parse checkGenotype output of glf

=head1 SYNOPSIS

use VertRes::Parser::glfCheckGenotype;

# create object, supplying glf checkGenotype output file or filehandle
my $pars = VertRes::Parser::glfCheckGenotype->new(file => 'glf.out');

# get basic summary info:
my $entropy = $pars->entropy();

# get the array reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

my $entropy = $pars->entropy;

# loop through the output, getting results
while ($pars->next_result()) {
    my $individual = $result_holder->[0];
    my $likelihood = $result_holder->[1];
    my $num_of_sites = $result_holder->[2];
    my $avg_depth = $result_holder->[3];
}

=head1 DESCRIPTION

A straightforward parser for checkGenotype output of glf.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::glfCheckGenotype;

use strict;
use warnings;

use base qw(VertRes::Parser::ParserI);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::glfCheckGenotype->new(file => 'filename');
 Function: Build a new VertRes::Parser::glfCheckGenotype object.
 Returns : VertRes::Parser::glfCheckGenotype object
 Args    : file => filename -or- fh => filehandle

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 entropy

 Title   : entropy
 Usage   : my $entropy = $obj->entropy();
 Function: Get the entropy.
 Returns : number
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub entropy {
    my $self = shift;
    
    $self->_get_header;
    unless (defined $self->{_entropy}) {
        $self->throw("entropy unknown - did you supply a file(/handle) and was it the output of 'glf checkGentotype'?");
    }
    
    return $self->{_entropy};
}

sub _get_header {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    return 1 if $self->{'_got_header'.$fh_id};
    
    my $saw = 0;
    while (<$fh>) {
        if (/^entropy (\d+\.\d+)$/) {
            $self->{_entropy} = $1;
            $saw++;
            last;
        }
    }
    
    if ($saw == 1) {
        $self->{'_got_header'.$fh_id} = 1;
        return 1;
    }
    return;
}

=head2 result_holder

 Title   : result_holder
 Usage   : my $result_holder = $obj->result_holder()
 Function: Get the data structure that will hold the last result requested by
           next_result()
 Returns : array ref, where the elements are:
           [0] the snp file (individual) being compared against
           [1] negative log-likelihood that the glf being checked has the same
               genotype as the individual in [0] (closer to 0 is more likely)
           [2] number of sites [1] was calculated over
           [3] average depth
 Args    : n/a

=cut

=head2 next_result

 Title   : next_result
 Usage   : while ($obj->next_result()) { # look in result_holder }
 Function: Parse the next base position line from the fastqcheck output.
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    my $self = shift;
    
    # just return if no file set
    my $fh = $self->fh() || return;
    
    # make sure we've gotten our header first
    $self->_get_header() || $self->throw("Unable to parse header before first result - is this the output of 'glf checkGentotype'?");
    
    # get the next line
    my $line = <$fh> || return;
    
    # sample /nfs/users/nfs_p/pd3/g1k-sandbox/test-genotypes/snps-hapmap3/NA19468.snp likelihood 341694 over 76226 sites, avg depth 0.101695
    my @data = $line =~ /^sample (\S+) likelihood (\d+) over (\d+) sites, avg depth (\S+)$/;
    @data == 4 or return;
    
    for my $i (0..$#data) {
        $self->{_result_holder}->[$i] = $data[$i];
    }
    
    return 1;
}

1;
