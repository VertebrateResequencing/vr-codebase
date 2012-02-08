=head1 NAME

VertRes::Parser::bamcheck - parse bamcheck files

=head1 SYNOPSIS

    use VertRes::Parser::bamcheck;

    # create object, supplying bamcheck output file or filehandle
    my $pars = VertRes::Parser::bamcheck->new(file => 'file.bam.bc');
    
    # get basic summary info:
    my $total_length_of_all_sequences = $pars->get('total_length');
    my $average_length_of_a_sequence  = $pars->get('avg_length');
    my $length_of_longest_sequence    = $pars->get('max_length');
    my $num_of_sequences = $pars->get('sequences');
    my $average_quality  = $pars->get('avg_qual');
    my $unmapped_reads   = $pars->get('reads_unmapped');
    my $is_sorted        = $pars->get('is_sorted');
    my $is_paired        = $pars->get('is_paired');

=head1 DESCRIPTION

A straightforward parser for bamcheck files.

=head1 AUTHOR

petr.danecek@sanger

=cut

package VertRes::Parser::bamcheck;

use strict;
use warnings;

use base qw(VertRes::Parser::ParserI);


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::bamcheck->new(file => 'filename');
 Function: Build a new VertRes::Parser::bamcheck object.
 Returns : VertRes::Parser::bamcheck object
 Args    : file => filename -or- fh => filehandle

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 get

 Title   : get
 Usage   : my $num_of_sequences = $obj->get('sequences');
 Function: Get the number of sequences that the bamcheck file is summarising.
 Returns : Numeric value or array pointer
 Args    : One of the following numeric values:
                '1st_fragments',
                'avg_insert_size',
                'avg_length',
                'avg_qual',
                'bases_duplicated',
                'bases_mapped',
                'bases_mapped_cigar',
                'bases_trimmed',
                'error_rate',
                'inward_pairs',
                'is_paired',
                'is_sorted',
                'last_fragments',
                'max_length',
                'mismatches',
                'other_pairs',
                'outward_pairs',
                'reads_duplicated',
                'reads_mapped',
                'reads_mq0',
                'reads_unmapped',
                'reads_unpaired',
                'reads_paired,
                'sd_insert_size',
                'sequences',
                'total_length',

            or one of the arrays:
                'coverage',
                'first_fragment_gc',
                'first_fragment_qualities',
                'gc_depth',
                'gc_content_per_cycle',
                'indel_cycles',
                'indel_dist',
                'insert_size',
                'last_fragment_gc',
                'last_fragment_qualities',
                'mismatches_per_cycle',
                'read_lengths',

=cut

sub get
{
    my ($self,$key) = @_;
    return $self->_getter($key);
}


sub _getter 
{
    my ($self, $var) = @_;
    
    $self->_get_header;
    unless (defined $self->{"_$var"}) {
        $self->throw("$var unknown - did you supply a file and was it a bamcheck file?");
    }
    return $self->{"_$var"};
}

sub _get_header {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    return 1 if $self->_header_parsed();

    my %sn_mapping = 
    ( 
        'sequences'         => 'sequences',
        '1st fragments'     => '1st_fragments',
        'last fragments'    => 'last_fragments',
        'reads mapped'      => 'reads_mapped',
        'reads MQ0'         => 'reads_mq0',
        'reads unmapped'    => 'reads_unmapped',
        'reads unpaired'    => 'reads_unpaired',
        'reads paired'      => 'reads_paired',
        'reads duplicated'  => 'reads_duplicated',
        'total length'      => 'total_length',
        'bases mapped'      => 'bases_mapped',
        'bases mapped (cigar)' => 'bases_mapped_cigar',
        'bases trimmed'     => 'bases_trimmed',
        'bases duplicated'  => 'bases_duplicated',
        'mismatches'        => 'mismatches',
        'error rate'        => 'error_rate',
        'average length'    => 'avg_length',
        'maximum length'    => 'max_length',
        'average quality'   => 'avg_qual',
        'insert size average' => 'avg_insert_size',
        'insert size standard deviation' => 'sd_insert_size',
        'is sorted'         => 'is_sorted',
        'is paired'         => 'is_paired',
        'inward oriented pairs'         => 'inward_pairs',
        'outward oriented pairs'        => 'outward_pairs',
        'pairs with other orientation'  => 'other_pairs',
    );
    my %mapping = 
    (
        'FFQ'   => 'first_fragment_qualities',
        'LFQ'   => 'last_fragment_qualities',
        'GCF'   => 'first_fragment_gc',
        'GCL'   => 'last_fragment_gc',
        'IC'    => 'indel_cycles',
        'ID'    => 'indel_dist',
        'IS'    => 'insert_size',
        'GCD'   => 'gc_depth',
        'COV'   => 'coverage',
        'MPC'   => 'mismatches_per_cycle',
        'GCC'   => 'gc_content_per_cycle',
        'RL'    => 'read_lengths',
    );
    
    while (my $line=<$fh>) 
    {
        if ( $line=~/^#/ ) { next; }

        #   SN  sequences:  73509458
        #   SN  1st fragments:  36754729
        #   ...
        if ( $line=~/^SN\t/ )
        {
            if ( !($line=~/^SN\t([^:]+):\t(\S+)/) ) { $self->throw("Could not parse: $line"); }
            if ( !exists($sn_mapping{$1}) ) { $self->throw("FIXME: no mapping for [$1] [$line]?\n"); }
            my $key = '_' . $sn_mapping{$1};
            $$self{$key} = $2;
            next;
        }

        my ($key,@items) = split(/\t/,$line);
        if ( !exists($mapping{$key}) ) { $self->throw("FIXME: no mapping for [$key] [$line]?"); }
        $key = '_'.$mapping{$key};

        chomp($items[-1]);
        push @{$$self{$key}}, \@items;
    }
    
    $self->_set_header_parsed();
    return 1;
}

=head2 result_holder

 Title   : result_holder
 Usage   : my $result_holder = $obj->result_holder()
 Function: Get the data structure that will hold the last result requested by
           next_result()
 Returns : array ref, where the elements are:
           [0] base position (starting at position 0 == 'total')
           [6] count of bases at position [0] with quality 0
           [7] count of bases at position [0] with quality 1
           ... etc. for quality up to the highest quality in the bam file
 Args    : n/a

=cut

=head2 next_result

 Title   : next_result
 Usage   : while ($obj->next_result()) { # look in result_holder }
 Function: Parse the next base position line from the bamcheck output.
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    
    $self->_get_header() || $self->throw("Unable to parse header before first result - is this a bamcheck file?");
    
    # get the next line
    my $line = <$fh> || return;
    
    my @data = split(qr/\s+/, $line);
    @data || return;

    chomp($data[-1]);
    $self->{_result_holder} = \@data;
       
    return 1;
}

1;
