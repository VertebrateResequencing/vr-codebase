=head1 NAME

VertRes::Utils::Mappers::bwa - mapping utility functions, bwa-specific

=head1 SYNOPSIS

use VertRes::Utils::Mappers::bwa;

my $mapping_util = VertRes::Utils::Mappers::bwa->new();

# use any of the utility functions described here, eg.
$mapping_util->do_mapping(ref => 'ref.fa',
                          read1 => 'reads_1.fastq',
                          read2 => 'reads_2.fastq',
                          output => 'output.sam',
                          insert_size => 2000);

=head1 DESCRIPTION

bwa-specific mapping functions, for selexa (illumina) lanes.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Utils::Mappers::bwa;

use strict;
use warnings;
use VertRes::Wrapper::bwa;

use base qw(VertRes::Utils::Mapping);

our %do_mapping_args = (insert_size => 'sampe_a');


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Mappers::bwa->new();
 Function: Create a new VertRes::Utils::Mappers::bwa object.
 Returns : VertRes::Utils::Mappers::bwa object
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

sub _bsub_opts {
    my ($self, $lane_path, $action) = @_;
    
    my %bsub_opts = (bsub_opts => '');
    
    if ($action eq 'map') {
        $bsub_opts{bsub_opts} = '-q normal -M5000000 -R \'select[mem>5000] rusage[mem=5000]\'';
    }
    else {
        return $self->SUPER::_bsub_opts($lane_path, $action);
    }
    
    return \%bsub_opts;
}

=head2 wrapper

 Title   : wrapper
 Usage   : my $wrapper = $obj->wrapper();
 Function: Get a bwa wrapper to actually do some mapping with.
 Returns : VertRes::Wrapper::bwa object (call do_mapping() on it)
 Args    : n/a

=cut

sub wrapper {
    my $self = shift;
    return VertRes::Wrapper::bwa->new(verbose => $self->verbose);
}

=head2 split_fastq

 Title   : split_fastq
 Usage   : $obj->split_fastq(read1 => 'reads_1.fastq',
                             read2 => 'reads_2.fastq',
                             split_dir => '/path/to/desired/split_dir',
                             chunk_size => 1000000);
 Function: Split the fastq(s) into multiple smaller files. This is just a
           convienience alias to VertRes::Utils::FastQ::split, with syntax
           more similar to do_mapping().
 Returns : int (the number of splits created)
 Args    : read1 => 'reads_1.fastq', read2 => 'reads_2.fastq'
           -or-
           read0 => 'reads.fastq'

           split_dir => '/path/to/desired/split_dir'
           chunk_size => int (max number of bases per chunk, default 1000000)

=cut

sub split_fastq {
    my ($self, %args) = @_;
    unless ($args{chunk_size}) {
        $args{chunk_size} = 1000000;
    }
    
    return $self->SUPER::split_fastq(%args);
}

=head2 do_mapping

 Title   : do_mapping
 Usage   : $obj->do_mapping(ref => 'ref.fa',
                            read1 => 'reads_1.fastq',
                            read2 => 'reads_2.fastq',
                            output => 'output.sam',
                            insert_size => 2000);
 Function: A convienience method that calls do_mapping() on the return value of
           wrapper(), translating generic options to those suitable for the
           wrapper.
 Returns : boolean (true on success)
 Args    : required options:
           ref => 'ref.fa'
           output => 'output.sam'

           read1 => 'reads_1.fastq', read2 => 'reads_2.fastq'
           -or-
           read0 => 'reads.fastq'

           and optional generic options:
           insert_size => int (default 2000)

=cut

sub do_mapping {
    my $self = shift;
    
    my @args = $self->_do_mapping_args(\%do_mapping_args, @_);
    
    my $wrapper = $self->wrapper;
    $wrapper->do_mapping(@args);
    
    # bwa directly generates sam files, so nothing futher to do
    
    return $wrapper->run_status >= 1;
}

1;
