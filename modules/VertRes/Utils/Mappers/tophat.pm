=head1 NAME

VertRes::Utils::Mappers::tophat - mapping utility functions

=head1 SYNOPSIS

use VertRes::Utils::Mappers::tophat;

my $mapping_util = VertRes::Utils::Mappers::tophat->new();

# use any of the utility functions described here, eg.
$mapping_util->do_mapping(ref => 'ref.fa',
                          read1 => 'reads_1.fastq',
                          read2 => 'reads_2.fastq',
                          output => 'output.sam',
                          );

=head1 DESCRIPTION

tophat-specific mapping functions, for selexa (illumina) lanes.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk
Thomas Keane: tk2@sanger.ac.uk

=cut

package VertRes::Utils::Mappers::tophat;

use strict;
use warnings;
use VertRes::Wrapper::tophat;

use base qw(VertRes::Utils::Mapping);


our %do_mapping_args = (insert_size => 'i');


=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Mappers::tophat->new();
 Function: Create a new VertRes::Utils::Mappers::tophat object.
 Returns : VertRes::Utils::Mappers::tophat object
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
        $bsub_opts{bsub_opts} = '-q long -M12000000 -R \'select[mem>12000] rusage[mem=12000]\'';
    }
    else {
        return $self->SUPER::_bsub_opts($lane_path, $action);
    }
    
    return \%bsub_opts;
}

=head2 wrapper

 Title   : wrapper
 Usage   : my $wrapper = $obj->wrapper();
 Function: Get a tophat wrapper to actually do some mapping with.
 Returns : VertRes::Wrapper::tophat object (call do_mapping() on it)
 Args    : n/a

=cut

sub wrapper {
    my $self = shift;
    my $exe = $self->{exe} || 'tophat';
    return VertRes::Wrapper::tophat->new(verbose => $self->verbose, exe => $exe);
}

=head2 name

 Title   : name
 Usage   : my $name = $obj->name();
 Function: Returns the program name.
 Returns : string representing name of the program 
 Args    : n/a

=cut

sub name {
    my $self = shift;
    return 'tophat';
}

=head2 split_fastq

 Title   : split_fastq
 Usage   : $obj->split_fastq(read1 => 'reads_1.fastq',
                             read2 => 'reads_2.fastq',
                             split_dir => '/path/to/desired/split_dir',
                             chunk_size => 1000000000);
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
        $args{chunk_size} = 1000000000;
    }
    
    return $self->SUPER::split_fastq(%args);
}

=head2 do_mapping

 Title   : do_mapping
 Usage   : $obj->do_mapping(ref => 'ref.fa',
                            read1 => 'reads_1.fastq',
                            read2 => 'reads_2.fastq',
                            output => 'output.sam',
                            );
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
=cut

sub do_mapping {
    my ($self, %input_args) = @_;
    
   # if (defined $input_args{insert_size} && $input_args{insert_size} != 2000) {
   #     # in samlt, the insert_size parameter is the maximum, not the average,
   #     # so we multiply by 3
   #     $input_args{insert_size} *= 3;
   # }
    
    my %args = $self->_do_mapping_args(\%do_mapping_args, %input_args);
    $args{insert_size} = $input_args{insert_size} if(defined($input_args{insert_size} ));
    
    if(defined($input_args{max_intron_length} ))
    {
      $args{max_intron_length} = $input_args{max_intron_length} ;
    }
    else
    {
        $args{max_intron_length} = 10000; 
    }
    
    my $wrapper = $self->wrapper;
    $wrapper->do_mapping(%args);
    $self->_add_command_line($wrapper->command_line());

    # tophat directly generates sam files, so nothing futher to do
    
    return $wrapper->run_status >= 1;
}

=head2 check_total_reads

 Title   : check_total_reads
 Usage   : $self->check_total_reads();
 Function: Flag which tells the pipeline if it should check the total number of input reads against the 
           number of reads in the output file. This check needs to be disabled for some specialised mappers 
           (e.g. tophat) which dont include the unmapped reads in the output file
 Example : $self->check_total_reads();
 Returns : 0 or 1
 Args    : n/a

=cut

sub  check_total_reads
{
   my $self = shift;
   return 0;
}

1;
