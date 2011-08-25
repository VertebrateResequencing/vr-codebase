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
        $bsub_opts{bsub_opts} = '-q long -M7900000 -R \'select[mem>7900] rusage[mem=7900]\'';
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
    my $exe = $self->{exe} || 'bwa';
    return VertRes::Wrapper::bwa->new(verbose => $self->verbose, exe => $exe);
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
    return 'bwa';
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
                            insert_size => 2000);
 Function: A convienience method that calls do_mapping() on the return value of
           wrapper(), translating generic options to those suitable for the
           wrapper. insert_size here is the expected/average insert size and
           will be ajusted as appropriate for bwa.
 Returns : boolean (true on success)
 Args    : required options:
           ref => 'ref.fa'
           output => 'output.sam'

           read1 => 'reads_1.fastq', read2 => 'reads_2.fastq'
           -or-
           read0 => 'reads.fastq'

           optionally, to check an error output file for common problems and
           delete certain files before reattempting:
           error_file => '/path/to/STDERR/output/of/a/previous/call'
           
           to set the aln q parameter (quality threshold for read trimming):
           aln_q => int (default 15)

           and optional generic options:
           insert_size => int (default 2000)

=cut

sub do_mapping {
    my ($self, %input_args) = @_;
    
    if (defined $input_args{insert_size} && $input_args{insert_size} != 2000) {
        # in bwa, the insert_size parameter is the maximum, not the average,
        # so we multiply by 3
        $input_args{insert_size} *= 3;
    }
    
    my @args = $self->_do_mapping_args(\%do_mapping_args, %input_args);
    
    my $error_file = delete $input_args{error_file};
    if ($error_file && -s $error_file && ! -s $input_args{output}) {
        # check the error file for common problems and delete various files as
        # necessary before reattempting
        open(my $efh, $error_file) || $self->throw("Could not open error file '$error_file'");
        my $delete_sais = 0;
        while (<$efh>) {
            if (/weird pairing/) {
                $delete_sais = 1;
                last;
            }
        }
        close($efh);
        
        if ($delete_sais) {
            my %args = @args;
            while (my ($key, $val)  = each %args) {
                if ($key =~ /^read/) {
                    my $sai = $val;
                    $sai =~ s/\.f[^.]+(?:\.gz)?$/.sai/;
                    unlink($sai);
                }
            }
            
            $self->warn("deleted existing sai files due to weird pairing in previous mapping attempt");
        }
    }
    
    my $aln_q = delete $input_args{aln_q};
    unless (defined $aln_q) {
        $aln_q = 15;
    }
    
    my $wrapper = $self->wrapper;
    $wrapper->do_mapping(@args, aln_q => $aln_q);
    push((@{$self->{command_list}}), @{$wrapper->{command_list}});
    
    # bwa directly generates sam files, so nothing futher to do
    
    return $wrapper->run_status >= 1;
}

1;
