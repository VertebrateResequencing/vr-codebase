=head1 NAME

VertRes::Parser::flagstat - parse samtools flagstat output

=head1 SYNOPSIS

use VertRes::Parser::flagstat;

# create object, supplying samtools flagstat output or filehandle
my $pars = VertRes::Parser::flagstat->new(file => 'bam.flagstat');

# get the result of a particular stat
my $mapped_reads = $pars->mapped_reads;
# etc.

=head1 DESCRIPTION

A straightforward parser for samtools flagstat output

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::flagstat;

use strict;
use warnings;

use base qw(VertRes::Parser::ParserI);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::flagstat->new(file => 'filename');
 Function: Build a new VertRes::Parser::flagstat object.
 Returns : VertRes::Parser::flagstat object
 Args    : file => filename -or- fh => filehandle

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}


=head2 total_reads

 Title   : total_reads
 Usage   : my $total_reads = $obj->total_reads();
 Function: Get the number of reads in the bam file the flagstat is summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub total_reads {
    my $self = shift;
    return $self->_getter('total_reads');
}

=head2 qc_failures

 Title   : qc_failures
 Usage   : my $qc_failures = $obj->qc_failures();
 Function: Get the number of qc failures in the bam file the flagstat is
           summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub qc_failures {
    my $self = shift;
    return $self->_getter('qc_failures');
}

=head2 duplicates

 Title   : duplicates
 Usage   : my $duplicates = $obj->duplicates();
 Function: Get the number of duplicates in the bam file the flagstat is
           summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub duplicates {
    my $self = shift;
    return $self->_getter('duplicates');
}

=head2 mapped_reads

 Title   : mapped_reads
 Usage   : my $mapped_reads = $obj->mapped_reads();
 Function: Get the number of mapped reads in the bam file the flagstat is
           summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub mapped_reads {
    my $self = shift;
    return $self->_getter('mapped_reads');
}

=head2 paired_reads

 Title   : paired_reads
 Usage   : my $paired_reads = $obj->paired_reads();
 Function: Get the number of paired reads in the bam file the flagstat is
           summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub paired_reads {
    my $self = shift;
    return $self->_getter('paired_reads');
}

=head2 read1_reads

 Title   : read1_reads
 Usage   : my $read1_reads = $obj->read1_reads();
 Function: Get the number of forward reads in the bam file the flagstat is
           summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub read1_reads {
    my $self = shift;
    return $self->_getter('read1_reads');
}

=head2 read2_reads

 Title   : read2_reads
 Usage   : my $read2_reads = $obj->read2_reads();
 Function: Get the number of reverse reads in the bam file the flagstat is
           summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub read2_reads {
    my $self = shift;
    return $self->_getter('read2_reads');
}

=head2 mapped_proper_paired_reads

 Title   : mapped_proper_paired_reads
 Usage   : my $mapped_proper_paired_reads = $obj->mapped_proper_paired_reads();
 Function: Get the number of mapped, properly paired reads in the bam file the
           flagstat is summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub mapped_proper_paired_reads {
    my $self = shift;
    return $self->_getter('mapped_proper_paired_reads');
}

=head2 mapped_paired_reads

 Title   : mapped_paired_reads
 Usage   : my $mapped_paired_reads = $obj->mapped_paired_reads();
 Function: Get the number of mapped, paired reads in the bam file the flagstat
           is summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub mapped_paired_reads {
    my $self = shift;
    return $self->_getter('mapped_paired_reads');
}

=head2 singletons

 Title   : singletons
 Usage   : my $singletons = $obj->singletons();
 Function: Get the number of singletons in the bam file the flagstat is
           summarising.
 Returns : int
 Args    : n/a (new() or file() must allready have been supplied with a
                filename or filehandle)

=cut

sub singletons {
    my $self = shift;
    return $self->_getter('singletons');
}

=head2 mate_mapped_to_other_chr

 Title   : mate_mapped_to_other_chr
 Usage   : my $mate_mapped_to_other_chr = $obj->mate_mapped_to_other_chr();
 Function: Get the number of reads where the mate is mapped to another chr in
           the bam file the flagstat is summarising.
 Returns : int
 Args    : boolean: false/undef for all, true to only report on mates mapped
           with Q 5 or higher.

=cut

sub mate_mapped_to_other_chr {
    my $self = shift;
    my $q5 = shift;
    if ($q5) {
        return $self->_getter('mate_mapped_to_other_chr_q5');
    }
    else {
       return $self->_getter('mate_mapped_to_other_chr');
    }
}

sub _getter {
    my ($self, $var) = @_;
    
    $self->_get_header || return;
    my $fh_id = $self->_fh_id;
    my $data = $self->{'_data'.$fh_id};
    unless (defined $data->{"_$var"}) {
        $self->throw("$var unknown - did you supply a file and was it the output of samtools flagstat?");
    }
    return $data->{"_$var"};
}

sub _get_header {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    return 1 if $self->{'_got_header'.$fh_id};
    
    my %data;
    
    while (<$fh>) {
        if (/^(\d+) in total$/) {
            $data{_total_reads} = $1;
        }
        elsif (/^(\d+) QC failure$/) {
            $data{_qc_failures} = $1;
        }
        elsif (/^(\d+) duplicates$/) {
            $data{_duplicates} = $1;
        }
        elsif (/^(\d+) mapped/) {
            $data{_mapped_reads} = $1;
        }
        elsif (/^(\d+) paired in sequencing$/) {
            $data{_paired_reads} = $1;
        }
        elsif (/^(\d+) read1$/) {
            $data{_read1_reads} = $1;
        }
        elsif (/^(\d+) read2$/) {
            $data{_read2_reads} = $1;
        }
        elsif (/^(\d+) properly paired/) {
            $data{_mapped_proper_paired_reads} = $1;
        }
        elsif (/^(\d+) with itself and mate mapped$/) {
            $data{_mapped_paired_reads} = $1;
        }
        elsif (/^(\d+) singletons/) {
            $data{_singletons} = $1;
        }
        elsif (/^(\d+) with mate mapped to a different chr$/) {
            $data{_mate_mapped_to_other_chr} = $1;
        }
        elsif (/^(\d+) with mate mapped to a different chr \(/) {
            $data{_mate_mapped_to_other_chr_q5} = $1;
        }
    }
    
    if (keys %data >= 11) {
        $self->{'_got_header'.$fh_id} = 1;
        $self->{'_data'.$fh_id} = \%data;
        return 1;
    }
    return;
}

=head2 result_holder

 Title   : result_holder
 Usage   : my $result_holder = $obj->result_holder()
 Function: Get the data structure that will hold the last result requested by
           next_result()
 Returns : array ref, where the array is empty
 Args    : n/a

=cut

=head2 next_result

 Title   : next_result
 Usage   : n/a; no need to use this
 Function: n/a
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    return;
}

1;
