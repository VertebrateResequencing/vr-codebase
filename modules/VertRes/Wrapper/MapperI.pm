=head1 NAME

VertRes::Wrapper::MapperI - interface common to mapper wrappers

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

[stub]

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::MapperI;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;
use Time::Format;
use VertRes::Utils::Sam;
use VertRes::Parser::sam;

use base qw(VertRes::Wrapper::WrapperI);


=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::MapperI->new();
 Function: Create a VertRes::Wrapper::MapperI object.
 Returns : VertRes::Wrapper::MapperI object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
}

=head2 version

 Title   : version
 Usage   : my $version = $obj->version();
 Function: Returns the program version.
 Returns : string representing version of the program 
 Args    : n/a

=cut

sub version {
    return 0;
}

=head2 setup_reference

 Title   : setup_reference
 Usage   : $obj->setup_reference($ref_fasta);
 Function: Do whatever needs to be done with the reference to allow mapping.
 Returns : boolean
 Args    : n/a

=cut

sub setup_reference {
    my ($self, $ref) = @_;
    return 1;
}

=head2 setup_fastqs

 Title   : setup_fastqs
 Usage   : $obj->setup_fastqs($ref_fasta, @fastqs);
 Function: Do whatever needs to be done with the fastqs to allow mapping.
 Returns : boolean
 Args    : n/a

=cut

sub setup_fastqs {
    my ($self, $ref, @fqs) = @_;
    return 1;
}

=head2 generate_sam

 Title   : generate_sam
 Usage   : $obj->generate_sam($out_sam, $ref_fasta, @fastqs);
 Function: Do whatever needs to be done with the reference and fastqs to
           complete mapping and generate a sam/bam file.
 Returns : boolean
 Args    : n/a

=cut

sub generate_sam {
    my ($self, $out, $ref, @fqs) = @_;
    return 0;
}

=head2 add_unmapped

 Title   : add_unmapped
 Usage   : $obj->add_unmapped($sam_file, @fastqs);
 Function: Do whatever needs to be done with the sam file to add in unmapped
           reads.
 Returns : boolean
 Args    : n/a

=cut

sub add_unmapped {
    my $self = shift;
    my $su = VertRes::Utils::Sam->new();
    return $su->add_unmapped(@_);
}

=head2 add_unmapped_from_fastq

 Title   : add_unmapped_from_fastq
 Usage   : $obj->add_unmapped_from_fastq($sam_file, $unmapped_fastq);
 Function: Do whatever needs to be done with the sam file to add in unmapped
           reads, found in a fastq containing only unmapped reads.
 Returns : boolean
 Args    : n/a

=cut

sub add_unmapped_from_fastq {
    my ($self, $sam, $unmapped) = @_;
    
    if (-s $unmapped) {
        open(my $sfh, '>>', $sam) || $self->throw("Could not append to $sam");
        my $fqp = VertRes::Parser::fastq->new(file => $unmapped);
        my $rh = $fqp->result_holder;
        
        # normally we'd work out the correct sam flags based on if the mate
        # was also unmapped, but since that requires we know what all the
        # unmapped reads are beforehand (difficult or high memory), screw it;
        # for mapper comparison the flags don't have to be perfect - just has
        # to have unmapped flag
        my $su = VertRes::Utils::Sam->new();
        my $flags = $su->calculate_flag(self_unmapped => 1);
        
        while ($fqp->next_result) {
            print $sfh "$rh->[0]\t$flags\t*\t0\t0\t*\t*\t0\t0\t$rh->[1]\t$rh->[2]\n";
        }
        close($sfh);
        
        return 1;
    }
    else {
        return 1;
    }
}

=head2 add_unmapped_from_fasta

 Title   : add_unmapped_from_fasta
 Usage   : $obj->add_unmapped_from_fasta($sam_file, $unmapped_fasta, @fastqs);
 Function: Do whatever needs to be done with the sam file to add in unmapped
           reads, found in a fasta containing only unmapped reads.
 Returns : boolean
 Args    : n/a

=cut

sub add_unmapped_from_fasta {
    my ($self, $sam, $unmapped, @fastqs) = @_;
    
    if (-s $unmapped) {
        open(my $sfh, '>>', $sam) || $self->throw("Could not append to $sam");
        # need a fasta parser...
        open(my $ufh, $unmapped) || $self->throw("Could not open $unmapped");
        
        # normally we'd work out the correct sam flags based on if the mate
        # was also unmapped, but since that requires we know what all the
        # unmapped reads are beforehand (difficult or high memory), screw it;
        # for mapper comparison the flags don't have to be perfect - just has
        # to have unmapped flag
        my $su = VertRes::Utils::Sam->new();
        my $flags = $su->calculate_flag(self_unmapped => 1);
        
        # and we should use the fasta to find the unmapped ids, then go through
        # the original whole fastqs to get the quality strings... but screw it;
        # for mapper comparison the quality string doesn't matter
        
        my $id;
        while (<$ufh>) {
            chomp;
            if (/^>(.+)/) {
                $id = $1;
            }
            else {
                my $seq = $_;
                my $length = length $seq;
                my $qual = '?' x $length;
                print $sfh "$id\t$flags\t*\t0\t0\t*\t*\t0\t0\t$seq\t$qual\n";
            }
        }
        close($sfh);
        
        return 1;
    }
    else {
        return 1;
    }
}

=head2 do_mapping

 Title   : do_mapping
 Usage   : $wrapper->do_mapping(ref => 'ref.fa',
                                read1 => 'reads_1.fastq',
                                read2 => 'reads_2.fastq',
                                output => 'output.sam');
 Function: Run mapper on the supplied files, generating a sam file of the
           mapping. Checks the sam file isn't truncated.
 Returns : n/a
 Args    : required options:
           ref => 'ref.fa'
           output => 'output.sam'

           read1 => 'reads_1.fastq', read2 => 'reads_2.fastq'
           -or-
           read0 => 'reads.fastq'

=cut

sub do_mapping {
    my ($self, %args) = @_;
    
    my $orig_run_method = $self->run_method;
    $self->run_method('system');
    
    my $ref_fa = delete $args{ref} || $self->throw("ref is required");
    my @fqs;
    if (defined $args{read1} && defined $args{read2} && ! defined $args{read0}) {
        push(@fqs, delete $args{read1});
        push(@fqs, delete $args{read2});
    }
    elsif (defined $args{read0} && ! defined $args{read1} && ! defined $args{read2}) {
        push(@fqs, delete $args{read0});
    }
    else {
        $self->throw("Bad read args: need read1 and read2, or read0");
    }
    my $out_sam = delete $args{output};
    
    # setup reference-related files
    $self->setup_reference($ref_fa) || $self->throw("failed during the reference step");
    
    # setup fastq-related files (may involve alignment)
    $self->setup_fastqs($ref_fa, @fqs) || $self->throw("failed during the fastq step");
    
    # run the alignment/ combine alignments into final sam file
    unless (-s $out_sam) {
        my $tmp_sam = $out_sam.'_tmp';
        
        if (@fqs == 2) {
            unless (-s $tmp_sam) {
               $self->generate_sam($tmp_sam, $ref_fa, @fqs, %args) || $self->throw("failed during the alignment step");
            }
        }
        else {
            unless (-s $tmp_sam) {
               $self->generate_sam($tmp_sam, $ref_fa, $fqs[0], undef, %args) || $self->throw("failed during the alignment step");
            }
        }
        
        # add in unmapped reads if necessary
        $self->add_unmapped($tmp_sam, @fqs) || $self->throw("failed during the add unmapped step");
        
        # check the sam file isn't truncted
        my $num_reads = 0;
        my $io = VertRes::IO->new();
        foreach my $fastq_file (@fqs) {
            $io->file($fastq_file);
            if ($fastq_file =~ /\.fastq|\.fq/) {
                my $fastq_lines = $io->num_lines();
                $num_reads += $fastq_lines / 4;
            }
            else {
                # blindly assume it's a fasta file
                my $fh = $io->fh;
                while (<$fh>) {
                    $num_reads++ if /^>/;
                }
                close($fh);
            }
        }
        $io->file($tmp_sam);
        my $sam_lines = $io->num_lines();
        
        if ($sam_lines >= $num_reads) {
            move($tmp_sam, $out_sam) || $self->throw("Failed to move $tmp_sam to $out_sam: $!");
        }
        else {
            $self->warn("a sam file was made by do_mapping(), but it only had $sam_lines lines, not the expected $num_reads...");
            $self->_set_run_status(-1);
            $self->warn("... moving it to $out_sam.bad");
            move($tmp_sam, "$out_sam.bad");
        }
    }
    else {
        $self->_set_run_status(1);
    }
    
    $self->run_method($orig_run_method);
    return;
}

=head2 run

 Title   : run
 Usage   : Do not call directly: use one of the other methods instead.
 Function: n/a
 Returns : n/a
 Args    : paths to input/output files

=cut

sub simple_run {
    my ($self, $cmd) = @_;
    my $exe = $self->exe;
    $cmd = $exe.' '.$cmd;
    $self->debug("[$time{'yyyy/mm/dd hh:mm:ss'}] will run command '$cmd'");
    system($cmd) && $self->throw("command '$cmd' failed");
}

1;
