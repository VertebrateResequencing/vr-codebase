=head1 NAME

VertRes::Wrapper::tophat - wrapper for tophat

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

bowtie-build ref.fa out_base

tophat ref read_1.fastq read_2.fastq

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::tophat;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;
use VertRes::Parser::fastqcheck;
use VertRes::Wrapper::bowtie;

use base qw(VertRes::Wrapper::MapperI);

my %e_parameter = (
    37 	=> 70,
    63	=> 120,
    76  => 140,
    92 	=> 160,
    108 => 160,
    123	=> 200,
    156	=> 240
);

=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::bowtie->new();
 Function: Create a VertRes::Wrapper::bowtie object.
 Returns : VertRes::Wrapper::bowtie object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args, exe => '/software/pathogen/external/apps/usr/bin/tophat');
    
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
    my $self = shift;
    
    my $exe = $self->{orig_exe};
    open(my $fh, "$exe 2>&1 |") || $self->throw("Could not start $exe");
    my $version = 0;
    while (<$fh>) {
        if (/TopHat v(\S+)/) {
            $version = $1;
            last;
        }
    }
    close($fh);
    
    return $version;
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
    
    my $bowtie = VertRes::Wrapper::bowtie->new();
    return $bowtie->setup_reference($ref);
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
    
    foreach my $fq (@fqs) {
        if ($fq =~ /\.gz$/) {
            my $fq_new = $fq;
            $fq_new =~ s/\.gz$//;
            
            unless (-s $fq_new) {
                my $i = VertRes::IO->new(file => $fq);
                my $o = VertRes::IO->new(file => ">$fq_new");
                my $ifh = $i->fh;
                my $ofh = $o->fh;
                while (<$ifh>) {
                    print $ofh $_;
                }
                $i->close;
                $o->close;
            }
        }
    }
    
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
    my ($self, $out, $ref, @fqs,%other_args) = @_;
    
    unless (-s $out) {
        my $e = 70;
        
        my $longest_read = 0;
        foreach my $fq (@fqs) {
            my $pars = VertRes::Parser::fastqcheck->new(file => "$fq.fastqcheck");
            my $length = $pars->max_length();
            if ($length > $longest_read) {
                $longest_read = $length;
            }
        }
        
        foreach (sort { $a <=> $b } keys %e_parameter) {
            if ($longest_read <= $_) {
                $e = $e_parameter{$_};
            }
        }
        
        foreach my $fq (@fqs) {
            $fq =~ s/\.gz$//;
        }
        
        my $X = 1000;
        my $max_intron_length_str = $self->_max_intron_length_option(%other_args);
        my $inner_mate_str = $self->_insert_size_option($longest_read, %other_args);
        
        # parameters needed
        #-I to specify sensible max intron length
        #-r for inner-mate distance  ( insert size - readlength*2
        # index files
        $self->simple_run(" $max_intron_length_str $inner_mate_str $ref $fqs[0] $fqs[1]");
    }
    
    return -s $out ? 1 : 0;
}

sub _max_intron_length_option
{
  my ($self, %other_args) = @_;
  
  $other_args{max_intron_length}  ||=10000;
  return "-I ".$other_args{max_intron_length};
}

sub _insert_size_option
{
  my ($self, $longest_read, %other_args) = @_;
  
  $other_args{insert_size}  ||=500;
  my $inner_mate_distance = $other_args{insert_size} - (2*$longest_read);
  if($inner_mate_distance <0)
  {
    $inner_mate_distance = 50;
  }
  return "-r $inner_mate_distance";
}

=head2 add_unmapped

 Title   : add_unmapped
 Usage   : $obj->add_unmapped($sam_file, $ref_fasta, @fastqs);
 Function: Do whatever needs to be done with the sam file to add in unmapped
           reads.
 Returns : boolean
 Args    : n/a

=cut

sub add_unmapped {
    my ($self, $sam, $ref, @fqs) = @_;
    return 1;
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

=head2 run

 Title   : run
 Usage   : Do not call directly: use one of the other methods instead.
 Function: n/a
 Returns : n/a
 Args    : paths to input/output files

=cut

1;
