package NeedleWrapper;

=pod

=head1 NAME

 NeedleWrapper.pm


=head1 DESCRIPTION

 This module is a wrapper for the EMBOSS needle sequence
 global alignment tool.


=head1 ALGORITHM HOMEPAGE

 http://emboss.sourceforge.net/apps/release/6.0/emboss/apps/needle.html


=head1 SYNOPSIS

 use NeedleWrapper;

 # creates a new NeedleWrapper object with default 
 # parameter settings
 my $cmw = NeedleWrapper->new();

 # runs the algorithm and returns the results as an
 # array of hashes
 # assumes that the sequences are presented in the order
 # of the trace first (longest, probably, sequence) and
 # the tag sequence (shorter) second

 my @alignments = @{$cmw->do_pairwise_alignment($trace_seq_name,
                                                $trace_seq,
                                                $tag_seq_name,
                                                $tag_seq};

 # e.g.
 #  my @alignments = @{$cmw->do_pairwise_alignment(
 #                                 'PLK_1',
 #                                 'AGCTGCTGCTGCTGCTGCTGTCC',
 #                                 'PiggyBAC',
 #                                 'TATCTTTCTAGGGTTAA')};

=head1 AUTHOR

 Alistair G. Rust, <ar12@sanger.ac.uk>
 Modified from Ssaha_wrapper.pm written by sr7

=cut


BEGIN {
  unshift(@INC, "/software/pubseq/PerlModules/BioPerl/1_5_2");
}


use strict;
use warnings;
use Carp;
use Bio::AlignIO;

my $MODULE_NAME = __PACKAGE__;


=head1 METHODS

=head2 new

  Example    : my $cmw->new()
  Description: Creates a new instance of a NeedleWrapper object
  Args       : None
  Returns    : NeedleWrapper object

=cut

sub new {
    my ($class,@args)= @_;

    my $self= {};

    # set up some default parameters
    $self->{'params'}={
	'gapopen' => 10,
	'gapext'  => 1,
    };

    $self->{'BIN_DIR'} =
	'/software/pubseq/bin/EMBOSS-5.0.0/bin/';

    bless($self, $class);


    # if an argument string was provided, override the default parameters where
    # applicable
    if ( @args ) {
      $self->_processCallerArguments(@args);
    }

    return $self;
}



=head2 do_pairwise_alignment

  Example    : my @alignments = 
                 @{$cmw->do_pairwise_alignment($seq_a_name,
					       $seq_a,
					       $seq_b_name,
					       $seq_b)};
  Description: Performs a pairwise cross_match comparison on
               two sequences
  Arg[1]     : Name of sequence A
  Arg[2]     : Sequence A
  Arg[3]     : Name of sequence B
  Arg[4]     : Sequence B
  Returns    : On success returns a hash of results
               E.g.
                %needle_results = {
                     score => 14,
		     start => 100,
		     end   => 117}
                }

               On failure to find any alignments returns undef

=cut

sub do_pairwise_alignment{

    my ($self,$seq_a_name,$seq_a,$seq_b_name,$seq_b)= @_;

    # create the individual sequence files required by cross_match
    my $seq_a_file = $self->_create_fasta_file($seq_a_name,$seq_a);
    my $seq_b_file = $self->_create_fasta_file($seq_b_name,$seq_b);

    # create a temporary output file
    my $tmp_op_file = "/tmp/" . int(rand(1000000)) . 
	"-" . int(rand(1000000)) . "-needle-output.txt";


    # check the global alignment on the forward strand
    # returned alignment is a Bio::SimpleAlign object
    my $alignment = $self->_do_pairwise_alignment($seq_a_file,
						  $seq_b_file,
						  $tmp_op_file);

    # process the alignment
    my $align_start;
    my $align_end;

    # assuming that sequence b is the shortest sequence
    # run along the tag and see where the bases match
    my $tag_length = length($seq_b);

    for my $curr_pos ( 1...$tag_length ) {
      my $pos = $alignment->column_from_residue_number($seq_b_name, $curr_pos);

      if ( $curr_pos == 1 ) {
	$align_start = $pos;
      }
      else {
	$align_end = $pos;
      }
    }

    my %alignment;
    $alignment{start} = $align_start;
    $alignment{end} = $align_end;
    $alignment{strand} = '+';
    $alignment{score} = $alignment->score();
    $alignment{identity} = $alignment->percentage_identity();

    # clean up
    unlink($seq_a_file);
    unlink($seq_b_file);
    unlink($tmp_op_file);

    return ( %alignment );
}



=head2 do_multistrand_pairwise_alignment

  Example    : my @alignments = 
                 @{$cmw->do_pairwise_alignment($seq_a_name,
					       $seq_a,
					       $seq_b_name,
					       $seq_b)};
  Description: Performs two pairwise cross_match comparisons on
               two sequences by considering sequence B both on
               the forward and reverse strands
  Arg[1]     : Name of sequence A
  Arg[2]     : Sequence A
  Arg[3]     : Name of sequence B
  Arg[4]     : Sequence B
  Returns    : On success returns a hash of results
               E.g.
                %needle_results = {
                     score => 14,
		     start => 100,
		     end   => 117}
                }

               On failure to find any alignments returns undef

=cut

sub do_multistrand_pairwise_alignment{

    my ($self,$seq_a_name,$seq_a,$seq_b_name,$seq_b)= @_;

    my $revcomp_seq_b = $self->rev_comp($seq_b);

    # create the individual sequence files required by cross_match
    my $seq_a_file = $self->_create_fasta_file($seq_a_name,$seq_a);
    my $seq_b_file = $self->_create_fasta_file($seq_b_name,$seq_b);
    my $revcomp_seq_b_file = $self->_create_fasta_file($seq_b_name,$revcomp_seq_b);

    # create a temporary output file
    my $tmp_op_file = "/tmp/" . int(rand(1000000)) . 
	"-" . int(rand(1000000)) . "-needle-output.txt";


    # check the global alignment on the forward strand
    # returned alignment is a Bio::SimpleAlign object
    my $fwd_alignment = $self->_do_pairwise_alignment($seq_a_file,
						      $seq_b_file,
						      $tmp_op_file);

    # now check with a reverse complemented version of sequence b
    my $revcomp_alignment = $self->_do_pairwise_alignment($seq_a_file,
							  $revcomp_seq_b_file,
							  $tmp_op_file);


    # determine the highest scoring alignment
    my $best_alignment;
    my $strand = '+';

    if ( $fwd_alignment->score() == $revcomp_alignment->score()) {
      carp("Needle alignments both have similar scores => " . 
	   $fwd_alignment->score() . "\n");

      # take the forward strand match
      $best_alignment = $fwd_alignment;
    }
    elsif( $fwd_alignment->score() > $revcomp_alignment->score() ) {
      $best_alignment = $fwd_alignment;
    } # $fwd_alignment->score() < $revcomp_alignment->score()
    else {
      $best_alignment = $revcomp_alignment;
      $strand = '-';
    }


    # process the best alignment
    my $align_start;
    my $align_end;

    # assuming that sequence b is the shortest sequence
    # run along the tag and see where the bases match
    my $tag_length = length($seq_b);

    for my $curr_pos ( 1...$tag_length ) {
      my $pos = $best_alignment->column_from_residue_number($seq_b_name, $curr_pos);

      if ( $curr_pos == 1 ) {
	$align_start = $pos;
      }
      else {
	$align_end = $pos;
      }
    }

    # results hash
    my %alignment;
    $alignment{start} = $align_start;
    $alignment{end} = $align_end;
    $alignment{strand} = $strand;
    $alignment{score} = $best_alignment->score();
    $alignment{identity} = $best_alignment->percentage_identity();


    # clean up
    unlink($seq_a_file);
    unlink($seq_b_file);
    unlink($tmp_op_file);
    unlink($revcomp_seq_b_file);

    return ( %alignment );
}



sub _do_pairwise_alignment{

    my ($self,$fasta_seq_a,$fasta_seq_b,$op_file)= @_;

    # build up the cross_match command
    my $sys_call= 
	$self->{'BIN_DIR'} . "needle " . $fasta_seq_a . " " 
	  . $fasta_seq_b . " -outfile " . $op_file;

    foreach my $option(sort keys %{$self->{params}} ){
        $sys_call .= " -" . $option . " " . $self->{params}->{$option};
    }

    # redirect STDERR to suppress verbose messages
    $sys_call .= " 2>/dev/null";

    eval{
	system($sys_call);
    };

    if($@){
        #error occurred
	warn "EMBOSS needle: $@";
        return undef;
    }

    my $op_alignment = Bio::AlignIO->new(-file   => $op_file,
					 -format => 'emboss');

    my $num_alignments = 0;

    my $aln;
    while( my $curr_aln = $op_alignment->next_aln ) {
      $aln = $curr_aln;
      $num_alignments++;
    }

    if ( $num_alignments > 1 ) {
      carp "More than one alignment found for a global pairwise alignment\n";
    }

    return ($aln);
}



##########################
#
# GETTERS/SETTERS METHODS
#
##########################

=head2 gap_ext

  Example    : $cmw->gap_ext(-1)
  Description: A get/set method for setting the gap_ext value
  Arg        : int less than zero
  Returns    : int

=cut

sub gapext {

    my $self = shift;
    my $curr_gapext = shift;

    if ( $curr_gapext ) {
	if ( $curr_gapext >= 0 ) {
	    croak "gapext term should be negative\n";
	}
	else {
	    $self->{'params'}{'gapext'} = $curr_gapext;
	}
    }
    return($self->{'params'}{'gapext'});
}



=head2 gap_init

  Example    : $cmw->gap_init(-1)
  Description: A get/set method for setting the gap_init value
  Arg        : int less than zero
  Returns    : int

=cut

sub gap_init {

    my $self = shift;
    my $curr_gap_init = shift;

    if ( $curr_gap_init ) {
	if ( $curr_gap_init >= 0 ) {
	    croak "gap_init term should be negative\n";
	}
	else {
	    $self->{'params'}{'gap_init'} = $curr_gap_init;
	}
    }
    return($self->{'params'}{'gap_init'});
}



=head2 minscore

  Example    : $cmw->minscore(20)
  Description: A get/set method for setting the minscore value
  Arg        : int greater than 0
  Returns    : int

=cut

sub minscore {

    my $self = shift;
    my $curr_minscore = shift;

    if ( $curr_minscore ) {
	if ( $curr_minscore <= 0 ) {
	    croak "minscore term should be greater than 0\n";
	}
	else {
	    $self->{'params'}{'minscore'} = int($curr_minscore);
	}
    }
    return($self->{'params'}{'minscore'});
}



=head2 penalty

  Example    : $cmw->penalty(-2)
  Description: A get/set method for setting the penalty value
  Arg        : int less than zero
  Returns    : int

=cut

sub penalty {

    my $self = shift;
    my $curr_penalty = shift;

    if ( $curr_penalty ) {
	if ( $curr_penalty >= 0 ) {
	    croak "Penalty term should be negative\n";
	}
	else {
	    $self->{'params'}{'penalty'} = int($curr_penalty);
	}
    }
    return($self->{'params'}{'penalty'});
}



##################
#
# PRIVATE METHODS
#
##################

=head1 PRIVATE METHODS

=head2 _create_fasta_file

  Example    : $self->_create_fasta_file('mySeq','acctgctgctcgt')
  Description: A 'private' method that creates a temporary FASTA file
  Arg[1]     : Sequence name
  Arg[2]     : DNA sequence (string)
  Returns    : Name of temporary file

=cut

sub _create_fasta_file {

    my ($self, $seq_name, $seq) = @_;

    my $tmp_fasta_file = "/tmp/" . int(rand(1000000)) .
	"-" . int(rand(1000000)) . ".fa";

    open (FASTA, "> $tmp_fasta_file") ||
	croak("Cannot open $tmp_fasta_file\n");

    print FASTA ">$seq_name\n";
    print FASTA "$seq\n";

    close(FASTA);

    return ($tmp_fasta_file);
}



=head2 rev_comp

  Arg [1]    : string sequence
  Example    : rev_comp($mySeq)
  Description: Reverse complements a DNA sequence
  Returntype : string

=cut


sub rev_comp {
    my ($self,$seq) = @_;
    my $revcomp=reverse($seq);
    $revcomp=~ tr/atgcATGC/tacgTACG/;
    return $revcomp;
}


=head2 _processCallerArguments

  Example    : $self->_processCallerArguments(@pars)
  Description: A 'private' method that processes the options passed
               when a new object is created
  Args       : array of parameters in the form -param value
  Returns    : nothing

=cut

sub _processCallerArguments {

    my ($self,@params) = @_;

    while (@params) {
	my $switch = shift(@params);
	my $val = shift(@params);

	if ( $switch !~ /^-/ ) {
	    die("Unknown option when processing caller arguments: $switch");
	}

	# chop off the - from the switch
	$switch =~ s/^-//;

	# overwrite the default value
	if ( exists $self->{$switch} ) {
	  $self->{$switch} = $val;
	}
	else {
	  die("$MODULE_NAME Trying to set an unknown command line option: $switch");
	}
    }
}

1;
