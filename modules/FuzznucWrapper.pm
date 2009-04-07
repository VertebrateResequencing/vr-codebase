package FuzznucWrapper;

use vars qw(@ISA $VERSION);
use strict;

=pod

=head1 NAME

 FuzznucWrapper

=head1 SYNOPSIS

 This module is a wrapper for the EMBOSS fuzznuc tool that
 searches for a specified PROSITE-style pattern in nucleotide
 sequences.


=head1 DESCRIPTION

 my $fw = FuzznucWrapper->new();

 my @matches = @{$fw->do_fuzznuc($my_sequence_name,
                                 $my_dna_sequence,
                                 $search_pattern)};

 # sets the allowed number of mismatches against the pattern
 $fw->mismatches(1);

 # enables searching on the reverse complement of the
 # input sequence
 $fw->complement(1);


=head1 ALGORITHM HOMEPAGE

 http://emboss.sourceforge.net/apps/release/6.0/emboss/apps/fuzznuc.html

=cut


my $MODULE_NAME = __PACKAGE__;
$VERSION = '1.0';


=head1 METHODS

=head2 new

  Title      : new
  Example    : $fw->new(@pars)
  Description: Creates and initializes a new FuzznucWrapper object
  Args       : [Optional] List of parameters
  Returns    : new FuzznucWrapper object

=cut

sub new {
  my ($class,@args) = @_;

  my $self;

  # initialise the algorithm with its default parameters and settings
  $self->{'mismatch'} = 3;
  $self->{'complement'} = 0;

  $self->{'BIN_DIR'} =
    '/software/pubseq/bin/EMBOSS-4.1.0/bin/';

  bless($self,$class);

  # if an argument string was provided, override the default parameters where
  # applicable
  if ( @args ) {
      $self->_processCallerArguments(@args);
  }

  return $self;
}


=head2 do_fuzznuc

  Title      : do_fuzznuc
  Example    : my @matches =
                @{$self->do_fuzznuc('mySeq','agatgatgatgat','acgt')};
  Description: Runs fuzznuc to perform a search of a DNA sequence
  Arg[1]     : Sequence name
  Arg[2]     : DNA sequence
  Arg[3]     : DNA pattern to search for
  Returns    : Arrayref of matches

=cut

sub do_fuzznuc {

    my ($self,$seq_name,$seq,$search_pattern) = @_;

    my $tmp_seq_file = $self->_create_fasta_file($seq_name,$seq);

    # create a temporary output file
    my $tmp_op_file = "/tmp/" . int(rand(1000000)) . 
      "-" . int(rand(1000000)) . "-fuzznuc-output.txt";

    my $complement_call = "";
    if ( $self->{'complement'} == 1 ) {
	$complement_call = " -complement";
    }

    # fix to get fuzznuc to play nicely with bash
    $search_pattern =~ s/\(/\\(/g;
    $search_pattern =~ s/\)/\\)/g;

    my $fuzz_call = $self->{'BIN_DIR'} . "fuzznuc " .
      $tmp_seq_file . " -pattern " . $search_pattern .
	" -pmismatch " . $self->{mismatch} .
	  $complement_call .
	    " -outfile " . $tmp_op_file;

    # redirect STDERR to suppress verbose messages
    $fuzz_call .= " 2>/dev/null";

    #print STDERR $fuzz_call . "\n";


    system($fuzz_call) and
      die("Failure in fuzznuc: $fuzz_call");

      eval{
	system($fuzz_call);
    };

    if($@){
        #error occurred
	warn "EMBOSS fuzznuc: $@";
        return undef;
    }


    my @match_results;

    # only process if there are hits
    if ( -s $tmp_op_file ) {

      # now open the dat file and parse the results
      open ( RES, "< $tmp_op_file" ) ||
	die("Cannot open  $tmp_op_file");

      while ( <RES> ) {

	if ( /^\#/ || /^\s*$/ ) { 
	  next;
	}
	if ( /Start/ ) {
	  next;
	}

	s/^\s+//;  # chop off the initial white spaces
	chomp;

	my ( $match_start,
	     $match_end,
	     $pattern,
	     $mismatch,
	     $subseq ) =  split(/\s+/,$_);

	my $is_complement_match = 0;

	if ( $match_start > $match_end ) {
	  my $tmp = $match_start;
	  $match_start = $match_end;
	  $match_end = $tmp;
	  $is_complement_match = 1;
	}

	my $match = {};
	$match->{start} = $match_start;
	$match->{end} = $match_end;
	$match->{complement} = $is_complement_match;
	$match->{mismatch} = $mismatch;
	$match->{sequence} = $subseq;

	push @match_results, $match;
      }
      close RES;
    }

    unlink($tmp_op_file) ||
      die ("Unable to remove the temporary output file: $tmp_op_file\n");
    unlink($tmp_seq_file) ||
      die ("Unable to remove the temporary output file: $tmp_op_file\n");



    return ( \@match_results );
}



#
# Algorithm-specific parameters
#

=head1 MODULE-SPECIFIC PARAMETERS

=head2 mismatch

  Title      : mismatch
  Example    : $obj->mismatch(1)
  Description: A get/set method for setting the mismatch
  Arg        : int
  Returns    : int

=cut

sub mismatch {
    my $self = shift;
    if (@_) { 
	my $curr_mismatch = shift;
	if ( $curr_mismatch < 0 ) {
	    die("$MODULE_NAME Probably trying to set an incorrect mismatch: $curr_mismatch");
	}
	$self->{mismatch} = $curr_mismatch;
    }
    return $self->{mismatch};
}



=head2 complement

  Title      : complement
  Example    : $obj->complement(1)
  Description: A get/set method for setting the complement
  Arg        : int
  Returns    : int

=cut

sub complement {
    my $self = shift;
    if (@_) {
	my $curr_complement = shift;
	if ( $curr_complement != 1 && $curr_complement != 0 ) {
	    die("$MODULE_NAME Probably trying to set an incorrect complement: $curr_complement\n" .
		"Should be 0 or 1");
	}
	$self->{complement} = $curr_complement;
    }
    return $self->{complement};
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
