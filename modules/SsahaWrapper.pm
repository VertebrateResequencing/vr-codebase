package SsahaWrapper;

=pod

=head1 NAME

 SsahaWrapper.pm


=head1 DESCRIPTION

 This is a wrapper module around the ssaha binary being used
 for the 1KGenomes project.


=head1 ALGORITHM HOMEPAGE

 http://www.sanger.ac.uk/Software/analysis/SSAHA2/


=head1 SYNOPSIS

 use SsahaWrapper;

 # creates a new SsahaWrapper object with default
 # parameter settings
 my $sw = SsahaWrapper->new();

 my %alignments = %{$sw->do_ssaha('species_hash_table_path',
                                  $my_fasta_file)};

=head1 AUTHOR

 Alistair G. Rust, <ar12@sanger.ac.uk>
 Modified from Ssaha_wrapper.pm written by sr7

=cut


use strict;
use warnings;
no warnings 'uninitialized';


# class variables
my $MODULE_NAME = __PACKAGE__;
my $SSAHA2_BIN;


BEGIN {
  # 1KGenomes ssaha executable
  $SSAHA2_BIN = '/lustre/sf4/1kgenomes/bin/ssaha2';

  if ( !-e $SSAHA2_BIN ) {
    die ("Unable to find defined Ssaha binary: $SSAHA2_BIN\n");
  }
}


=head1 METHODS

=head2 new

  Title      : new
  Example    : $fw->new(@pars)
  Description: Creates and initializes a new SsahaWrapper object
  Args       : [Optional] List of parameters as -param value
  Returns    : new SsahaWrapper object

=cut

sub new {
    my ($class,@args)= @_;

    my $self= {};

    # set up default Ssaha parameters
    $self->{'params'}={'-seeds'    => 2,
		       '-depth'    => 5,
		       '-score'    => 20,
		       '-identity' => 30.0,
		       '-align'    => 0,
		       '-output'   => 'ssaha2',
		    };

    bless($self, $class);

    # if an argument string was provided, override the default 
    # parameters where applicable
    if ( @args ) {
      $self->_processCallerArguments(@args);
    }

    return $self;
}


=head2 do_ssaha

  Title      : do_ssaha
  Example    : $sw->do_ssaha('species_hash_table_path',$my_fasta_file)
  Description: Runs ssaha for a given file of multiple sequences in
               FASTA format
  Arg[1]     : Path to pre-computed Ssaha hash tables
  Arg[2]     : FASTA file
  Returns    : hashref to hash of arrays or undef if no matches are
               identified
               %alignments = (
                  seq_one => (
                    [score     => $score,     q_name    => $q_name,
                     s_name    => $s_name,    q_start   => $q_start,
                     q_end     => $q_end,     s_start   => $s_start,
                     s_end     => $s_end,     direction => $direction,
                     num_bases => $num_bases, identity  => $identity]
                    [score => $score,         q_name    => $q_name, ...]
                   )
                  seq_two => (
                    [score => $score, q_name => $q_name, ...}]
                    []
                   )
                  ...
               )

=cut

sub do_ssaha{

  # input the species name to search against
  # and a array-ref of sequences in fasta format
  # on success returns a hashref to alignments as
  # a hash array
  #
  # %alignments = (
  #      seq_one => (
  #                  [{score => $score,         q_name => $q_name,
  #                    s_name => $s_name,       q_start => $q_start,
  #                    q_end  => $q_end,        s_start => $s_start,
  #                    s_end => $s_end,         direction => $direction,
  #                    num_bases => $num_bases, identity => $identity}]
  #                  [{score => $score,         q_name => $q_name, ...}]
  #                 )
  #      seq_two => (
  #                  [{score => $score, q_name => $q_name, ...}]
  #                  []
  #                 )
  #      ...
  #   )
  #
  # each sequence in the defined FASTA sequence file has an entry
  # in the results hash.  if no alignment is found the results
  # entry is 'undef'
  #
  # on failure returns undef

  my ($self, $species_hash_tables, $fasta_file)= @_;


  # some simple defensive programming
  if ( !defined $species_hash_tables ) {
    die("Species hash tables not defined\n");
  }

  my $head_file = $species_hash_tables . ".head";
  if ( !-e $head_file ) {
    die("Unable to find species hash tables body\n$head_file");
  }

  if ( !-e $fasta_file ) {
    die("Unable to find specified FASTA file\n$fasta_file");
  }



  # construct the ssaha command
  my $cmd= $SSAHA2_BIN;

  # tag on any specified params
  foreach my $option(keys %{$self->{params}} ){
    $cmd .= " " . $option . " " . $self->{params}->{$option};
  }

  $cmd .= " -save " . $species_hash_tables;
  $cmd .= " " . $fasta_file;

  #print $cmd . "\n";


  # hash to store alignment results
  my %ssaha_results;

  eval{
    open ( RES, "$cmd |");

    while (<RES>) {
      my $line = $_;
      chomp $line;

      my $q_name;

      if($line =~ /^Matches For Query .+: (\S+)/){
	$q_name= $1;
	$ssaha_results{$q_name}= undef;
      }

      if($line =~/^ALIGNMENT/) {
	# data line seems to start with ALIGNMENT
	#print "LINE: $line\n";

	my @fields= split(/\s+/, $line);

	my $curr_qname= $fields[2];

	# check that the query already has an entry
	# in the results hash, which should have
	# come from the previous check on
	# Matches For Query

	if ( !exists $ssaha_results{$curr_qname} ) {
	  die("Error processing Ssaha ouput\n" .
	     "Trying to store an unknown query\n");
	}

	my $href= {'score'     => $fields[1],
		   'q_name'    => $fields[2],
		   's_name'    => $fields[3],
		   'q_start'   => $fields[4],
		   'q_end'     => $fields[5],
		   's_start'   => $fields[6],
		   's_end'     => $fields[7],
		   'direction' => $fields[8],
		   'num_bases' => $fields[9],
		   'identity'  => $fields[10]};

	push @{$ssaha_results{$curr_qname}}, $href;
      }
    }
    close(RES);
  };

  if($@){
    #error occurred
    warn "SSAHA: $@";
    return undef;
  }
  else{
    return \%ssaha_results;
  }
}


#
# Getters / setters
#

=head2 ckmer

  Example    : $sWrapper->ckmer(10)
  Description: A get/set method for specifying the word size for 
               cross_match matching
  Arg        : int
  Returns    : int

=cut

sub ckmer {
  my $self = shift;
  if (@_) { 
    my $curr_ckmer = shift;
    if ( $curr_ckmer < 1 ) {
      die("$MODULE_NAME Trying to set an incorrect ckmer: $curr_ckmer");
    }
    $self->{'ckmer'} = int($curr_ckmer);
  }
  return $self->{'ckmer'};
}



=head2 cmatch

  Example    : $sWrapper->cmatch(10)
  Description: A get/set method for specifying the minimum match
               length for cross_match matching
  Arg        : int
  Returns    : int

=cut

sub cmatch {
  my $self = shift;
  if (@_) { 
    my $curr_cmatch = shift;
    if ( $curr_cmatch < 1 ) {
      die("$MODULE_NAME Trying to set an incorrect cmatch: $curr_cmatch");
    }
    $self->{'cmatch'} = int($curr_cmatch);
  }
  return $self->{'cmatch'};
}



=head2 depth

  Example    : $sWrapper->depth(42)
  Description: A get/set method for the number of hits to consider
               for alignment
  Arg        : int
  Returns    : int

=cut

sub depth {
  my $self = shift;
  if (@_) {
    my $curr_depth = shift;
    if ( $curr_depth < 1 ) {
      die("$MODULE_NAME Trying to set an incorrect depth: $curr_depth");
    }
    $self->{'depth'} = int($curr_depth);
  }
  return $self->{'depth'};
}



=head2 identity

  Example    : $sWrapper->identity(42)
  Description: A get/set method for specifying the minimum identity
               for a match to be reported
  Arg        : float
  Returns    : float

=cut

sub identity {
    my $self = shift;
    if (@_) { 
	my $curr_identity = shift;
	if ( $curr_identity < 0 || $curr_identity > 100) {
	  die("$MODULE_NAME Trying to set an incorrect identity: $curr_identity");
	}
	$self->{'identity'} = $curr_identity;
      }
    return $self->{'identity'};
}



=head2 seeds

  Example    : $sWrapper->seeds(10)
  Description: A get/set method for specifying the number of kmer
               matches required to flag a hit
  Arg        : int
  Returns    : int

=cut

sub seeds {
    my $self = shift;
    if (@_) { 
	my $curr_seeds = shift;
	if ( $curr_seeds < 0 ) {
	  die("$MODULE_NAME Trying to set an incorrect seeds: $curr_seeds");
	}
	$self->{'seeds'} = int($curr_seeds);
      }
    return $self->{'seeds'};
}



=head2 score

  Example    : $sWrapper->score(10)
  Description: A get/set method for specifying the minimum score for 
               a match to be reported
  Arg        : float
  Returns    : float

=cut

sub score {
    my $self = shift;
    if (@_) { 
	my $curr_score = shift;
	if ( $curr_score < 0 ) {
	  die("$MODULE_NAME Trying to set an incorrect score: $curr_score");
	}
	$self->{'score'} = $curr_score;
      }
    return $self->{'score'};
}




#
# Private Methods
#

=head2 _processCallerArguments

  Example    : $self->_processCallerArguments(@pars)
  Description: A 'private' method that processes the options passed
               when a new object is created
  Args       : array of parameters
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
	  die("Trying to set an unknown command line option: $switch");
	}
    }
}

1;
