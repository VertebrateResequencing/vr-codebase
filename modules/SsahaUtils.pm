package SsahaUtils;

use strict;

# defines
my $SSAHA_SCORE_THRESHOLD = 20;
my $SSAHA_FOLD_DIFFERENCE = 1.5;

my $DEBUG = 0;


sub parseForMouseInsertionSites {

  my $ssaha_hits_hash_ref = shift;

  my %results_hash;

  foreach my $curr_seq ( keys %$ssaha_hits_hash_ref ) {

    $results_hash{$curr_seq}{status} = 'NA';
    $results_hash{$curr_seq}{chr} = 'NA';
    $results_hash{$curr_seq}{genome_start} = 'NA';
    $results_hash{$curr_seq}{genome_end} = 'NA';
    $results_hash{$curr_seq}{alignment_score} = 'NA';
    $results_hash{$curr_seq}{mapping_score} = 'NA';
    $results_hash{$curr_seq}{strand} = 'NA';

    # if ssaha returns no alignments for a trace
    if ( !defined $$ssaha_hits_hash_ref{$curr_seq} ) {
      $results_hash{$curr_seq}{status} = 'no_genome_alignment';
      next;
    }


    # check that the best alignment starts from the very first nt
    # of the trace.  if not it implies that the insertion is not
    # actually at the real biological site (TA in the current case)

    if ( $$ssaha_hits_hash_ref{$curr_seq}[0]->{'direction'} eq 'F' ) {
      if ( $$ssaha_hits_hash_ref{$curr_seq}[0]->{'q_start'} != 1 ) {
	$results_hash{$curr_seq}{status} = 'distant_from_insertion_site';
	next;
      }
    }
    elsif ( $$ssaha_hits_hash_ref{$curr_seq}[0]->{'direction'} eq 'C' ) {
      if ( $$ssaha_hits_hash_ref{$curr_seq}[0]->{'q_end'} != 1 ) {
	$results_hash{$curr_seq}{status} = 'distant_from_insertion_site';
	next;
      }
    }


    # if the best alignment is above the score threshold
    if ($$ssaha_hits_hash_ref{$curr_seq}[0]->{score} >= $SSAHA_SCORE_THRESHOLD){

      my $best_hit = $$ssaha_hits_hash_ref{$curr_seq}[0];
      my $chr = '';
      my $mapping_score = -1;

      # very clunky ...

      ################################
      # check for multiple alignments
      ################################
      if(defined $$ssaha_hits_hash_ref{$curr_seq}[1]) {

	#print "Number of alignments => " . scalar @{$$ssaha_hits_hash_ref{$curr_seq}} . "\n";
	my $num_of_alignments = scalar @{$$ssaha_hits_hash_ref{$curr_seq}};

	#print $$ssaha_hits_hash_ref{$curr_seq}[0]->{score}  . 
	#" " .  $$ssaha_hits_hash_ref{$curr_seq}[1]->{score} . " ";
	#print $num_of_alignments . "\n";

	$mapping_score = 
	  ( $$ssaha_hits_hash_ref{$curr_seq}[0]->{score} - $$ssaha_hits_hash_ref{$curr_seq}[1]->{score} )
	    - ( 2 * scalar @{$$ssaha_hits_hash_ref{$curr_seq}} );


	# loop through all the reported alignments to find one that also
	# maps from the first nt of the trace (i.e. should ensure that the
	# neighbouring nts are a TA site).

	# heavily assumes that the ssaha hits are returned sorted on score
	my $best_matching_hit = 0;

	for ( my $i = 1; $i < $num_of_alignments; $i++ ) {
	  if ( $$ssaha_hits_hash_ref{$curr_seq}[$i]->{'direction'} eq 'F' &&
	       $$ssaha_hits_hash_ref{$curr_seq}[$i]->{'q_start'} == 1 ) {
	    $best_matching_hit = $i;
	    last;
	  }
	  elsif ( $$ssaha_hits_hash_ref{$curr_seq}[$i]->{'direction'} eq 'C' &&
		  $$ssaha_hits_hash_ref{$curr_seq}[$i]->{'q_end'} == 1 ) {
	    $best_matching_hit = $i;
	    last;
	  }
	}


	# if another alignment begins from nt 1 then check how good that hit is
	if ( $best_matching_hit != 0 ) {
	  # if the second best score is less than the threshold score
	  # flag the best score as being unique
	  if ( $$ssaha_hits_hash_ref{$curr_seq}[$best_matching_hit]->{score} < $SSAHA_SCORE_THRESHOLD ) {
	    $chr = $best_hit->{'s_name'};
	    $results_hash{$curr_seq}{status} = 'unique_alignment';
	  }

	  # if the difference between the best and second best hit is
	  # above some threshold, again flag the hit as being unique
	  elsif ( $$ssaha_hits_hash_ref{$curr_seq}[0]->{score} >=
		  $SSAHA_FOLD_DIFFERENCE * $$ssaha_hits_hash_ref{$curr_seq}[$best_matching_hit]->{score}) {

	    $chr = $best_hit->{'s_name'};
	    $results_hash{$curr_seq}{status} = 'unique_alignment';
	  }

	  # if the second best hit, is 'close' then flag the trace as
	  # having multiple alignments
	  else {
	    $chr = "MULTI";
	    $results_hash{$curr_seq}{status} = 'multiple_alignments';
	  }
	}
	# otherwise flag the best alignment as being unique
	else {
	    $chr = $best_hit->{'s_name'};
	    $results_hash{$curr_seq}{status} = 'unique_alignment';
	}
      }

      ###############
      # a unique hit
      ###############
      else {
	$chr = $best_hit->{'s_name'};
	$results_hash{$curr_seq}{status} = 'unique_alignment';
	$mapping_score = 255;
      }

      # store the alignment result into the results hash
      $results_hash{$curr_seq}{chr} = $chr;
      $results_hash{$curr_seq}{genome_start} = $best_hit->{'s_start'};
      $results_hash{$curr_seq}{genome_end} = $best_hit->{'s_end'};
      $results_hash{$curr_seq}{alignment_score} = $best_hit->{'score'};
      $results_hash{$curr_seq}{mapping_score} = $mapping_score;

      # sort out the strand
      if ( $best_hit->{'direction'} eq 'F' ) {
	$results_hash{$curr_seq}{strand} = '+';
      }
      elsif ( $best_hit->{'direction'} eq 'C' ) {
	$results_hash{$curr_seq}{strand} = '-';
      }
    }
    else {
      $results_hash{$curr_seq}{status} = 'weak_alignment';
    }
  }

  return (\%results_hash);
}

1;
