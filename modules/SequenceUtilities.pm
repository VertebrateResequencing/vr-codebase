package SequenceUtilities;

use strict;
use warnings;
use Carp;

=head2 read_fasta

  Arg [1]    : file of multiple sequences in FASTA 
  Example    : read_fasta('insertion_reads.fa');
  Description: Reads a FASTA file and returns a hash where the keys
               are the sequence names read
  Returntype : hash ref of sequences read, keyed on sequence names

=cut


sub read_fasta {
    my ($fasta_file) = shift;

    my (%name2seq, $name, $seq);
    my $concatenated_seq;

    open (FASTA, "< $fasta_file") ||
	croak("Cannae open $fasta_file\n");

    while (<FASTA>) {
	$concatenated_seq .= $_;
    }
    close(FASTA);

    # guess line ending and go with that:
    $concatenated_seq =~ /(\r\n|\r|\n)/;
    my $line_end = $1;

    my @lines = split $line_end, $concatenated_seq;
    
    foreach my $line(@lines){

	if ($line =~ /^>/) {

	   $line =~ s/>//;
	    if (defined $name) {
		$name2seq{$name} = $seq;
	    }
	    $name = $line;

	    if ($name =~ /^__/){
		# keys starting __ are reserved, so prepend a '!'
		$name = "!$name";
	    }

	    if ($name =~ /\s/) {
		my @arr = split(" ",$line);
		$name = $arr[0];
	    }
	    $seq = "";
	}
	else{
	    $line =~ s/[^ATGCNatgcn]//g;  # clean up sequence
	    $seq .= uc($line);	# make sequence uppercase
	}
    }
    # store the last read sequence
    $name2seq{$name} = uc($seq);
    
    return (\%name2seq);
}



=head2 read_tag_file

  Arg [1]    : file of a single sequence in FASTA 
  Example    : read_tag_file('piggyBAC.fa');
  Description: Reads a FASTA file and the name and sequence of the 
               tag found
  Returntype : Strings of the tag name and its sequence

=cut

sub read_tag_file {
    my ($fasta_file) = shift;

    my ($tag_name, $tag_seq);
    my $concatenated_seq;

    open (FASTA, "< $fasta_file") ||
	croak("Cannae open $fasta_file\n");

    while (<FASTA>) {
	$concatenated_seq .= $_;
    }
    close(FASTA);

    # guess line ending and go with that:
    $concatenated_seq =~ /(\r\n|\r|\n)/;
    my $line_end = $1;

    my @lines = split $line_end, $concatenated_seq;

    my $tag_cnt = 0;
    foreach my $line(@lines){

	if ($line =~ /^>/) {
	    $tag_cnt++;

	    if ( $tag_cnt > 2 ) {
		croak("Multiple tags found in tag file when expecting one.\n");
	    }
	    $line =~ s/>//;
	    $tag_name = $line;
	}
	else{
	    $line =~ s/[^ATGCNatgcn]//g;  # clean up sequence
	    $tag_seq .= uc($line);	  # make sequence uppercase
	}
    }
    # store the last read sequence
    $tag_seq = uc($tag_seq);
    
    return ($tag_name,$tag_seq);
}



=head2 filter_by_sequence

  Arg [1]    : a string sequence to filter upon
  Arg [2]    : hash ref of sequence names and sequence
  Example    : filter_by_sequence('acgt',\%mySequences)
  Description: Filters a hash ref of sequences for the supplied sequence
               string
  Returntype : hash ref of filtered sequences

=cut

sub filter_by_sequence($$) {
    my ($raw_filter_seq,$seq_hash,) = @_;

    my (%filtered_seqs);

    # clean up the sequence to filter against
    my $filter_seq = uc($raw_filter_seq);

    # and sanity-check it that it is a DNA sequence
    if ( $filter_seq !~ /^[ACGTN]+$/ ) {
	croak "Non-DNA bases in filter sequence: $filter_seq\n";
    }

    # loop over the sequences and reject those containing
    # the supplied sequence / pattern

    foreach my $curr_seq ( keys %{$seq_hash} ) {
	
	# retain those sequences that DON'T have the filter sequence
	if ( $seq_hash->{$curr_seq} !~ /$filter_seq/ ) {
	    $filtered_seqs{$curr_seq} = $seq_hash->{$curr_seq};
	}
    }
    
    return (\%filtered_seqs);
}



=head2 match_by_sequence

  Arg [1]    : a string sequence to filter upon
  Arg [2]    : hash ref of sequence names and sequence
  Example    : match_by_sequence('acgt',\%mySequences)
  Description: Filters a hash ref of sequences for the supplied sequence
               string
  Returntype : hash ref of filtered sequences

=cut

sub match_by_sequence($$) {
    my ($raw_filter_seq,$seq_hash,) = @_;

    my (%filtered_seqs);

    # clean up the sequence to filter against
    my $filter_seq = uc($raw_filter_seq);

    # and sanity-check it that it is a DNA sequence
    if ( $filter_seq !~ /^[ACGTN]+$/ ) {
	croak "Non-DNA bases in filter sequence: $filter_seq\n";
    }

    # loop over the sequences and reject those containing
    # the supplied sequence / pattern

    foreach my $curr_seq ( keys %{$seq_hash} ) {
	
	# retain those sequences that DO have the filter sequence
	if ( $seq_hash->{$curr_seq} =~ /$filter_seq/ ) {
	    $filtered_seqs{$curr_seq} = $seq_hash->{$curr_seq};
	}
    }
    
    return (\%filtered_seqs);
}



=head2 clip_sequence

  Arg [1]    : a string sequence to filter upon
  Arg [2]    : hash ref of sequence names and sequence
  Example    : match_by_sequence('acgt',\%mySequences)
  Description: Filters a hash ref of sequences for the supplied sequence
               string
  Returntype : hash ref of filtered sequences

=cut

sub clip_sequence($$$) {
    my ($seq,$clip_point,$is_complement) = @_;

    my $clipped_seq;
    
    if ( $is_complement == 1 ) {

	#              Clip point
	#          <---*
	# AGATGATAGATGATAGATGAGAGCGCGCGC
	#               -------------
	#               Match

	$clipped_seq = substr($seq,0,$clip_point-1);
    }
    else {

	#              Clip point
	#              *--->
	# AGATGATAGATGATAGATGAGAGCGCGCGC
	# -------------
	# Match

	$clipped_seq = substr($seq,$clip_point-1);
    }
    
    return ($clipped_seq);
}



=head2 slice_sequence

  Arg [1]    : a string sequence to filter upon
  Arg [2]    : hash ref of sequence names and sequence
  Example    : match_by_sequence('acgt',\%mySequences)
  Description: Filters a hash ref of sequences for the supplied sequence
               string
  Returntype : hash ref of filtered sequences

=cut

sub slice_sequence($$$) {
    my ($seq,$start,$end) = @_;

    if ( $end < $start ) {
      die("Trying to clip a sequence where the end is less than the start\n");
    }

    my $clipped_seq;

    my $slice_length = $end - $start + 1;
    $clipped_seq = substr($seq,$start-1,$slice_length);

    return ($clipped_seq);
}



=head2 rev_comp

  Arg [1]    : string sequence
  Example    : rev_comp($mySeq)
  Description: Reverse complements a DNA sequence
  Returntype : string

=cut

sub rev_comp {
    my ($seq) = @_;
    my $revcomp=reverse($seq);
    $revcomp=~ tr/atgcATGC/tacgTACG/;
    return $revcomp;
}



sub find_exact_matches{
    my ($seq, $seq2match) = @_;

    my @matches;

    while ($seq =~ /$seq2match/g){ 
	my $start_pos = pos($seq) - length($seq2match) + 1;
	my $end_pos = pos($seq);

	my $rec = {};
	$rec->{start} = $start_pos;
	$rec->{end} = $end_pos;

        push @matches, $rec;
    }
    return \@matches;
}




=head2 read_barcodes

  Arg [1]    : file barcodes
  Example    : read_barcodes($my_barcodes);
  Description: Reads a text file of barcodes, where each barcode
               is on a separate line in the file
  Returntype : hash ref of barcodes read, keyed on generated barcode
               names (BC..)

=cut

sub read_barcodes {
    my ($barcodes_file) = shift;


    # some simple error checking
    if ( !defined $barcodes_file ) {
      die ("No barcodes file specified\n");
    }
    if ( !-e $barcodes_file ) {
      die ("File of barcodes does not appear to exist\n");
    }


    my (%barcodes);

    my $barcode_cnt = 0;
    open (BAR, "< $barcodes_file") ||
	croak("Cannae open $barcodes_file\n");

    while (<BAR>) {
      if ( /^\n/ || /^\#/) {
	next;
      }

      chomp;
      $barcode_cnt++;

      my $formatted_bc = sprintf("%03d",$barcode_cnt);
      my $barcode_key = 'BC' . $formatted_bc;

      $barcodes{$barcode_key} = $_;
    }
    close(BAR);

    return (\%barcodes);
}


1;
