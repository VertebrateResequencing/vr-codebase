#!/usr/local/bin/perl

=head1 NAME

DASGFF.pm - module to write GFF of DAS data suitable for pro_load.pl

=head1 DESCRIPTION

Provides writeGFF which takes a reference to an array of hashrefs.  Writes the data from the hashes to STDOUT in GFF format suitable for reading by pro_load.pl to load into a ProServer database.

=head1 AUTHOR

Jim Stalker (jws@sanger.ac.uk)

=cut

package ProServer::DASGFF;

use strict;
use vars qw(@ISA @EXPORT);
require  Exporter;
@ISA     = qw(Exporter);
@EXPORT  = qw(writeGFF);


$|++;

sub writeGFF {
    my $featref = shift;

    foreach my $feature (@$featref){
	# Get & check the essential GFF fields
	my $segment = $feature->{'segment'};
	die "No segment defined\n" unless defined $segment;

	my $start   = $feature->{'start'};
	die "Start ($start) not a number\n" unless $start =~ /^\d+$/;

	my $end	    = $feature->{'end'};
	die "End ($end) not a number\n" unless $end =~ /^\d+$/;
	
	($start, $end) = ($end, $start) unless $start <= $end;

	my $method = $feature->{'method'};

	my $type = $feature->{'type'} || $feature->{'type_id'};
	die "feature does not contain a type or type_id field\n" unless $type;

	my $score   = $feature->{'score'};
	$score = "." unless defined $score;

	my $orient  = $feature->{'orient'};
	$orient = "+" if $orient eq "+1";
	$orient = "-" if $orient eq "-1";
	$orient = "." unless $orient =~ /^[+-]$/;

	my $phase   = $feature->{'phase'};
	$phase = "." unless defined $phase and $phase =~ /^[012]$/;

	$feature->{'type_category'} ||= 'similarity';

	my @other_das_fields = qw(id label type_id type_category group_id
	group_label group_type group_link_url group_link_text target_id
	target_start target_end link_url link_text note group_note);
	
	my @attribs;
	foreach (@other_das_fields){
	    my $value = $feature->{$_};
	    next unless defined $value;
	    $value =~ s/\t/\\t/g;
	    $value =~ s/\r/\\r/g;
	    $value =~ s/\n/\\n/g;
	    push @attribs, qq($_ "$value");
	}

	# GFF2 spec specifies semicolon separators, not space-semicolon-space,
	# but the examples use it, and so apparently does everyone else.
	# Makes parsing easier, anyway...
	my $attribs = join (" ; ", @attribs);

	print STDOUT "$segment\t$method\t$type\t$start\t$end\t$score\t$orient\t$phase\t$attribs\n";

    }
}

1;

