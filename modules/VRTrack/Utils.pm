package VRTrack::Utils;
=head1 NAME

VRTrack::Utils

=head1 SYNOPSIS
    my $enumvals = VRTrack::Utils::list_enum_vals($dbh,$table,$column);

=head1 DESCRIPTION

Some utility methods for VRTrack to share

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;
no warnings 'uninitialized';

=head2 list_enum_vals

  Arg [1]    : database handle
  Arg [2]    : table name
  Arg [3]    : column name
  Example    : my $vals = VRTrack::Utils::list_enum_vals($dbh,'library','qc_status');
  Description: retrieves the list of allowed enum values for a column in lowercase.  Dies if the column is not of type enum
  Return_type_id : array ref

=cut

sub list_enum_vals {
    my ($dbh, $table, $column) = @_;
    my $row = $dbh->selectrow_hashref("SHOW COLUMNS FROM $table LIKE ?", undef,$column);
    my $type = lc($row->{Type});
    unless ($type =~ /^enum/){
        die "$table:$column is not of type enum";
    }
    $type =~ s/^enum\('//;
    $type =~ s/'\)$//;
    my @vals = split /','/, $type;
    return \@vals;
}

1;
