package Pathogens::Variant::Sample;
use Moose;


has 'name' => (is => 'rw', isa => 'Str', required => 1);
has 'alleles' => (is => 'rw', isa => sub {[]}) ;
has 'homozygous' => (is => 'rw', isa => 'Bool');
has 'genotype' => (is => 'rw', isa => 'Str');
 
=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 USAGE

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

    Feyruz Yalcin
    CPAN ID: FYALCIN

=head1 SEE ALSO

=cut

1;
