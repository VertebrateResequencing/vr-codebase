package Pathogens::Variant::Utils::DP4Parser;
use Moose;
use namespace::autoclean;


sub parse {

    my ($self, $dp4string) = @_;

    my (
          $ratio_forward_reference_bases
        , $ratio_forward_alternative_bases
        , $ratio_reverse_reference_bases
        , $ratio_reverse_alternative_bases
        , $ref_total
        , $alt_total
    );

    #format of the $dp4string:
    #ref-forward bases, ref-reverse bases, alt-forward bases, alt-reverse bases
    my (   $ref_forward                                 #0
         , $ref_reverse                                 #1
         , $alt_forward                                 #2
         , $alt_reverse) = split(',', $dp4string);      #3
    
    $ref_total = $ref_forward + $ref_reverse;
    $alt_total = $alt_forward + $alt_reverse;

    if (($ref_forward + $alt_forward) > 0 ) {
        $ratio_forward_reference_bases = $ref_forward/($ref_forward + $alt_forward); # 1st div 1st + 3rd
        $ratio_forward_alternative_bases = $alt_forward/($ref_forward + $alt_forward); # 3rd div 1st + 3rd
    } else {
        $ratio_forward_reference_bases = $ratio_forward_alternative_bases = 0; #all zero
    }
    
    if ( ($ref_reverse+$alt_reverse) > 0 ) {
        $ratio_reverse_reference_bases = $ref_reverse/($ref_reverse + $alt_reverse); # 2nd div 2nd + 4th
        $ratio_reverse_alternative_bases = $alt_reverse/($ref_reverse + $alt_reverse); # 4th div 2nd + 4th
    } else {
        $ratio_reverse_reference_bases = $ratio_reverse_alternative_bases = 0;
    }
    
    my %parse_results = (
          count_reference_forward_bases   => $ref_forward
        , count_alternative_forward_bases => $alt_forward
        , count_reference_reverse_bases   => $ref_reverse
        , count_alternative_reverse_bases => $alt_reverse
        , count_alternative_bases         => $alt_total
        , count_referecence_bases         => $ref_total
        , ratio_forward_reference_bases   => $ratio_forward_reference_bases
        , ratio_forward_alternative_bases => $ratio_forward_alternative_bases
        , ratio_reverse_reference_bases   => $ratio_reverse_reference_bases
        , ratio_reverse_alternative_bases => $ratio_reverse_alternative_bases

    );

    return (\%parse_results);

}

__PACKAGE__->meta->make_immutable;
1;