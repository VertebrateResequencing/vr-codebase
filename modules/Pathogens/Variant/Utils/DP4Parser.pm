package DP4Parser;
use Moose;
use namespace::autoclean;


sub create_dp4_object {
    
    my ($self, $dp4string) = @_;

    #format of the $dp4string:
    #ref-forward bases, ref-reverse bases, alt-forward bases, alt-reverse bases
    my ($ref_forward, $ref_reverse, $alt_forward, $alt_reverse) = split(',', $dp4string);
        
    my (
        , $ratio_forward_reference_bases
        , $ratio_forward_alternative_bases
        , $ratio_reverse_reference_bases
        , $ratio_reverse_alternative_bases
    );
    
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
    return Pathogens::Variant::Util::DP4->new(
                                                forward_reference_strand_depth => $ref_forward
                                              , forward_alternative_strand_depth => $alt_forward
                                              , reverse_reference_strand_depth => $ref_reverse
                                              , reverse_alternative_strand_depth => $alt_reverse
                                              , 
                                             )
}

1;