package Pathogens::Variant::Evaluator::Pseudosequence;
use Moose;
extends 'Pathogens::Variant::Root';
with 'Pathogens::Variant::Role::Evaluator';#forces this class to implement "evaluate" subroutine

use Pathogens::Variant::Utils::DP4Parser;
use Pathogens::Variant::EvaluationReporter;
use Pathogens::Variant::Utils::EventManipulator;

use Log::Log4perl qw(get_logger);
use Data::Dumper;

use namespace::autoclean;





has 'minimum_depth'          => ( is => 'ro', isa => 'Int', default => 4 );
has 'minimum_depth_strand'   => ( is => 'ro', isa => 'Int', default => 2 );
has 'minimum_ratio'          => ( is => 'ro', isa => 'Num', default => 0.8 );
has 'minimum_quality'        => ( is => 'ro', isa => 'Int', default => 0 );
has 'minimum_map_quality'    => ( is => 'ro', isa => 'Int', default => 1 );
has 'minimum_af1'            => ( is => 'ro', isa => 'Num', default => 0.95 );
has 'af1_complement'         => (is => 'ro',  isa => 'Num', default => 0.05);
has 'minimum_ci95'           => ( is => 'ro', isa => 'Num', default => 0.0 );
has 'minimum_strand_bias'    => ( is => 'ro', isa => 'Num', default => 0.001 );
has 'minimum_map_bias'       => ( is => 'ro', isa => 'Num', default => 0.001 );
has 'minimum_tail_bias'      => ( is => 'ro', isa => 'Num', default => 0.001 );
has 'minimum_base_quality_bias' => ( is => 'ro', isa => 'Num', default => 0.0 );
has 'reporter'            => ( is => 'ro', isa => 'Pathogens::Variant::EvaluationReporter', lazy => 1, default => sub { return Pathogens::Variant::EvaluationReporter->new } );

has '_event'               => ( is => 'rw', isa => 'Pathogens::Variant::Event');

has '_dp4_parser'          => ( is => 'ro', isa => 'Pathogens::Variant::Utils::DP4Parser', lazy => 1, default => sub { return Pathogens::Variant::Utils::DP4Parser->new } );
has '_event_manipulator'   => ( is => 'ro', isa => 'Pathogens::Variant::Utils::EventManipulator', lazy => 1, default => sub { return Pathogens::Variant::Utils::EventManipulator->new } );


sub evaluate {
    my ($self, $event) = @_;

    #logger is for DEBUGging purposes, little effect on overall performance
    my $logger = get_logger("Pathogens::Variant::Evaluator::Pseudosequence");
    
    #marks non polymorphic events
    $self->_mark_the_event_if_not_polymorphic($event);
    
    #Set the _event for this evaluation round
    $self->_event($event);
    
    #increments total number of event evaluations counter
    $self->reporter->inc_total_number_of_event_evaluations;
    
    #sub-evaluation of the values in the VCF's "INFO" field
    my $good_values_in_info_field = $self->_passed_vcf_info_field_evaluation;
    
    #sub-evaluation of indels(we will discard these)
    my $is_not_an_indel = $self->_is_not_an_indel;
    
    #sub-evaluation of heterozygous variants (we will remove secondary alleles later on)
    my $is_heterozygous =0;
    if ($self->_event->polymorphic and $self->_has_secondary_heterozygous_alternative_alleles) {
        $is_heterozygous = 1;
    }
    
    #sub-evaluation of VCF's "QUAL" field
    my $has_good_quality = $self->_passed_quality;


    #Bringing all sub-evaluation results together, and taking here the final decision on this variant
    if ( $good_values_in_info_field and $is_not_an_indel and $has_good_quality ) {

        #make it homozygous
        if ( $is_heterozygous ) {
            $self->_event_manipulator->remove_secondary_alternative_heterozygous_alleles($event);
        }

        $self->_event->passed_evaluation(1);
        
        
        if (not $self->_event->polymorphic) {
            $self->reporter->inc_counter_accepted_reference_calls; #increments the counter called "accepted_reference_calls" by 1
        } else {
            $self->reporter->inc_counter_accepted_snp_calls        #increments the counter called "accepted_snp_calls" by 1
        }

        $logger->is_debug() && $logger->debug("Event dump after passing the evaluation:...\n". Dumper($event) . "\nReporter dump after passing the evaluation:...\n". Dumper($self->reporter) );

        return 1; #PASSED OVERALL EVALUATION

    } else {

        $self->_event->passed_evaluation(0); #i.e. the event failed to pass the filters

        $logger->is_debug() && $logger->debug("Event dump after failing the evaluation:...\n". Dumper($event) . "\nReporter dump after failing the evaluation:...\n". Dumper($self->reporter) );

        return 0; #FAILED OVERALL EVALUATION 

    }
}

sub _mark_the_event_if_not_polymorphic {
    my ($self, $event) = @_;

    if ($event->alternative_allele eq '.') {
        $event->polymorphic(0); #set this one to false (0). Note default value is true
    }

}
sub _has_secondary_heterozygous_alternative_alleles {
    my ($self) = @_;

    my @num_alleles = split(',', $self->_event->alternative_allele);
    
    #seeing more than 1 element after splitting on comma, implies heterozygousity
    if (scalar @num_alleles > 1) {
        $self->reporter->inc_counter_heterozygous_calls;
        return 1;
    } else {
        return 0;
    }
}

sub _passed_quality {
    my ($self) = @_;

    if ($self->_event->quality < $self->minimum_quality) {
        $self->reporter->inc_counter_failed_quality;
        return 0;
    } else {
        return 1;
    }
}

sub _passed_vcf_info_field_evaluation {
    my ($self) = @_;

    my $logger = get_logger("Pathogens::Variant::Evaluator::Pseudosequence");

    $logger->is_debug() && $logger->debug("Evaluating vcf info field values...\n" . $self->_event->info);
    
    my %param;
    my $evaluation_status = 1; #default is 1, i.e. passed the evaluation
   
    foreach my $param_value ( split(";", $self->_event->info)) {
        my ($parameter, $value) = split('=', $param_value);
        $param{$parameter} = $value if ($value);
    }

    if ( exists $param{'MQ'} and $param{'MQ'} < $self->minimum_map_quality) {
        $evaluation_status = 0;
        $self->reporter->inc_counter_failed_map_quality;
    }
    
    if ( exists $param{'AF1'} ) {
        if (not $self->_event->polymorphic) {
            if ($param{'AF1'} >= $self->af1_complement) { #Why is this not ">=" in SH's script?
                $evaluation_status = 0;
                $self->reporter->inc_counter_failed_af1_allele_frequency;
            }
        } else {
            if ($param{'AF1'} < $self->minimum_af1) {
                $evaluation_status = 0;
                $self->reporter->inc_counter_failed_af1_allele_frequency;
            }
        }
    }

    if (exists $param{'PV4'} and not $self->_passed_pv4_evaluation($param{'PV4'}) ) {
        $evaluation_status = 0 
    }

    if (exists $param{'DP4'} and not $self->_passed_dp4_evaluation($param{'DP4'}) ) {
        $evaluation_status = 0 
    }

    return $evaluation_status; 

}

sub _passed_pv4_evaluation {
    my ($self, $pv4string) = @_;

    #Do not test on pv4 if this is not a polymorphic site: 
    return 1 if (not $self->_event->polymorphic);
    
 
    my ( $strand_bias
        , $base_quality_bias
        , $map_bias
        , $tail_distance_bias ) = split(',', $pv4string);
    
    my $evaluation_status = 1;
    
    if ( $strand_bias < $self->minimum_strand_bias ) {
        $evaluation_status = 0;
        $self->reporter->inc_counter_failed_strand_bias;
    }
    if ( $base_quality_bias < $self->minimum_base_quality_bias) {
        $evaluation_status = 0;
        $self->reporter->inc_counter_failed_base_quality_bias;
    }
    if ( $map_bias < $self->minimum_map_bias) {
        $evaluation_status = 0;
        $self->reporter->inc_counter_failed_map_bias;
    }
    if ( $tail_distance_bias < $self->minimum_tail_bias) {
        $evaluation_status = 0;
        $self->reporter->inc_counter_failed_tail_distance_bias;
    }
    
    return $evaluation_status;
}

sub _passed_dp4_evaluation {
    my ($self, $dp4string) = @_;

    $self->_dp4_parser->parse( $dp4string );
    
    my $evaluation_status = 1;
    
    
    if (not $self->_event->polymorphic) { #dealing with NON-polymorphic site (i.e. Vcf's ALT field equals '.')

        #reference depth test for the reference allele
        if ( $self->_dp4_parser->count_referecence_bases < $self->minimum_depth) {
            $self->reporter->inc_counter_failed_depth;
            $evaluation_status = 0;
        }
        
        #forward strand depth test for the reference allele
        if ( $self->_dp4_parser->count_reference_forward_bases < $self->minimum_depth_strand) {
            $self->reporter->inc_counter_failed_depth_forward;
            $evaluation_status = 0;
        }
        
        #reverse strand depth test for the reference allele
        if ( $self->_dp4_parser->count_reference_reverse_bases < $self->minimum_depth_strand) {
            $self->reporter->inc_counter_failed_depth_reverse;
            $evaluation_status = 0;
        }
        
        #forward ratio test for the reference allele
        if ( $self->_dp4_parser->ratio_forward_reference_bases < $self->minimum_ratio) {
            $self->reporter->inc_counter_failed_ratio_forward;
            $evaluation_status = 0;
        }
        #reverse ratio test for the reference allele
        if ( $self->_dp4_parser->ratio_reverse_reference_bases < $self->minimum_ratio) {
            $self->reporter->inc_counter_failed_ratio_reverse;
            $evaluation_status = 0;
        }
    } else { #dealing with polymorphic site
    
        #alternative allele depth test
        if ( $self->_dp4_parser->count_alternative_bases < $self->minimum_depth) {
            $self->reporter->inc_counter_failed_depth;
            $evaluation_status = 0;
        }
        
        #forward strand depth test for alternative allele
        if ( $self->_dp4_parser->count_alternative_forward_bases < $self->minimum_depth_strand) {
            $self->reporter->inc_counter_failed_depth_forward;
            $evaluation_status = 0;
        }
        
        #reverse strand depth test for alternative allele
        if ( $self->_dp4_parser->count_alternative_reverse_bases < $self->minimum_depth_strand) {
            $self->reporter->inc_counter_failed_depth_reverse;
            $evaluation_status = 0;
        }
        
        #forward ratio test for alternative allele
        if ( $self->_dp4_parser->ratio_forward_alternative_bases < $self->minimum_ratio) {
            $self->reporter->inc_counter_failed_ratio_forward;
            $evaluation_status = 0;
        }
        #reverse ratio test for alternative allele
        if ( $self->_dp4_parser->ratio_reverse_alternative_bases < $self->minimum_ratio) {
            $self->reporter->inc_counter_failed_ratio_reverse;
            $evaluation_status = 0;
        }
    }
    return $evaluation_status;
}

sub _is_not_an_indel {
    my ($self) = @_;

    my $logger = get_logger("Pathogens::Variant::Evaluator::Pseudosequence");
    
    my ($ref_allele) = split(',', $self->_event->reference_allele);   #(as SH) taking only the 1st element
    my ($alt_allele) = split(',', $self->_event->alternative_allele); #(as SH) taking only the 1st element
    
    $logger->is_debug() && $logger->debug("Checking if the event is an indel...");
    
    if ( length($ref_allele) > 1
            or
         length($alt_allele) > 1
            or 
         $self->_event->info =~ /INDEL/
       )
    {
        $self->reporter->inc_counter_skipped_indel;

        $logger->is_debug() && $logger->debug("Event was classified as indel!...");

        return 0;
    } else {
        $logger->is_debug() && $logger->debug("Event was classified as NON-indel!...");

        return 1;
    }
}

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

=head1 SEE ALSO

=cut

__PACKAGE__->meta->make_immutable;

1;