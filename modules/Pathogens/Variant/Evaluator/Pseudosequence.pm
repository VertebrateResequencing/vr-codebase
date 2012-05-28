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





has 'minimum_depth'        => ( is => 'rw', isa => 'Int', default => 4 );
has 'minimum_depth_strand' => ( is => 'rw', isa => 'Int', default => 2 );
has 'minimum_ratio'        => ( is => 'rw', isa => 'Num', default => 0.8 );
has 'minimum_quality'      => ( is => 'rw', isa => 'Int', default => 0 );
has 'minimum_map_quality'  => ( is => 'rw', isa => 'Int', default => 1 );
has 'minimum_af1'          => ( is => 'rw', isa => 'Num', default => 0.95 );
has 'minimum_ci95'         => ( is => 'rw', isa => 'Num', default => 0.0 );
has 'minimum_strand_bias'  => ( is => 'rw', isa => 'Num', default => 0.001 );
has 'minimum_base_quality_bias' => ( is => 'rw', isa => 'Num', default => 0.0 );
has 'minimum_map_bias'     => ( is => 'rw', isa => 'Num', default => 0.001 );
has 'minimum_tail_bias'    => ( is => 'rw', isa => 'Num', default => 0.001 );

has '_reporter'             => ( is => 'rw', isa => 'Pathogens::Variant::EvaluationReporter', default => sub { return Pathogens::Variant::EvaluationReporter->new } );
has '_dp4_parser'          => ( is => 'ro', isa => 'Pathogens::Variant::Utils::DP4Parser', lazy => 1, default => sub { return Pathogens::Variant::Utils::DP4Parser->new } );
has '_event_manipulator'   => ( is => 'ro', isa => 'Pathogens::Variant::Utils::EventManipulator', lazy => 1, default => sub { return Pathogens::Variant::Utils::EventManipulator->new } );

sub evaluate {

    my ($self, $event) = @_;
    my $logger = get_logger("Pathogens::Variant::Evaluator::Pseudosequence");

    my $passed_vcf_info_evaluation =  $self->_passed_vcf_info_field_evaluation($event);
    my $is_not_indel               =  $self->_is_not_an_indel($event);
    my $is_heterozygous            =  $self->_has_secondary_heterozygous_alternative_alleles($event);
    
    if ( $passed_vcf_info_evaluation and $is_not_indel ) {
        if ( $is_heterozygous ) {
            #the line below, modifies the heterozygous event from its original state
            $self->_event_manipulator->remove_secondary_alternative_heterozygous_alleles($event);
        }
        
        $event->passed_evaluation(1);
        $self->_reporter->inc_counter_accepted_snp_calls;
        
        $logger->debug("Event dump after passing the evaluation:...\n". Dumper $event);
        $logger->debug("Reporter dump after passing the evaluation:...\n". Dumper $self->_reporter);
       
        return 1;

    } else {

        $event->passed_evaluation(0); #i.e. the event failed to pass the filters

        $logger->debug("Event dump after failing the evaluation:...\n". Dumper $event);
        $logger->debug("Reporter dump after failing the evaluation:...\n". Dumper $self->_reporter);

        return 0;

    }
}

sub _has_secondary_heterozygous_alternative_alleles {
    
    my ($self, $event) = @_;
    
    my @num_alleles = split(',', $event->alternative_allele);
    
    if (scalar @num_alleles > 1) {
        $self->_reporter->inc_counter_heterozygous_calls;
        return 1;
    } else {
        return 0;
    }

}
sub _passed_vcf_info_field_evaluation {

    my ($self, $event) = @_;
    my $logger = get_logger("Pathogens::Variant::Evaluator::Pseudosequence");

    my %param;
    my $evaluation_status = 1; #default is 1, i.e. passed the evaluation
   
    foreach my $param_value ( split(";", $event->info)) {
        my ($parameter, $value) = split('=', $param_value);
        $param{$parameter} = $value;
    }

    $logger->debug("Evaluating vcf info field values...\n" . Dumper %param);

    if ( exists $param{'MQ'} and $param{'MQ'} < $self->minimum_map_quality) {
        $evaluation_status = 0;
        $self->_reporter->inc_counter_failed_map_quality;
    }
    
    if ( exists $param{'AF1'} and $param{'AF1'} < $self->minimum_af1) {
        $evaluation_status = 0;
        $self->_reporter->inc_counter_failed_af1_allele_frequency;
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

    my ( $strand_bias
        , $base_quality_bias
        , $map_bias
        , $tail_distance_bias ) = split(',', $pv4string);
    my $evaluation_status = 1;
    
    if ( $strand_bias < $self->minimum_strand_bias ) {
        $evaluation_status = 0;
        $self->_reporter->inc_counter_failed_strand_bias;
    }
    if ( $base_quality_bias < $self->minimum_base_quality_bias) {
        $evaluation_status = 0;
        $self->_reporter->inc_counter_failed_base_quality_bias;
    }
    if ( $map_bias < $self->minimum_map_bias) {
        $evaluation_status = 0;
        $self->_reporter->inc_counter_failed_map_bias;
    }
    if ( $tail_distance_bias < $self->minimum_tail_bias) {
        $evaluation_status = 0;
        $self->_reporter->inc_counter_failed_tail_distance_bias;
    }
    
    return $evaluation_status;
}

sub _passed_dp4_evaluation {
    
    my ($self, $dp4string) = @_;
    
    my $parsed = $self->_dp4_parser->parse( $dp4string );
    my $evaluation_status = 1;
    
    #total depth tests
    if ( $parsed->{count_alternative_bases} < $self->minimum_depth) {
        $self->_reporter->inc_counter_failed_depth;
        $evaluation_status = 0;
    }
    
    #strand depth tests
    if ( $parsed->{count_alternative_forward_bases} < $self->minimum_depth_strand) {
        $self->_reporter->inc_counter_failed_depth_forward;
        $evaluation_status = 0;
    }
    
    if ( $parsed->{count_alternative_reverse_bases} < $self->minimum_depth_strand) {
        $self->_reporter->inc_counter_failed_depth_reverse;
        $evaluation_status = 0;
    }
    
    #ratios tests
    if ( $parsed->{ratio_forward_alternative_bases} < $self->minimum_ratio) {
        $self->_reporter->inc_counter_failed_ratio_forward;
        $evaluation_status = 0;
    }
    if ( $parsed->{ratio_reverse_alternative_bases} < $self->minimum_ratio) {
        $self->_reporter->inc_counter_failed_ratio_reverse;
        $evaluation_status = 0;
    }
    
    return $evaluation_status;
}

sub _is_not_an_indel {
    
    my ($self, $event) = @_;
    my $logger = get_logger("Pathogens::Variant::Evaluator::Pseudosequence");
    
    my ($ref_allele) = split(',', $event->reference_allele);   #taking only the 1st element
    my ($alt_allele) = split(',', $event->alternative_allele); #taking only the 1st element
    
    $logger->debug("Checking if the event is an indel...");
    
    if ( length($ref_allele) > 1
            or
         length($alt_allele) > 1
            or
         $event->info =~ /INDEL/ )
    {
        $self->_reporter->inc_counter_skipped_indel;

        $logger->debug("Event was classified as indel!...");

        return 0;
    } else {
        $logger->debug("Event was classified as NON-indel!...");

        return 1;
    }
}

sub dump_evaluation_statistics {
    my ($self) = @_;
    return Dumper $self->_reporter;
}

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

=head1 SEE ALSO

=cut

__PACKAGE__->meta->make_immutable;

1;