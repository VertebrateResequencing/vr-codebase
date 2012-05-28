package Pathogens::Variant::EvaluationReporter;
use Moose;
extends 'Pathogens::Variant::Root';
use namespace::autoclean;



#Calling "inc_counter_failed_base_quality" on the object will increment 
#the counter of the attribute "failed_base_quality"
#( http://search.cpan.org/dist/Moose/lib/Moose/Meta/Attribute/Native/Trait/Counter.pm )
has 'failed_base_quality' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_base_quality => 'inc'  
    }
);

has 'failed_map_quality' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_map_quality => 'inc'
    }
);

has 'failed_depth' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_depth => 'inc'
    }
);

has 'failed_depth_forward' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_depth_forward => 'inc'
    }
);

has 'failed_depth_reverse' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_depth_reverse => 'inc'
    }
);

has 'failed_ratio_forward' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_ratio_forward => 'inc'
    }
);

has 'failed_ratio_reverse' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_ratio_reverse => 'inc'
    }
);

has 'failed_af1_allele_frequency' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_af1_allele_frequency => 'inc'
    }
);

has 'failed_heterozygous_snp' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_heterozygous_snp => 'inc'
    }
);

has 'failed_strand_bias' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_strand_bias => 'inc'
    }
);

has 'failed_base_quality_bias' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_base_quality_bias => 'inc'
    }
);

has 'failed_map_bias' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_map_bias => 'inc'
    }
);

has 'failed_tail_distance_bias' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_failed_tail_distance_bias => 'inc'
    }
);


has 'skipped_indel' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_counter_skipped_indel => 'inc'
    }
);
__PACKAGE__->meta->make_immutable;

1;