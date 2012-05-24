package Pathogens::Variant::EvaluationReporter;
use Moose;
extends 'Pathogens::Variant::Root';
use namespace::autoclean;

#traits + handles make use of the "Moose::Meta::Attribute::Native::Trait::Counter"
#( http://search.cpan.org/dist/Moose/lib/Moose/Meta/Attribute/Native/Trait/Counter.pm )
has 'failed_base_quality' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_bq_fail => 'inc'
    }
);

has 'failed_map_quality' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_map_quality_fail => 'inc'
    }
);

has 'failed_depth' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_depth_fail => 'inc'
    }
);

has 'failed_depth_forward' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_depth_forward_fail => 'inc'
    }
);

has 'failed_depth_reverse' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_depth_reverse_fail => 'inc'
    }
);

has 'failed_ratio_forward_snp' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_ratio_forward_fail => 'inc'
    }
);

has 'failed_ratio_reverse_snp' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_ratio_reverse_fail => 'inc'
    }
);

has 'failed_allele_frequency' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_af1_fail => 'inc'
    }
);

has 'failed_heterozygous_snp' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_heterozygous_fail => 'inc'
    }
);

has 'failed_strand_bias' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_strand_bias_fail => 'inc'
    }
);

has 'failed_base_quality_bias' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_base_quality_bias_fail => 'inc'
    }
);

has 'failed_map_bias' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_map_bias_fail => 'inc'
    }
);

has 'failed_tail_distance_bias' => ( 
      is => 'ro' 
    , isa => 'Int', 
    , default => 0
    , traits  => ['Counter']
    , handles => {
          inc_tail_bias_fail => 'inc'
    }
);

__PACKAGE__->meta->make_immutable;

1;