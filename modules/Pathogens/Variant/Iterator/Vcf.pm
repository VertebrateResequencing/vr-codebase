=head1 NAME

Pathogens::Variant::Iterator::Vcf  - A Vcf file traverser that returns an iterator. The iterator delivers "Pathogens::Variant::Event" objects at each iteration

=head1 SYNOPSIS

#iterates through a vcf
$object = Pathogens::Variant::Iterator::Vcf->new (vcf_file => "$dir/data/single.snp.vcf.gz");

while( my $next_event = $object->next_event ) {
   #do something with "$next_event"...
}

=cut

package Pathogens::Variant::Iterator::Vcf;
use Moose;
extends 'Pathogens::Variant::Root';

use Vcf; # Petr Danecek's Vcf module
use Pathogens::Variant::Event;
use Pathogens::Variant::Sample;

use namespace::autoclean;






has 'vcf_file' => (is => 'rw', isa => 'Str', required => 1, trigger => \&_vcf_set);
has '_vcf_obj'  => (is => 'rw', isa => 'Vcf');


#update the internal object 
#and parse the vcf header
sub _vcf_set {
    my ($self, $vcf_file) = @_;
    #we will wrap our iterator on top of Petr's Vcf module
    my $vcf = Vcf->new( file => $vcf_file );
    $vcf->recalc_ac_an(0); #turn this off. See "perldoc Vcf"
    $vcf->parse_header; #required by Petr's module after new(), See "perldoc Vcf"
    $self->_vcf_obj($vcf);  
}

sub next_event {

    my ($self) = @_;
    
    my $data_array_ref = $self->_vcf_obj->next_data_array;

    if ($data_array_ref) {
        
        my $event = Pathogens::Variant::Event->new();
        
        return $self->_populate_event_attributes($event, $data_array_ref);
        
    } else {
        
        $self->_vcf_obj->close();
        return undef;
        
    }
}

sub _populate_event_attributes {
    my ($self, $event, $data_array_ref) = @_;

    #there are 8 fixed fields in a vcf (http://www.1000genomes.org/node/101)
    $event->chromosome($$data_array_ref[0]);
    $event->position($$data_array_ref[1]);
    $event->id($$data_array_ref[2]);
    $event->reference_allele($$data_array_ref[3]);
    $event->alternative_allele($$data_array_ref[4]);
    $event->quality($$data_array_ref[5]);
    $event->filter($$data_array_ref[6]);
    $event->info($$data_array_ref[7]);
    $event->format($$data_array_ref[8]);

    #Vcf.pm's get_samples returns a sorted list of the sample names as they appear
    #starting on the 9th (0-based count) field of a Vcf file
    my $data_array_counter = 9; 
    foreach my $sample_name ($self->_vcf_obj->get_samples) {
        my $sample_object = Pathogens::Variant::Sample->new(name => $sample_name);
        $sample_object->genotype($$data_array_ref[$data_array_counter]);
        $event->add_sample($sample_object);
        $data_array_counter++;
    }
    return $event;
}

__PACKAGE__->meta->make_immutable;

1;
