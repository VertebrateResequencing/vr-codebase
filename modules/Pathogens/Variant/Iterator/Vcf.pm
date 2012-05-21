package Pathogens::Variant::Iterator::Vcf;
use Moose;

use Pathogens::Variant::Exception;
use Vcf; # Petr Danecek's Vcf module
use Switch;
use Pathogens::Variant::Event::Snp;
use Pathogens::Variant::Event::Dip;
use Pathogens::Variant::Sample;

use namespace::autoclean;






has 'vcf_file' => (is => 'rw', isa => 'Str', required => 1, trigger => \&_vcf_set);
has '_vcf_obj'  => (is => 'rw', isa => 'Vcf');

#update the internal object 
#and parse the vcf header
sub _vcf_set {
    my ($self, $vcf_file) = @_;
    
    my $vcf = Vcf->new( file => $vcf_file );   
    $vcf->parse_header; #required by Petr's module after new()
    $self->_vcf_obj($vcf);  
}

sub next_event {
    my ($self) = @_;
    
    while (my $data_array_ref = $self->_vcf_obj->next_data_array) {
        my $event;  
        if ( $$data_array_ref[7] =~ /INDEL/) {
            $event = Pathogens::Variant::Event::Dip->new(); 
        } else {
            $event = Pathogens::Variant::Event::Snp->new; 
        }
        return $self->_populate_event_attributes($event, $data_array_ref);	
    }     
    return undef;
}

sub close_vcf {
    my $self = shift;
    $self->_vcf_obj->close();
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
