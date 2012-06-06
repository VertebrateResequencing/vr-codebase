    package Pathogens::Variant::Utils::PseudoReferenceMaker;

use Moose;
use Pathogens::Variant::Iterator::Vcf;
use Pathogens::Variant::Evaluator::Pseudosequence;
use Pathogens::Variant::EvaluationReporter;
use Pathogens::Variant::Utils::BamParser;
use Pathogens::Variant::Exception;
use namespace::autoclean;

use Log::Log4perl qw(get_logger);







has 'arguments'       => (is => 'rw', isa => 'HashRef', required => 1, trigger => \&_initialise);

#all attributes below are initialised by the trigger in the "arguments" attribute above
has '_iterator'       => (is => 'rw', isa => 'Pathogens::Variant::Iterator::Vcf', init_arg => undef );
has '_reporter'       => (is => 'rw', isa => 'Pathogens::Variant::EvaluationReporter', init_arg => undef );
has '_evaluator'      => (is => 'rw', isa => 'Pathogens::Variant::Evaluator::Pseudosequence', init_arg => undef );
has '_bam_parser'     => (is => 'rw', isa => 'Pathogens::Variant::Utils::BamParser', init_arg => undef );
has '_out_filehandle' => (is => 'rw', isa => 'FileHandle', init_arg => undef );





#################################################
#initialise all the objects/settings for this run
#################################################
sub _initialise {
    
    my ($self) = @_;
    my $logger = get_logger("Pathogens::Variant::Utils::PseudoReferenceMaker");
     
    #complement af1 is needed for non variant site evaluations
    $self->arguments->{af1complement} = 1 - $self->arguments->{af1};
    
    #total depth has to be at least twice larger than strand-depth
    my $doubled_arg_D_depth_strand = $self->arguments->{depth_strand} * 2;
    if ( $self->arguments->{depth} < $doubled_arg_D_depth_strand ) {
        
        $logger->is_warn() && $logger->warn("Argument value for '-d' (depth) must be larger than '-D' (depth_strand)! Automatically increasing it to " .$doubled_arg_D_depth_strand);

        $self->arguments->{depth} = $doubled_arg_D_depth_strand;
    }
    
    #this will create the vcf file traverser
    my $iterator = Pathogens::Variant::Iterator::Vcf->new(vcf_file => $self->arguments->{vcf_file});
    
    #reporter will be used by the evaluator object below to keep counters on evaluation statistics
    my $reporter = Pathogens::Variant::EvaluationReporter->new;
    
    #create an evaluator object to mark the bad calls and to convert heterozygous into homozygous calls 
    my $evaluator = Pathogens::Variant::Evaluator::Pseudosequence->new(
          minimum_depth             => $self->arguments->{depth}
        , minimum_depth_strand      => $self->arguments->{depth_strand}
        , minimum_ratio             => $self->arguments->{ratio}
        , minimum_quality           => $self->arguments->{quality}
        , minimum_map_quality       => $self->arguments->{map_quality}
        , minimum_af1               => $self->arguments->{af1}
        , af1_complement            => $self->arguments->{af1_complement}
        , minimum_ci95              => $self->arguments->{ci95}
        , minimum_strand_bias       => $self->arguments->{strand_bias}
        , minimum_base_quality_bias => $self->arguments->{base_quality_bias}
        , minimum_map_bias          => $self->arguments->{map_bias}
        , minimum_tail_bias         => $self->arguments->{tail_bias}
        , reporter                  => $reporter
    );

    open( my $outfhd, ">", $self->arguments->{out} ) or throw Pathogens::Variant::Exception::File->new({text => "Couldn't open filehandle to write pseudo sequence: " . $!});

    my $bam_parser = Pathogens::Variant::Utils::BamParser->new(bam => $self->arguments->{bam});
    $self->_bam_parser($bam_parser);
    $self->_evaluator($evaluator);
    $self->_iterator($iterator);
    $self->_reporter($reporter);
    $self->_out_filehandle($outfhd);

}


############################################################
#Iterate through the VCF file, judge the variant/non-variant
#and write accordingly a new "pseudoReference" to the disk
#################################################### #######
sub create_pseudo_reference {

    my ($self) = @_;
    my $logger = get_logger("Pathogens::Variant::Utils::PseudoReferenceMaker");
    
    

    my $processed_entries = 0;;
    my $filehandle = $self->_out_filehandle;


    my $event = $self->_iterator->next_event();
    $self->_evaluator->evaluate($event);
    my $last_allele = $self->_get_allele_of_evaluated_event($event);
    my $last_chr = $event->chromosome;
    my $last_pos = $event->position;
    my %seen_chromosomes;
    
    print $filehandle ">" . $last_chr . "\n"; 
    
    while( $event = $self->_iterator->next_event() ) {

        my $curr_chr = $event->chromosome; 
        my $curr_pos = $event->position;

        #indication of duplicate entries in the VCF
        if ( ($last_chr eq $curr_chr) and ($last_pos == $curr_pos) ) {

            #increments duplicate artifact counter
            $self->_reporter->inc_skipped_vcf_duplicate_entry_artifact;

        } else {
 
            $self->_evaluator->evaluate($event);
            my $curr_allele = $self->_get_allele_of_evaluated_event($event);
 
            $self->_write_next_pseudo_reference_allele(   $event, $last_pos, $curr_pos, $last_chr, $curr_chr, $last_allele, $curr_allele);
            
            $last_chr = $curr_chr;
            $last_pos = $curr_pos;
            $last_allele = $curr_allele;
            
        }

        $seen_chromosomes{$last_chr} = 1;
        $processed_entries++;
        $logger->info("Processed $processed_entries sites") unless $processed_entries % 10000;

    }


    #having finished looping through all the variants in the VCF file, see if it is necessary to
    #pad the last chromosome with Ns up to the end of its original size
    my $pad_size = $self->_bam_parser->get_chromosome_size($last_chr) - $last_pos;
    $self->_pad_chromosome_file_with_Ns($pad_size); #pad the end with "N" if necessary
    print $filehandle "\n";


    #For all chromosomes that were not present in the VCF file, we fill 
    #N's in the pseudoreference 
    $self->_fill_unseen_chromosomes_with_Ns(\%seen_chromosomes);

}

sub _fill_unseen_chromosomes_with_Ns {

    my ($self, $seen_chromosomes) = @_;
    my $filehandle = $self->_out_filehandle;
    
    foreach my $chr_name ( $self->_bam_parser->get_chromosome_names ) {
        if (not exists $seen_chromosomes->{$chr_name} ) {
            
            print $filehandle ">" . $chr_name . "\n";
            my $chromosome_length =$self->_bam_parser->get_chromosome_size($chr_name);
            $self->_pad_chromosome_file_with_Ns($chromosome_length);
            print $filehandle "\n";
        }
    }
}

sub _write_next_pseudo_reference_allele {
    
    my (  $self
        , $event
        , $last_pos
        , $curr_pos
        , $last_chr
        , $curr_chr
        , $last_allele
        , $curr_allele
       ) = @_;

    my $filehandle = $self->_out_filehandle;
    my $pad_size   = 0;

    if ($curr_chr eq $last_chr) {  #still at the previous chromosome
        
        #if the last position is not followed immediately by the current position, 
        #pad the file with N's before writing the current allele
        $pad_size = $curr_pos - $last_pos - 1;
        $self->_pad_chromosome_file_with_Ns($pad_size) if ($pad_size > 0);

        #write the current allele to the chromosome
        print $filehandle $curr_allele;

    } else { #at a new Chromosome

        #finish writing to the previous chromosome after the last edits
        print $filehandle $last_allele;
        $pad_size =  $self->_bam_parser->get_chromosome_size($last_chr) - $last_pos;
        #pad the end with "N" if necessary
        $self->_pad_chromosome_file_with_Ns($pad_size) if ($pad_size > 0);
        print $filehandle "\n";


        #start writing to the new chromosome
        print $filehandle ">" . $curr_chr . "\n";
        $pad_size = $curr_pos - 1;
        #pad the begin with "N" if necessary
        $self->_pad_chromosome_file_with_Ns($pad_size) if ($pad_size > 0);
        print $filehandle $curr_allele;

    }
}


sub _pad_chromosome_file_with_Ns {
    
    my ($self, $pad_size) = @_;
    my $fhd = $self->_out_filehandle;
    while ($pad_size > 0) {
        print $fhd "N";
        $pad_size--;
    }
}

sub _get_allele_of_evaluated_event {
    
    my ($self, $event) = @_;
    
    if ( $event->passed_evaluation ) {
        if (not $event->polymorphic) {
            return $event->reference_allele;
        } else {
            return $event->alternative_allele;
        }
    } else {
        if ( $event->was_evaluated ) {
            return 'N';
        } else {
            throw Pathogens::Variant::Exception::ObjectUsage->new({text => 'An event must be evaluated before using _get_allele_of_evaluated_event function on it.'})
        }
    }
}



#####################################################################
#To dump evaluation statistics after the sub 'create_pseudo_reference'
#####################################################################
sub get_statistics_dump {
    my ($self) = @_;
    return $self->_reporter->dump;
}

__PACKAGE__->meta->make_immutable;
1;