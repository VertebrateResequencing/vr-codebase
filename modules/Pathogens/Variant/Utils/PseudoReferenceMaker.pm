=head1 NAME

Pathogens::Variant::Utils::PseudoReferenceMaker  - Creates a Pseudoreference. The module required a bam file, a reference fasta file, a lane name and a output file name. It can take a list of filtering options for evaluation of each snp or reference base. The results of the evaluation will shape the pseudo-reference.

=head1 SYNOPSIS

#requires SAMTOOLS shell environment variable to point to samtools compiled dir
  
my %args = (

    #REQUIRED
     out            => $pseudo_reference_output
   , bam            => "$dir/data/test.bam"
   , lane_name      => "test_lane"
   , reference      => "$dir/data/reference.snpAt10.insertionAt30.deletionAt50.fasta"
      
   #OPTIONALY CHANGE TO OTHER VALUES FOR MORE/LESS STRINGENT FILTERING
   , depth          => 4
   , depth_strand   => 2
   , ratio          => 0.75
   , quality        => 50
   , map_quality    => 30
   , af1            => 0.95
   , ci95           => 0.0
   , strand_bias    => 0.001
   , base_quality_bias => 0.0
   , map_bias       => 0.001
   , tail_bias      => 0.001

);

my $object = Pathogens::Variant::Utils::PseudoReferenceMaker->new(arguments => \%args );

#writes the pseudoreference (in fasta format) to a file named "$pseudo_reference_output"
$object->create_pseudo_reference;

=cut

package Pathogens::Variant::Utils::PseudoReferenceMaker;

use Moose;
use Pathogens::Variant::Iterator::Vcf;
use Pathogens::Variant::Evaluator::Pseudosequence;
use Pathogens::Variant::EvaluationReporter;
use Pathogens::Variant::Utils::BamParser;
use Utils;
use Pathogens::Variant::Exception;
use Log::Log4perl qw(get_logger);
use File::Basename;
use namespace::autoclean;





has 'arguments'       => (is => 'rw', isa => 'HashRef', required => 1, trigger => \&_initialise);

#all attributes below are initialised by the trigger in the "arguments" attribute above
has '_reporter'       => (is => 'rw', isa => 'Pathogens::Variant::EvaluationReporter', init_arg => undef );
has '_evaluator'      => (is => 'rw', isa => 'Pathogens::Variant::Evaluator::Pseudosequence', init_arg => undef );
has '_bam_parser'     => (is => 'rw', isa => 'Pathogens::Variant::Utils::BamParser', init_arg => undef );
has '_output_filehandle_temporary_file' => (is => 'rw', isa => 'FileHandle', init_arg => undef );


######################################################################
#initialise some objects and settings according to the given arguments
######################################################################
sub _initialise {
    my ($self) = @_;

    my $logger = get_logger("Pathogens::Variant::Utils::PseudoReferenceMaker");
     
    #reporter will be used by the evaluator object below to hold evaluation statistics
    my $reporter = Pathogens::Variant::EvaluationReporter->new;
    
    #the evaluator will help us decide, whether to keep an allele in the new pseudo-genome
    #by looking if the quality of a position in the VCF file is good enough
    my $evaluator = Pathogens::Variant::Evaluator::Pseudosequence->new;

    $evaluator->reporter($reporter);
    
    $evaluator->minimum_depth($self->arguments->{depth})  
        if ( exists $self->arguments->{depth} ); 

    $evaluator->minimum_depth_strand($self->arguments->{depth_strand})
        if ( exists $self->arguments->{depth_strand} );
 
    $evaluator->minimum_ratio($self->arguments->{ratio})
        if ( exists $self->arguments->{ratio} );

    $evaluator->minimum_quality($self->arguments->{quality})
        if ( exists $self->arguments->{quality} );

    $evaluator->minimum_map_quality($self->arguments->{map_quality})
        if ( exists $self->arguments->{map_quality} );

    $evaluator->minimum_af1($self->arguments->{af1}) 
        if ( exists $self->arguments->{af1} );

    $evaluator->minimum_ci95($self->arguments->{ci95}) 
        if ( exists $self->arguments->{ci95} );

    $evaluator->minimum_strand_bias($self->arguments->{strand_bias}) 
        if ( exists $self->arguments->{strand_bias} );

    $evaluator->minimum_base_quality_bias($self->arguments->{base_quality_bias})
        if ( exists $self->arguments->{base_quality_bias} );

    $evaluator->minimum_map_bias($self->arguments->{map_bias})
        if ( exists $self->arguments->{map_bias} );

    $evaluator->minimum_tail_bias($self->arguments->{tail_bias})
        if ( exists $self->arguments->{tail_bias} );
    

    my $bam_parser = Pathogens::Variant::Utils::BamParser->new();
    $bam_parser->fetch_chromosome_size_into_hash($self->arguments->{bam});
    $self->_bam_parser($bam_parser);
    
    $self->_evaluator($evaluator);
    
    $self->_reporter($reporter);

}

############################################################
#Iterate through the VCF file, judge the variant/non-variant
#and write accordingly a new "pseudoReference" to the disk
#################################################### #######
sub create_pseudo_reference {
    my ($self) = @_;

    my $logger = get_logger("Pathogens::Variant::Utils::PseudoReferenceMaker");


    $logger->info("Creating a VCF file with all mapped sites");
    my $temporary_vcf_file_with_all_mapped_positions = $self->_generate_vcf_file_with_all_reference_sites;
    $logger->info("VCF has been created.");


    my $temporary_tab_delimited_pseudo_reference_file = $self->_generate_tab_delimited_pseudo_reference_file($temporary_vcf_file_with_all_mapped_positions);
    
    
    #The entries in the $temporary_tab_delimited_pseudo_reference_file are first sorted (dictionary sort) by sequence name
    #and the sequences are then concatenated end-to-end into a single string, and saved in fasta format.
    $self->_sort_and_merge_and_generate_final_output($temporary_tab_delimited_pseudo_reference_file);

    unlink($temporary_tab_delimited_pseudo_reference_file);
    unlink($temporary_vcf_file_with_all_mapped_positions);

}

sub _generate_tab_delimited_pseudo_reference_file {
    my ($self, $temporary_vcf_file_with_all_mapped_positions) = @_;

    my $logger = get_logger("Pathogens::Variant::Utils::PseudoReferenceMaker");
    
    #this will create the vcf file traverser
    my $iterator = Pathogens::Variant::Iterator::Vcf->new(vcf_file => $temporary_vcf_file_with_all_mapped_positions);


    #This is the temporary file where we will write the pseudo-sequences in a random order in first stage.
    #This temporary file will be tab delimited with two fields: "SequenceId" TAB "PseudosequenceThatBelongsToTheSequenceID"
    #The file will later be sorted by sequenceID and converted into a final fasta formatted file
    my ($file,$path,$suffix) = File::Basename::fileparse($self->arguments->{out}); # returns  ("baz", "/foo/bar/", ".txt") 
    my $temporary_tab_delimited_pseudo_reference_file = "$path/temporary.$file";

    #open the filehandle for the temporary file
    my $filehandle;
    unless ( open( $filehandle, ">" . $temporary_tab_delimited_pseudo_reference_file ) ) {
        $logger->error("Couldn't open filehandle to write the pseudoreference. " . $!); 
        throw Pathogens::Variant::Exception::File({text => "Couldn't open filehandle to write the pseudoreference. " . $!});
    }
    
    #set it here, so that we can access this file handle from other subs
    $self->_output_filehandle_temporary_file($filehandle);
   
    #this hash will be used to keep track of observed chromosomes in the VCF file
    my %seen_chromosomes;
    
    ######################################################
    #The very first event in the VCF file is handled here:
    my $event = $iterator->next_event();
    $self->_evaluator->evaluate($event);
    my $last_allele = $self->_get_allele_of_evaluated_event($event);
    my $last_chr = $event->chromosome;
    my $last_pos = $event->position;
    
    #write the name of the new chromosome to the file
    print $filehandle $last_chr . "\t";
    
    #we may not be at the first position of this chromosome, we fill any preeceding 
    #positions with N's in order to keep identical chromosome size to the reference
    my $initial_pad_size = $last_pos - 1;
    $self->_pad_chromosome_file_with_Ns($initial_pad_size) if ($initial_pad_size > 0); #writes indirectly to the $filehandle   
    
    #now we are ready to append the allele to the file:
    print $filehandle $last_allele;  
    
    #keep track of the chromosomes that we have seen in this vcf file
    $seen_chromosomes{$last_chr} = 1; 
    my $processed_entries = 1; 

    ######################################################
    #All other events in the VCF are handled here:########
    while( $event = $iterator->next_event() ) {

        my $curr_chr = $event->chromosome; 
        my $curr_pos = $event->position;

        #indication of duplicate entries in the VCF (this is strange, but is observed sometimes)
        if ( ($last_chr eq $curr_chr) and ($last_pos == $curr_pos) ) {

            #increments duplicate artifact counter
            $self->_reporter->inc_skipped_vcf_duplicate_entry_artifact;

        #a unique entry in the vcf file
        } else {

            $self->_evaluator->evaluate($event);
            my $curr_allele = $self->_get_allele_of_evaluated_event($event);
            $self->_write_next_pseudo_reference_allele($event, $last_pos, $curr_pos, $last_chr, $curr_chr, $last_allele, $curr_allele);
            $last_allele = $curr_allele;

        }

        #reset for the next round in the loop
        $last_chr = $curr_chr;
        $last_pos = $curr_pos;


        #so that we know which ones have been seen in the vcf file
        $seen_chromosomes{$last_chr} = 1; 

        #optinal to monitor the progress whenever log-level is set to "INFO" (or lower)
        $processed_entries++;
        $logger->info("Processed $processed_entries sites") unless $processed_entries % 10000;

    }

    #Pads the last sequence of pseudoreference with Ns to have identical size with the original reference
    my $pad_size = $self->_bam_parser->get_chromosome_size($last_chr) - $last_pos;
    $self->_pad_chromosome_file_with_Ns($pad_size); #pad the end with "N" if necessary
    print $filehandle "\n";


    #For all chromosomes that were not present in the VCF file, we fill 
    #N's in the pseudoreference to have identical chromosome set with the original reference
    $self->_fill_unseen_chromosomes_with_Ns(\%seen_chromosomes);

    close($filehandle);

    return $temporary_tab_delimited_pseudo_reference_file;

}

sub _fill_unseen_chromosomes_with_Ns {
    my ($self, $seen_chromosomes) = @_;

    my $filehandle = $self->_output_filehandle_temporary_file;
    
    foreach my $chr_name ( $self->_bam_parser->get_chromosome_names ) {
        if (not exists $seen_chromosomes->{$chr_name} ) {
            
            print $filehandle $chr_name . "\t";
            my $chromosome_length =$self->_bam_parser->get_chromosome_size($chr_name);
            $self->_pad_chromosome_file_with_Ns($chromosome_length);
            print $filehandle "\n";
        }
    }
}

sub _write_next_pseudo_reference_allele {
    my (  $self, $event, $last_pos, $curr_pos, $last_chr, $curr_chr, $last_allele, $curr_allele) = @_;

    my $filehandle = $self->_output_filehandle_temporary_file;
    my $pad_size   = 0;

    if ($curr_chr eq $last_chr) {  #still at the previous chromosome
        
        #if the last position is not followed immediately by the current position, 
        #pad the file with N's before writing the current allele
        $pad_size = $curr_pos - $last_pos - 1;
        $self->_pad_chromosome_file_with_Ns($pad_size) if ($pad_size > 0);

        #write the current allele to the chromosome
        print $filehandle $curr_allele;

    } else { #at a new Chromosome

        #pad the end with "N" if necessary
        $pad_size =  $self->_bam_parser->get_chromosome_size($last_chr) - $last_pos;
        $self->_pad_chromosome_file_with_Ns($pad_size) if ($pad_size > 0);
        print $filehandle "\n";


        #start writing to the new chromosome
        print $filehandle $curr_chr . "\t";
        
        #pad the begin with "N" if this is not the 1st position
        $pad_size = $curr_pos - 1;
        $self->_pad_chromosome_file_with_Ns($pad_size) if ($pad_size > 0);
        print $filehandle $curr_allele;

    }
}

sub _pad_chromosome_file_with_Ns {
    my ($self, $pad_size) = @_;

    my $filehandle = $self->_output_filehandle_temporary_file;
    while ($pad_size > 0) {
        print $filehandle "N";
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
            my $logger = get_logger("Pathogens::Variant::Utils::PseudoReferenceMaker");
            $logger->error('An un-evaluated reference position was asked to return an allele. You should first evaluate the position');
            throw Pathogens::Variant::Exception::ObjectUsage({text => 'An un-evaluated reference position was asked to return an allele. You should first evaluate the position'});
        }
    }
}

sub _sort_and_merge_and_generate_final_output {
    my ($self, $temporary_tab_delimited_pseudo_reference_file) = @_;
    my $logger = get_logger("Pathogens::Variant::Utils::PseudoReferenceMaker");

    my $final_output_file =  $self->arguments->{out};
    
    my $fhd_file;
    unless ( open($fhd_file, ">" . $final_output_file) ) {
        $logger->error("Could not open output filehandle to write pseudo_reference. $!");
        throw Pathogens::Variant::Exception::File({text => "Could not open output filehandle to write pseudo_reference" . $!});
    }

    my $pipecmd = "set -o pipefail; sort -d -k1,1 $temporary_tab_delimited_pseudo_reference_file | cut -f 2 | tr -d '\n' | fold -w 60 |";
    $logger->debug($pipecmd);
    
    my $fhd_pipe;
    unless (open($fhd_pipe, $pipecmd) ) {
         $logger->error("Pipe failed with run status ($?) when running command: $pipecmd $!");
         throw Pathogens::Variant::Exception::CommandExecution({text => "Pipe failed with run status ($?) when running command: $pipecmd $!"}); 
    }


    #final fasta file
    print $fhd_file ">" . $self->arguments->{lane_name} . "\n";
    while(<$fhd_pipe>) {
         print $fhd_file $_;
    }
    print $fhd_file "\n";
    close $fhd_file;
    
    #most pipe errors are only caught to the parent when closing the pipe handle
    unless (close $fhd_pipe) {
         $logger->error("Pipe failed with run status ($?) when running command: $pipecmd $!");
         throw Pathogens::Variant::Exception::CommandExecution({text => "Pipe failed with run status ($?) when running command: $pipecmd $!"});
    }
}

sub _generate_vcf_file_with_all_reference_sites {
    
    my ($self) = @_;
    my $logger = get_logger("Pathogens::Variant::Utils::PseudoReferenceMaker");
    
    my $reference_file = $self->arguments->{reference};
    
    my $bam_file   = $self->arguments->{bam};
    my ($file,$path,$suffix) = File::Basename::fileparse($bam_file); # returns  ("baz", "/foo/bar/", ".txt") 
    my $temporary_vcf_file_with_all_mapped_positions = "$path/temporary.$file.vcf";


    my $pipecmd = qq[samtools mpileup -d 1000 -DSugBf $reference_file $bam_file | bcftools view -cg - > $temporary_vcf_file_with_all_mapped_positions];
    $logger->debug($pipecmd);

    Utils::CMD($pipecmd);
    
    return $temporary_vcf_file_with_all_mapped_positions;

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
