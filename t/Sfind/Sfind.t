
#!/software/bin/perl
use strict;
use warnings;
no warnings 'uninitialized';

use Getopt::Long;
use Sfind::Sfind;
use Test::More qw( no_plan );

use_ok('Sfind::Sfind'); 

ok(testStudyTree(),'Sfind study->sample->library_request->library->sequencing_request->lane->fastq');

sub testStudyTree{
    
#CONNECT TO DATABASE: Warehouse 

my $strack = Sfind::Sfind->new();

unless ($strack){
    die "Can't connect to warehouse database! \n";
}    

#GET the following studies
#my @studies=("Targeted sequencing of plasmodium genes","Trypanosoma brucei");
my @studies=("Targeted sequencing of plasmodium genes");

#my @studies=("African invasive NTS isolates");

#my @studies=("Neospora caninum sequencing");
#my @studies=("P_berghei_ANKA");
#my @studies=("Babesia bigemina sequencing");
#my @studies=("Staphylococcus aureus EARSS1");
#my @studies=("Streptococcus pneumoniae PMEN1");
#my @studies=("Staphylococcus aureus ST239 diversity");
#And then for each study....

foreach my $pname (@studies){
	
    my $study = $strack->get_study_by_name($pname);
    print "study=".$study->id()." ".$study->name()."\n";
    
    my $samples;
    eval {
        $samples = $study->samples;
    };
    if ($@){
        print STDERR "Skipping $pname: Error getting samples: ",$@,"\n";
        next;
    }
	#Samples
    foreach my $sample(@$samples){
	my $sample_id=$sample->id();
	print STDOUT "\tsample_id=$sample_id\n";
	my $library_requests=$sample->library_requests();
	foreach my $lib_request(@$library_requests){
	    
	    print STDOUT "\t\tlibrary_request_id=".$lib_request->id()."\n";
	    print STDOUT "\t\tlibrary_request_status=".$lib_request->status()."\n";
	    
	    my $libraries=$lib_request->libraries();
	    foreach my $library(@$libraries){
		print STDOUT "\t\t\tlibrary_id=".$library->id()."\n";
		print STDOUT "\t\t\tis_tagged=".$library->is_tagged()."\n";
		print STDOUT "\t\t\ttype=".$library->type()."\n";
		my ($frag_from, $frag_to) = @{$library->fragment_size};
		print STDOUT "\t\t\tfragment_size=".$frag_from." - ".$frag_to."\n";
		if($library->is_tagged()){
		    my $multiplex_pool_asset_ids=$library->multiplex_pool_asset_ids();
		    my $tag_group_id=$library->tag_group_id();
		    my $tag_id=$library->tag_id();
		    my $tag_sequence=$library->tag_sequence();
		    foreach my $multiplex_pool_asset_id (@$multiplex_pool_asset_ids){
			print STDOUT "\t\t\tmultiplex_pool_asset_id".$multiplex_pool_asset_id."\ttag_group_id=".$tag_group_id."\ttag_id=".$tag_id."\ttag_sequence=".$tag_sequence."\n";
		    }
		}
		print STDOUT "\t\t\tlibrary_name=".$library->name()."\n";
		
		my $seqrequests=$library->seq_requests();
		foreach my $seqrequest(@$seqrequests){
		print STDOUT "\t\t\t\tseqrequest_id=".$seqrequest->id()."\n";
		print STDOUT "\t\t\t\tlibrary_id=".$seqrequest->library_id()."\n";
		print STDOUT "\t\t\t\tseqrequest_created=".$seqrequest->created()."\n";
		print STDOUT "\t\t\t\tseqrequest_type=".$seqrequest->type()."\n";
		print STDOUT "\t\t\t\tseqrequest_status=".$seqrequest->status()."\n";
		my $is_multiplexed_seq_request=($seqrequest->library_id==$library->id) ? 0 : 1; 
		print STDOUT "\t\t\t\tis_multiplexed_request=".$is_multiplexed_seq_request."\n";
		    
		    my $lanes=$seqrequest->lanes();
		    foreach my $lane(@$lanes){
			    print STDOUT "\t\t\t\t\tlane_id=".$lane->id()."\n";
			    print STDOUT "\t\t\t\t\tlane_name=".$lane->name()."\n";
			    print STDOUT "\t\t\t\t\tlane_created=".$lane->created()."\n";
			    print STDOUT "\t\t\t\t\tlane_run_name=".$lane->run_name()."\n";

			    my $fastqs=$lane->fastq();
			    foreach my $fastq(@$fastqs){
				print STDOUT "\t\t\t\t\t\tfastq_name=".$fastq->name()."\n";
				print STDOUT "\t\t\t\t\t\tfastq_location=".$fastq->location()."\n";
			    }
		    }
		    print STDOUT "\n";
		}
	    }
	  
	}
	print "\n";
    }
}
print "Done\n";
}




