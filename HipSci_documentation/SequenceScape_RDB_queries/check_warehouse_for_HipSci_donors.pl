#!/bin/perl -w

use strict;

my (
    @cohortIDs, 
    $thisCohortID
    );

my $cohortsFile = "currentHipSciCohorts.txt";

my $cohortsQuery = "mysql -u warehouse_ro --host mcs7 --port 3379 sequencescape_warehouse --default-character-set=utf8 -e \"
select distinct current_samples.donor_id from current_samples, current_studies, current_study_samples where current_studies.internal_id = 2624 and (current_samples.uuid = current_study_samples.sample_uuid) and (current_study_samples.study_uuid = current_studies.uuid);
\" > $cohortsFile
";
system("$cohortsQuery");

open(COHORTS, "$cohortsFile") || die "cannot open $cohortsFile ($!)";
<COHORTS>;
while(<COHORTS>){
    chomp;
    if($_ ne "NULL"){
        push(@cohortIDs, $_);
    }
}
close(COHORTS);

my $count = 0;

foreach(@cohortIDs){
    $count++;

    my $firstPartOfQuery = "mysql -u warehouse_ro --host mcs7 --port 3379 sequencescape_warehouse --default-character-set=utf8 -e \"
select current_samples.public_name as 'Friendly ID', current_samples.supplier_name as 'Supplier Name', current_samples.sanger_sample_id as 'Sanger Sample ID', current_samples.control as 'Control Status'
from current_samples, current_study_samples, current_studies 
where 
current_samples.donor_id = '$_' 
and 
( 
(current_samples.uuid = current_study_samples.sample_uuid) 
and 
(current_study_samples.study_uuid = current_studies.uuid)
) ";

    my $genotypingQuery = $firstPartOfQuery . "and 
(current_studies.internal_id = 2624) 
order by current_samples.public_name;
\"
";

    my $geneexpressionQuery = $firstPartOfQuery . "and 
(current_studies.internal_id = 2625) 
order by current_samples.public_name;
\"
";

    print "Donor #$count \"$_\"\n";
    print "Genotyping info:\n";
    #print("$genotypingQuery");
    system("$genotypingQuery");
    print "\n";

    print "Gene expression info:\n";
    #print("$geneexpressionQuery");
    system("$geneexpressionQuery");
    print "\n\n";
}
