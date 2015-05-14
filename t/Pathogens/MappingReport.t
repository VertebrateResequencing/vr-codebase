#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most;
    eval {
        require VRTrack::Testconfig;
    };
    if ($@) {
        plan skip_all => "Skipping all tests because VRTrack tests have not been configured";
    }
    
    # Reports modules
    use_ok('Pathogens::Reports::Mapping::Report');
    use_ok('Pathogens::Reports::Mapping::Row');
    use_ok('Pathogens::Reports::Mapping::Spreadsheet');

    # Database modules
    use VRTrack::VRTrack;
    use VRTrack::Project; # Project
    use VRTrack::Sample; # Sample
    use VRTrack::Library_request;
    use VRTrack::Library;
    use VRTrack::Seq_request;
    use VRTrack::Lane;
    use VRTrack::File;
    use VRTrack::Mapstats;
    use VRTrack::Mapper; # Mapper 
    use VRTrack::Assembly; # Assembly
    use VRTrack::Image; # Image
}

my $connection_details = { database => VRTrack::Testconfig->config('test_db'),
                           host     => VRTrack::Testconfig->config('host'),
                           port     => VRTrack::Testconfig->config('port'),
                           user     => VRTrack::Testconfig->config('user'),
                           password => VRTrack::Testconfig->config('password') };


# Create test database.

# confirm vrtrack set
ok my $vrtrack = VRTrack::VRTrack->new($connection_details), 'vrtrack set';
isa_ok $vrtrack, 'VRTrack::VRTrack';

# cleanup
clean_test_db($vrtrack);

# Test data.
my %project = ( ssid => 1000, name => 'My Pet Tyrannosaurus' );
my %sample  = ( ssid => 2000, name => 'Sparky' );
my %library = ( ssid => 3000, name => 'Spark_001' );
my %lane    = ( name => '1234_5#6', processed => 7, qc_status => 'passed',
		readlen => 76, paired => 1, raw_reads => 10066330, raw_bases => 765041080);
my %assembly = ( name => 'Sue_FMNH_PR_2081', reference_size => 23298589);
my %mapper   = ( name_a => 'velocimapper', version_a => '0.0.1',
		 name_b => 'dinomap',      version_b => '0.0.2' );

my %mapstats_qc  = (is_qc => 1, raw_reads => 1438046, raw_bases => 109291496, clip_bases => 108674266,
		    reads_mapped => 1423575, reads_paired => 1393135, bases_mapped  => 107477960, 
		    rmdup_reads_mapped => 1423083, rmdup_bases_mapped => 107443038,
		    adapter_reads => 83, error_rate => 0.00297524, mean_insert  => 315.6, sd_insert => 42.6,
		    target_bases_mapped => 17238015, mean_target_coverage => 4.35, target_coverage_sd => 5.37);

my %mapstats_map = (raw_reads => 10066330, raw_bases => 765041080, 
		    reads_mapped => 9613107, reads_paired => 9268974, bases_mapped  => 726737126,
		    mean_insert => 319, sd_insert => 50.99,
		    target_bases_mapped => 20546772, mean_target_coverage=> 30.16, target_coverage_sd => 37.12,
		    target_bases_1X   => 88.19, target_bases_2X  => 81.15, target_bases_5X  => 70.42, 
		    target_bases_10X  => 62.22, target_bases_20X => 52.09, target_bases_50X => 25.05,
		    target_bases_100X => 2.22);

my %vrobject;
$vrobject{project} = VRTrack::Project->create($vrtrack, $project{name});
$vrobject{project}->ssid($project{ssid});
$vrobject{project}->update();

$vrobject{sample} = VRTrack::Sample->create($vrtrack, $sample{name});
$vrobject{sample}->ssid($sample{ssid});
$vrobject{sample}->project_id($vrobject{project}->id);
$vrobject{sample}->update();

$vrobject{library} = VRTrack::Library->create($vrtrack, $library{name});
$vrobject{library}->ssid($library{ssid});
$vrobject{library}->sample_id($vrobject{sample}->id);
$vrobject{library}->update();

$vrobject{lane} = VRTrack::Lane->create($vrtrack, $lane{name});
$vrobject{lane}->library_id($vrobject{library}->id);
$vrobject{lane}->read_len($lane{readlen});
$vrobject{lane}->raw_reads($lane{raw_reads});
$vrobject{lane}->raw_bases($lane{raw_bases});
$vrobject{lane}->processed($lane{processed});
$vrobject{lane}->qc_status($lane{qc_status});
$vrobject{lane}->update();

$vrobject{assembly} = VRTrack::Assembly->create($vrtrack, $assembly{name});
$vrobject{assembly}->reference_size($assembly{reference_size});
$vrobject{assembly}->update();

$vrobject{mapper_a} = VRTrack::Mapper->create($vrtrack, $mapper{name_a}, $mapper{version_a});
$vrobject{mapper_b} = VRTrack::Mapper->create($vrtrack, $mapper{name_b}, $mapper{version_b});

$vrobject{mapstats_qc} = VRTrack::Mapstats->create($vrtrack, $vrobject{lane}->id()); # qc mapstats
$vrobject{mapstats_qc}->lane_id($vrobject{lane}->id());
$vrobject{mapstats_qc}->mapper_id($vrobject{mapper_a}->id);
$vrobject{mapstats_qc}->assembly_id($vrobject{assembly}->id);
for my $field (keys %mapstats_qc)
{
    $vrobject{mapstats_qc}->$field($mapstats_qc{$field});
}
$vrobject{mapstats_qc}->update();

$vrobject{image} = VRTrack::Image->create($vrtrack, 'graph', ':)'); # qc graph image
$vrobject{image}->mapstats_id($vrobject{mapstats_qc}->id());
$vrobject{image}->update();

$vrobject{mapstats_map} = VRTrack::Mapstats->create($vrtrack, $vrobject{lane}->id()); # mapstats
$vrobject{mapstats_map}->lane_id($vrobject{lane}->id());
$vrobject{mapstats_map}->mapper_id($vrobject{mapper_b}->id);
$vrobject{mapstats_map}->assembly_id($vrobject{assembly}->id);
for my $field (keys %mapstats_map)
{
    $vrobject{mapstats_map}->$field($mapstats_map{$field});
}
$vrobject{mapstats_map}->update();

$vrobject{mapstats_blank} = VRTrack::Mapstats->create($vrtrack, $vrobject{lane}->id()); # blank/failed
$vrobject{mapstats_blank}->lane_id($vrobject{lane}->id());
$vrobject{mapstats_blank}->mapper_id($vrobject{mapper_b}->id);
$vrobject{mapstats_blank}->assembly_id($vrobject{assembly}->id);
$vrobject{mapstats_blank}->update();


# Test Rows
# blank row
ok my $row_blank = Pathogens::Reports::Mapping::Row->new(vrtrack => $vrtrack, lane => $vrobject{lane}, mapstats => $vrobject{mapstats_blank}), 'Read lane and blank mapstats to row.';

is $row_blank->study_id,                1000,'study_id ok';
is $row_blank->sample,              'Sparky','sample ok';
is $row_blank->lanename,          '1234_5#6','lane ok';
is $row_blank->cycles,                    76,'cycles ok';
is $row_blank->reads,               10066330,'yield reads ok';
is $row_blank->bases,              765041080,'yield bases ok';
is $row_blank->map_type,           'Mapping','map_type ok';
is $row_blank->reference, 'Sue_FMNH_PR_2081','reference ok';
is $row_blank->reference_size,      23298589,'reference_size ok';
is $row_blank->mapper,             'dinomap','mapper ok';
is $row_blank->mapstats_id,                5,'mapstats_id ok'; # based on order mapstats created
is $row_blank->adapter_perc,           undef,'adapter ok';
is $row_blank->transposon_perc,        undef,'transposon ok';
is $row_blank->mapped_perc,            '0.0','mapped ok';
is $row_blank->paired_perc,            '0.0','paired ok';
is $row_blank->mean_insert_size,       undef,'mean_insert_size ok';
is $row_blank->genome_covered,         undef,'genome_covered ok';
is $row_blank->genome_covered_1x,      undef,'genome_covered_1x   ok';
is $row_blank->genome_covered_5x,      undef,'genome_covered_5x   ok';
is $row_blank->genome_covered_10x,     undef,'genome_covered_10x  ok';
is $row_blank->genome_covered_50x,     undef,'genome_covered_50x  ok';
is $row_blank->genome_covered_100x,    undef,'genome_covered_100x ok';
is $row_blank->depth_of_coverage,      undef,'depth_of_coverage ok';
is $row_blank->depth_of_coverage_sd,   undef,'depth_of_coverage_sd ok';
is $row_blank->duplication_rate,       undef,'duplication_rate ok';
is $row_blank->error_rate,             undef,'error_rate ok';
is $row_blank->npg_qc,             'pending','npg_qc ok';
is $row_blank->manual_qc,           'passed','manual_qc ok';
is $row_blank->is_mapping_complete,    undef,'correctly marked as mapping not run';
is $row_blank->is_qc_mapstats,             0,'correctly marked not qc mapstats';

# valid qc lane
ok my $row_qc = Pathogens::Reports::Mapping::Row->new(vrtrack => $vrtrack, lane => $vrobject{lane}, mapstats => $vrobject{mapstats_qc}), 'Read lane qc mapstats to row.';

# qc lane
is $row_qc->study_id,                1000,'study_id ok';
is $row_qc->sample,              'Sparky','sample ok';
is $row_qc->lanename,              '1234_5#6','lane ok';
is $row_qc->cycles,                    76,'cycles ok';
is $row_qc->reads,               10066330,'reads ok';
is $row_qc->bases,              765041080,'bases ok';
is $row_qc->map_type,                'QC','map_type ok';
is $row_qc->reference, 'Sue_FMNH_PR_2081','reference ok';
is $row_qc->reference_size,      23298589,'reference_size ok';
is $row_qc->mapper,        'velocimapper','mapper ok';
is $row_qc->mapstats_id,                1,'mapstats_id ok';
is $row_qc->adapter_perc,           '0.0','adapter ok'; # need larger number to check calc.
is $row_qc->transposon_perc,        undef,'transposon ok'; # need large number for calc.
is $row_qc->mapped_perc,            '99.0','mapped ok';
is $row_qc->paired_perc,             96.9,'paired ok';
is $row_qc->mean_insert_size,       315.6,'mean_insert_size ok';
is $row_qc->genome_covered,         73.99,'genome_covered ok';
is $row_qc->genome_covered_1x,      undef,'genome_covered_1x   ok';
is $row_qc->genome_covered_5x,      undef,'genome_covered_5x   ok';
is $row_qc->genome_covered_10x,     undef,'genome_covered_10x  ok';
is $row_qc->genome_covered_50x,     undef,'genome_covered_50x  ok';
is $row_qc->genome_covered_100x,    undef,'genome_covered_100x ok';
is $row_qc->depth_of_coverage,      30.45,'depth_of_coverage ok';
is $row_qc->depth_of_coverage_sd,   37.59,'depth_of_coverage_sd ok';
is $row_qc->duplication_rate,      0.0003,'duplication_rate ok';
is $row_qc->error_rate,             0.003,'error_rate ok';
is $row_qc->npg_qc,             'pending','npg_qc ok';
is $row_qc->manual_qc,           'passed','manual_qc ok';
is $row_qc->is_mapping_complete,        1,'correctly marked as mapping run';
is $row_qc->is_qc_mapstats,             1,'correctly marked as qc mapstats';

# valid mapping lane
ok my $row_map = Pathogens::Reports::Mapping::Row->new(vrtrack => $vrtrack, lane => $vrobject{lane}, mapstats => $vrobject{mapstats_map}), 'Read lane qc mapstats to row.';

# mapping lane
is $row_map->study_id,                1000,'study_id ok';
is $row_map->sample,              'Sparky','sample ok';
is $row_map->lanename,          '1234_5#6','lane ok';
is $row_map->cycles,                    76,'cycles ok';
is $row_map->reads,               10066330,'reads ok';
is $row_map->bases,              765041080,'bases ok';
is $row_map->map_type,           'Mapping','map_type ok';
is $row_map->reference, 'Sue_FMNH_PR_2081','reference ok';
is $row_map->reference_size,      23298589,'reference_size ok';
is $row_map->mapper,             'dinomap','mapper ok';
is $row_map->mapstats_id,                3,'mapstats_id ok';
is $row_map->adapter_perc,           undef,'adapter ok';
is $row_map->transposon_perc,        undef,'transposon ok';
is $row_map->mapped_perc,             95.5,'mapped ok';
is $row_map->paired_perc,             92.1,'paired ok';
is $row_map->mean_insert_size,         319,'mean_insert_size ok';
is $row_map->genome_covered,         undef,'genome_covered ok';
is $row_map->genome_covered_1x,       88.2,'genome_covered_1x   ok';
is $row_map->genome_covered_5x,       70.4,'genome_covered_5x   ok';
is $row_map->genome_covered_10x,      62.2,'genome_covered_10x  ok';
is $row_map->genome_covered_50x,      25.1,'genome_covered_50x  ok';
is $row_map->genome_covered_100x,      2.2,'genome_covered_100x ok';
is $row_map->depth_of_coverage,      30.16,'depth_of_coverage ok';
is $row_map->depth_of_coverage_sd,   37.12,'depth_of_coverage_sd ok';
is $row_map->duplication_rate,       undef,'duplication_rate ok';
is $row_map->error_rate,             undef,'error_rate ok';
is $row_map->npg_qc,             'pending','npg_qc ok';
is $row_map->manual_qc,           'passed','manual_qc ok';

ok open(my $csv_fh, ">/dev/null"),'Opened filehandle to /dev/null';
ok my $spreadsheet = Pathogens::Reports::Mapping::Spreadsheet->new( filehandle => $csv_fh, rows => [$row_qc, $row_map]), 'spreadsheet ok';

ok $spreadsheet->_output_headers,'wrote header';
ok $spreadsheet->_output_rows,   'wrote rows';
ok $spreadsheet->output_csv,     'wrote spreadsheet';
close $csv_fh;

my $report_output;
ok open(my $report_fh,">",\$report_output), 'opened report_fh to scalar';
ok my $report = Pathogens::Reports::Mapping::Report->new( vrtrack => $vrtrack, filehandle => $report_fh, lanes => [$vrobject{lane}]),'opened report ok';
ok $report->output_csv,'wrote report';

my $expected_report_output = qq["Study ID",Sample,Lane,Cycles,"Yield (Reads)","Yield (Bases)","Type (QC/Mapping)",Reference,"Reference Size",Mapper,"Mapstats ID","Adapter (%)","Transposon (%)","Mapped (%)","Paired (%)","Mean Insert Size","Genome Covered (%)","Genome Covered (% >= 1X)","Genome Covered (% >= 5X)","Genome Covered (% >= 10X)","Genome Covered (% >= 50X)","Genome Covered (% >= 100X)","Depth of Coverage (X)","Depth of Coverage StdDev (X)","Duplication Rate","Error Rate","NPG QC","Manual QC"\r
1000,Sparky,1234_5#6,76,10066330,765041080,QC,Sue_FMNH_PR_2081,23298589,velocimapper,1,0.0,NA,99.0,96.9,315.6,73.99,NA,NA,NA,NA,NA,30.45,37.59,0.0003,0.003,pending,passed\r
1000,Sparky,1234_5#6,76,10066330,765041080,Mapping,Sue_FMNH_PR_2081,23298589,dinomap,3,0.0,NA,95.5,92.1,319,NA,88.2,70.4,62.2,25.1,2.2,30.16,37.12,NA,NA,pending,passed\r
1000,Sparky,1234_5#6,76,10066330,765041080,Mapping,Sue_FMNH_PR_2081,23298589,dinomap,5,NA,NA,0.0,0.0,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,pending,passed\r
];
is $report_output,$expected_report_output,'report gives expected output';


# cleanup
clean_test_db($vrtrack);

done_testing();
exit;

sub clean_test_db
{
    my ($vrtrack) = @_;

    my $dbh = $vrtrack->{_dbh};
    foreach ($dbh->tables()){
    next if /\`latest_/;
    next if /schema_version/;
    $dbh->do("TRUNCATE TABLE $_");
    }
}
