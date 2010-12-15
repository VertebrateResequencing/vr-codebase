#!/usr/bin/env perl

use strict;
use warnings;
no warnings 'uninitialized';
use Getopt::Long;

my ($lanefile, $nosample, $nostudy, $libsource, $libselection, $help);

GetOptions(
    'l|lanes=s'     =>  \$lanefile,
    'nosample'	    =>  \$nosample,
    'nostudy'	    =>  \$nostudy,
    'source=s'	    =>  \$libsource,
    'selection=s'   =>  \$libselection,
    'h|help'        =>  \$help,
    );

$libsource ||= 'GENOMIC';
$libselection ||= 'RANDOM';

(-s $lanefile && !$help) or die <<USAGE;
    Usage: $0
                --lanes	    <tab-delimited file of lane info as described below>
		[--source    <library source, default: GENOMIC>]
		[--selection <library selection method, default: RANDOM>]
		[--nosample  <don't make sample XML - for existing projects>]
		[--nostudy   <don't make study XML - for existing projects>]
                --help	    <this message>

Generates a minimal set of submission files for the ERA/SRA.

Takes a tab-delimited file of lane information, one line/lane, with fields:

study (either an alias or an SRP/ERP accession for existing studies)
project (internal name for study, e.g. 'g1k','gorilla')
filename (e.g. 1234_s_1.srf)
filetype (e.g. srf)
md5 of srf
sample (e.g. ggo-1)
library (e.g. GGO-200bp-1)
species (e.g. Gorilla gorilla gorilla)
cycles (e.g. 37)
machine type (i.e. Illumina Genome Analyzer or Illumina Genome Analyzer II or unspecified)

Builds one set of submission files for each study in the file.

Note that this is extremely simple, and assumes paired-end Illumina
sequencing.  It should produce the basic & minimal study, sample, experiment,
run, and submission XMLs required to submit to the ERA.

USAGE

$|++;

our $CENTER = 'SC';
our $center = 'sc';
our $DOMAIN      = 'sanger.ac.uk';
our %MACHINECODE = ('Illumina Genome Analyzer'	    => '1',
		    'Illumina Genome Analyzer II'   => '2',
	            'unspecified'		    => '0',
		  );


open (my $LANE, $lanefile) or die "Can't open lane file: $!\n";

my %lane;
my %experiment;
my %sample;
my %study;
while (<$LANE>){
    chomp;
    my ($studyname, $projname, $file, $ftype, $md5, $sample, $library, $species, $cycles, $machine) = split "\t",$_; 
    unless ($studyname && $projname && $file && $ftype && $md5 && $sample && $library && $species && $cycles && $machine){
        die "Missing data in line:\n$_\n";
    }
    die "Duplicate: $file" if exists $lane{$file};
    $projname =~ s/ /_/g;
    $study{$studyname} = $projname;
    my $slxcode = $MACHINECODE{$machine};
    die "Unrecognised platform: $machine" unless defined $slxcode;

    # need library to match existing G1K libraries via David Carter's script,
    # e.g g1k-sc-NA18980-JPT-1
    $library =~ s/[_ ]/-/g;
    $library = sprintf("%s-%s-%s",$projname,$center,$library);

    my $experiment_id = sprintf("%s-%s-%s-%s", $library,$sample,$cycles,$slxcode);
    $lane{$studyname}{$file} = {'filetype'  => $ftype,
                    'md5'       => $md5,
		    'experiment'=> $experiment_id,
                    };

    $sample{$studyname}{$sample}{species} = $species;

    unless (defined $experiment{$studyname}{$experiment_id}){
	$experiment{$studyname}{$experiment_id} = {	
					'sample'    => $sample,
					'species'   => $species,
					'library'   => $library,
					'cycles'    => $cycles,
					'machine'   => $machine,
				      };
    }
}


foreach my $study(keys %study){
    # Submission name:
    my $projname = $study{$study};
    my $subname = lc(sprintf("%s-%s-%s", $projname, $center ,dateStr(time)));

    createStudyXML($study, $projname, $subname) unless $nostudy;
    createSampleXML($sample{$study}, $subname) unless $nosample;
    createExperimentXML($experiment{$study}, $study, $subname, $libsource, $libselection);
    createRunXML($lane{$study}, $subname);
    createSubXML($lane{$study}, $subname);
}


sub createStudyXML {
    my ($studyname, $projname, $subname) = @_;
    open (my $STUDY, ">$subname.study.xml") or die "Can't open $subname.study.xml : $!\n";
    my $studyxml = <<'EOXML';
<?xml version="1.0" encoding="UTF-8"?>
<STUDY_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
     <STUDY alias="%s">
          <DESCRIPTOR>
               <STUDY_TITLE></STUDY_TITLE>
               <STUDY_TYPE existing_study_type="Whole Genome Sequencing"/>
               <CENTER_NAME>%s</CENTER_NAME>
               <CENTER_PROJECT_NAME>%s</CENTER_PROJECT_NAME>
               <PROJECT_ID>0</PROJECT_ID>
          </DESCRIPTOR>
     </STUDY>
</STUDY_SET>
EOXML

    printf $STUDY $studyxml, ($studyname, $CENTER,$projname);
    close ($STUDY);
}

sub createSampleXML {
    my ($sampleref,  $subname) = @_;
    open (my $SAMPLE, ">$subname.sample.xml") or die "Can't open $subname.sample.xml : $!\n";
    print $SAMPLE qq(<?xml version="1.0" encoding="UTF-8"?>\n);
    print $SAMPLE qq(<SAMPLE_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n);

    my $sampxml = <<'EOXML';
    <SAMPLE alias="%s">
          <SAMPLE_NAME>
               <COMMON_NAME>%s</COMMON_NAME>
          </SAMPLE_NAME>
     </SAMPLE>
EOXML

    foreach my $samplename (sort keys %{$sampleref}){
	my $species = $sampleref->{$samplename}{species};
	printf $SAMPLE $sampxml, ($samplename,$species);
    }
    
    print $SAMPLE "</SAMPLE_SET>\n";

    close ($SAMPLE);
}

sub createExperimentXML {
    my ($expref, $studyname, $subname, $libsource, $libselection) = @_;
    open (my $EXP, ">$subname.experiment.xml") or die "Can't open $subname.experiment.xml : $!\n";
    print $EXP qq(<?xml version="1.0" encoding="UTF-8"?>\n);
    print $EXP qq(<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n);
    my $expxml = <<'EOXML';
     <EXPERIMENT alias="%s">
          <STUDY_REF %s="%s"/>
          <DESIGN>
               <DESIGN_DESCRIPTION>Illumina sequencing of %s %s library</DESIGN_DESCRIPTION>
               <SAMPLE_DESCRIPTOR refname="%s"/>
               <LIBRARY_DESCRIPTOR>
                    <LIBRARY_NAME>%s</LIBRARY_NAME>
                    <LIBRARY_SOURCE>%s</LIBRARY_SOURCE>
		    <LIBRARY_SELECTION>%s</LIBRARY_SELECTION>
                    <LIBRARY_LAYOUT>
                         <PAIRED />
                    </LIBRARY_LAYOUT>
               </LIBRARY_DESCRIPTOR>
               <SPOT_DESCRIPTOR>
                    <SPOT_DECODE_SPEC>
			 <NUMBER_OF_READS_PER_SPOT>2</NUMBER_OF_READS_PER_SPOT>
			 <READ_SPEC>
                              <READ_INDEX>0</READ_INDEX>
                              <READ_CLASS>Application Read</READ_CLASS>
                              <READ_TYPE>Forward</READ_TYPE>
                              <BASE_COORD>1</BASE_COORD>
                         </READ_SPEC>
                         <READ_SPEC>
                              <READ_INDEX>1</READ_INDEX>
                              <READ_CLASS>Application Read</READ_CLASS>
                              <READ_TYPE>Reverse</READ_TYPE>
                              <BASE_COORD>%s</BASE_COORD>
                         </READ_SPEC>
                    </SPOT_DECODE_SPEC>
               </SPOT_DESCRIPTOR>
          </DESIGN>
          <PLATFORM>
               <ILLUMINA>
                    <INSTRUMENT_MODEL>%s</INSTRUMENT_MODEL>
                    <CYCLE_SEQUENCE/>
                    <CYCLE_COUNT>%s</CYCLE_COUNT>
               </ILLUMINA>
          </PLATFORM>
          <PROCESSING>
               <BASE_CALLS>
                    <SEQUENCE_SPACE>Base Space</SEQUENCE_SPACE>
                    <BASE_CALLER>Solexa primary analysis</BASE_CALLER>
               </BASE_CALLS>
               <QUALITY_SCORES qtype="other">
                    <QUALITY_SCORER>Solexa primary analysis</QUALITY_SCORER>
                    <NUMBER_OF_LEVELS>80</NUMBER_OF_LEVELS>
                    <MULTIPLIER>1</MULTIPLIER>
               </QUALITY_SCORES>
          </PROCESSING>
     </EXPERIMENT>
EOXML

    foreach my $exp(sort keys %{$expref}){
	my $sample = $expref->{$exp}{'sample'};
	my $species = $expref->{$exp}{'species'};
	my $cycles = $expref->{$exp}{'cycles'};
	my $library = $expref->{$exp}{'library'};
	my $machine = $expref->{$exp}{'machine'};
	my $studylink;	# link by name or by accession?
	if ($studyname =~ /^[ES]RP\d+$/){
	    $studylink = 'accession';
	}
	else {
	    $studylink = 'refname';
	}

	printf $EXP $expxml, ($exp,$studylink,$studyname,$species,$library,$sample,$library,$libsource, $libselection, $cycles + 1, $machine, $cycles);
    }

    print $EXP "</EXPERIMENT_SET>";
    close ($EXP);
}


sub createRunXML {
    my ($laneref, $subname) = @_;
    open (my $RUN, ">$subname.run.xml") or die "Can't open $subname.run.xml : $!\n";
    print $RUN qq(<?xml version="1.0" encoding="UTF-8"?>\n);
    print $RUN qq(<RUN_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n);
    my $runxml = <<'EOXML';
     <RUN alias="%s" total_data_blocks="1">
          <EXPERIMENT_REF refname="%s"/>
          <DATA_BLOCK name="%s">
               <FILES>
                    <FILE filename="%s" filetype="%s"/>
               </FILES>
          </DATA_BLOCK>
     </RUN>
EOXML

    my %exp_count;
    foreach my $lane(sort keys %{$laneref}){
	my $filetype= $laneref->{$lane}{'filetype'};
	my $blockname = $lane;
	$blockname =~ s/\..*//;	# trim off type
	my $exp  = $laneref->{$lane}{'experiment'};
	$exp_count{$exp}++; 
	my $index = $exp_count{$exp}; 
	my $runname = "${exp}_$index";

	printf $RUN $runxml, ($runname,$exp,$blockname,$lane,$filetype);
    }
    print $RUN "</RUN_SET>";
    close ($RUN);
}


sub createSubXML {
    my ($laneref, $subname) = @_;
    open (my $SUB, ">$subname.submission.xml") or die "Can't open $subname.submission.xml : $!\n";

    # header
    print $SUB qq(<?xml version="1.0" encoding="UTF-8"?>\n);
    print $SUB sprintf (qq(<SUBMISSION center_name="%s" submission_date="%s" submission_id="%s" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n), $CENTER,isoDate(time),$subname);

    # Contact info
    my $conxml = <<'EOXML';
  <CONTACTS>
    <CONTACT inform_on_error="%s" inform_on_status="%s" name="%s" />
  </CONTACTS>
EOXML

    my $userdata = getUserData();
    my $email = $userdata->{mail};
    my $name = $userdata->{name};
    printf $SUB $conxml, ($email, $email, $name);

    # reference the other XML files
    my $actxml =  <<'EOXML';
    <ACTION>
      <ADD schema="%s" source="%s" />
    </ACTION>
EOXML
    print $SUB "  <ACTIONS>\n";
    unless($nostudy){
	printf $SUB $actxml, "study", "$subname.study.xml";
    }
    unless ($nosample){
	printf $SUB $actxml, "sample", "$subname.sample.xml";
    }
    foreach $_ qw(experiment run){
	printf $SUB $actxml, $_, "$subname.$_.xml";
    }

    print $SUB "  </ACTIONS>\n";

    # Files
    my $filexml =  <<'EOXML';
    <FILE filename="%s" checksum="%s" checksum_method="MD5" />
EOXML

    print $SUB "  <FILES>\n";
    foreach my $lane(sort keys %{$laneref}){
	my $md5= $laneref->{$lane}{'md5'};
	printf $SUB $filexml, ($lane,$md5);
    }
    print $SUB "  </FILES>\n";

    print $SUB "</SUBMISSION>\n";
    close $SUB;
}


sub dateStr {
    my ($t) = @_;
    my ($sec,$min,$hour,$mday,$mon,$year) = localtime($t);
    my $date = sprintf("%4d%.2d%.2d", $year+1900,$mon+1,$mday);
    return $date;
}


sub isoDate {
    my ($t) = @_;
    my ($sec,$min,$hour,$mday,$mon,$year) = localtime($t);
    my $date = sprintf("%4d-%2d-%2dT%2d:%2d:%2dZ",
                      $year+1900,$mon+1,$mday,$hour,$min,$sec);
    $date =~ tr/ /0/;
    return $date;
}


sub getUserData {
    my $u = $ENV{USER};
    my @a = getpwnam($u);
    my $name = $a[6];
    $name =~ s/,.*//;
    return {uid => $u,
            name => $name,
            mail => "$u\@$DOMAIN"};
}

