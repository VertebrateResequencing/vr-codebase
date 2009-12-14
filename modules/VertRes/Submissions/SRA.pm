=head1 NAME

VertRes::Submissions::SRA - SRA submission utility functions

=head1 SYNOPSIS

use VertRes::Submissions::SRA;

my $dbsnp_util = VertRes::Submissions::SRA->new();

=head1 DESCRIPTION

General utility functions for making dbSNP submission files - 1 method per dbSNP submission section

=head1 AUTHOR

Thomas Keane: thomas.keane@sanger.ac.uk

=cut

package VertRes::Submissions::SRA;

use strict;
use warnings;

use VRTrack::VRTrack;

use base qw(VertRes::Base);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Submissions::SRA->new(database=>'database_name',lanes=>'lanes file', holduntil=>'YYYY-MM-DD', project=>'project', contact_email=>'x@sanger.ac.uk', contact_name=>'Thomas',library_source=>'genomic', library_selection=>'random', machine_code=>'0', center=>'SC');
 Function: Create a new VertRes::Submissions::SRA object.
 Returns : VertRes::Submissions::SRA object
 Args    : 	database: VR Track db name
			lanes: file of lane names
			holdUntil: date at which the lanes become publicly visible (optional)
			project: project name to be used in the xml??
			contact_email: person to contact with submission status/problems
			contact_name: name of contact person
			library_selection: 'random' for whole genome shotgun
			machine_code: 0 for 'unspecified', 1 for 'Illumina Genome Analyzer', 2 for 'Illumina Genome Analyzer II'
			accessions: file lane/accession pairs (tab separated) - this is provided post submission to set the submitted status of the lanes post-submission
=cut

sub new
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
	bless($self,$class);
	
	if ( !$$self{'database'} ) { $self->throw("Expected a database name to be provided!");}
	if( $$self{'lanes'} )
	{
		if ( $$self{'lanes'} && ! -f $$self{'lanes'} ) {$self->throw("Cant find file of lane names: ".$$self{'lanes'})}
		if( !$$self{'contact_name'} || !$$self{'contact_email'} ){$self->thow("Insufficient contact information provided - contact_name and contact_email required");}
		$self->throw( "Library selection not specified" ) unless $$self{'library_selection'};
		$self->warn( "Library selection not set to RANDOM" ) unless $$self{ 'library_selection' } =~ /^RANDOM$/i;
		$self->throw( "Machine code not specified correctly (0,1,2)" ) unless $$self{'machine_code'} && $$self{'machine_code'} =~ /^[0-2]$/;
		$self->warn( "Library source not set to GENOMIC" ) unless $$self{'library_source'} =~ /GENOMIC/i;
		
		if( $$self{'holdUntil'} )
		{
			my $hold = $$self{'holduntil'};
			$self->throw('Invalid holduntil data provided (DD-MM-YYYY)') unless $hold =~ /^(\d\d)-(\d\d)-(\d\d\d\d)$/;
			my ($day, $month, $year) = split( /-/, $hold );
			$self->throw('Invalid holduntil data provided (DD-MM-YYYY)') unless $day > 0 && $day < 31 && $month > 0 && $month < 13 && $year > 2009;
		}
		
		#connect to the database
		print "Connecting to ".$$self{'database'}.".....\n";
		my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate($$self{'database'}, 'r');
		
		if( exists($$self{'lanes'}) )
		{
			_validateLanes($self, $vrtrack);
		}
		else
		{
			_gatherAllUnsubmittedLanes($self, $$self{'project'}, $vrtrack);
		}
		
		_gatherMetaInformation($self, $vrtrack);
	}
	elsif( $$self{ 'accessions' } )
	{
		#connect to the database
		print "Connecting to ".$$self{'database'}.".....\n";
		my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate($$self{'database'}, 'r');
		
		_readAccessionInformation($self, $vrtrack);
	}
	
    return $self;
}

=head2 printLanes

 Title   : printLanes
 Usage   : my $obj = VertRes::Submissions::SRA->printLanes();
 Function: Print the list of lanes to be submitted
 Returns : n/a
 Args    : n/a

=cut
sub printLanes
{
	my $self = shift;
	
	my @lanes = @{ $$self{'lanes'} };
	foreach( @lanes )
	{
		print $_."\n";
	}
}

=head2 writeXMLs

 Title   : writeXMLs
 Usage   : $obj->writeXMLs();
 Function: Goes through the studies corresponding to the lanes and writes the various XML submission files
 Returns : n/a
 Args    : n/a

=cut
sub writeXMLs
{
	my $self = shift;
	my %studyInfo = %{$$self{'studyinfo'}};
	my %subnames = %{$$self{'subnames'}};
	
	open( my $lfh, ">lane_submission_info.tab" ) or $self->throw("Cannot create lane_submission_info.tab");
	foreach my $study( keys %studyInfo )
	{
		_createStudyXML( $_, $studyInfo{ $_ }, $subnames{ $_ } );
		
		_createSampleXML( $$self{'sampleInfo'}, $subnames{ $_ } );
		
		_createExperimentXML( $$self{'experimentInfo'}, $_, $subnames{ $_ } );
		
		_createRunXML( $$self{'laneInfo'}, $subnames{ $_ } );
		
		_createSubXML( $$self{'laneInfo'}, $subnames{ $_ } );
		
		print $lfh
	}
	close( $lfh );
}

=head2 createStudyXML

 Title   : createStudyXML
 Usage   : $obj->createStudyXML();
 Function: Creates the study XML for a submission
 Returns : n/a
 Args    : the SRA study alias, the SRA project name, the center project name

=cut
sub _createStudyXML 
{
    my ( $self, $studyname, $projectName, $subname ) = @_;
	
    open (my $STUDY, ">$subname.study.xml") or $self->throw( "Can't open $subname.study.xml : $!" );
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

	printf $STUDY $studyxml, ($studyname, uc( $$self{'center'} ), $projectName);
    close ($STUDY);
}

=head2 createSampleXML

 Title   : createSampleXML
 Usage   : $obj->createSampleXML();
 Function: Creates the sample XML file for submission
 Returns : n/a
 Args    : the SRA study alias, the SRA project name, the center project name

=cut
sub _createSampleXML 
{
    my ($self, $sampleref,  $subname) = @_;
    open (my $SAMPLE, ">$subname.sample.xml") or $self->throw( "Can't open $subname.sample.xml : $!\n" );
    print $SAMPLE qq(<?xml version="1.0" encoding="UTF-8"?>\n);
    print $SAMPLE qq(<SAMPLE_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n);

    my $sampxml = <<'EOXML';
    <SAMPLE alias="%s">
          <SAMPLE_NAME>
               <COMMON_NAME>%s</COMMON_NAME>
          </SAMPLE_NAME>
     </SAMPLE>
EOXML

    foreach my $samplename ( sort keys %{$sampleref} )
	{
		my $species = $sampleref->{$samplename}{species};
		printf $SAMPLE $sampxml, ($samplename,$species);
    }
    
    print $SAMPLE "</SAMPLE_SET>\n";
	
    close ($SAMPLE);
}

=head2 createExperimentXML

 Title   : createExperimentXML
 Usage   : $obj->createExperimentXML();
 Function: Creates the experiment XML file for submission
 Returns : n/a
 Args    : the SRA study alias, the SRA project name, the center project name

=cut
sub _createExperimentXML 
{
    my ($self, $expref, $studyname, $subname) = @_;
    open (my $EXP, ">$subname.experiment.xml") or $self->throw( "Can't open $subname.experiment.xml : $!" );
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

    foreach my $exp( sort keys %{$expref} )
	{
		my $sample = $expref->{$exp}{'sample'};
		my $species = $expref->{$exp}{'species'};
		my $cycles = $expref->{$exp}{'cycles'};
		my $library = $expref->{$exp}{'library'};
		my $machine = $expref->{$exp}{'machine'};
		my $studylink = 'accession';
		
		printf $EXP $expxml, ($exp,$studylink,$studyname,$species,$library,$sample,$library,$$self{'library_source'}, $$self{'library_selection'}, $cycles + 1, $machine, $cycles);
    }
	
	print $EXP "</EXPERIMENT_SET>";
    close ($EXP);
}

sub _createRunXML 
{
	my ($self, $laneref, $subname) = @_;
	open (my $RUN, ">$subname.run.xml") or $self->throw( "Can't open $subname.run.xml : $!" );
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
    foreach my $lane(sort keys %{$laneref})
	{
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

sub _createSubXML
{
	my ($self, $laneref, $subname) = @_;
	open (my $SUB, ">$subname.submission.xml") or $self->throw( "Can't open $subname.submission.xml : $!" );

	# header
    print $SUB qq(<?xml version="1.0" encoding="UTF-8"?>\n);
    print $SUB sprintf (qq(<SUBMISSION center_name="%s" submission_date="%s" submission_id="%s" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n), $$self{'center'},_isoDate(time),$subname);
	
    # Contact info
    my $conxml = <<'EOXML';
  <CONTACTS>
    <CONTACT inform_on_error="%s" inform_on_status="%s" name="%s" />
  </CONTACTS>
EOXML

	printf $SUB $conxml, ($$self{'contact_email'}, $$self{'contact_email'}, $$self{'contact_name'});

    # reference the other XML files
    my $actxml =  <<'EOXML';
    <ACTION>
      <ADD schema="%s" source="%s" />
    </ACTION>
EOXML
    print $SUB "  <ACTIONS>\n";
	
	printf $SUB $actxml, "study", "$subname.study.xml";
	printf $SUB $actxml, "sample", "$subname.sample.xml";

    foreach $_ qw(experiment run){
	printf $SUB $actxml, $_, "$subname.$_.xml";
    }
	
	#if there's a hold date - then add it as a tag
	if( $$self{'holduntil'} )
	{
		my $date = $$self{'holduntil'};
		print $SUB qq[<ACTION>
						<HOLD HoldUntilDate="$date" 
					</ACTION>];
	}

    print $SUB "  </ACTIONS>\n";
	
    # Files
    my $filexml =  <<'EOXML';
    <FILE filename="%s" checksum="%s" checksum_method="MD5" />
EOXML

    print $SUB "  <FILES>\n";
    foreach my $lane(sort keys %{$laneref})
	{
		my $md5= $laneref->{$lane}{'md5'};
		printf $SUB $filexml, ($lane,$md5);
    }
    print $SUB "  </FILES>\n";

    print $SUB "</SUBMISSION>\n";
    close $SUB;
}

=head2 writeToDatabase

 Title   : writeToDatabase
 Usage   : $obj->writeToDatabase();
 Function: Writes the post submission accessions to the database (assuming this info was provided during the object construction)
 Returns : n/a
 Args    : n/a

=cut
sub writeToDatabase
{
	my $self = shift;
	my $vrtrack = shift;
	$self->throw("Lane accessions must be defined on object construction!") unless $$self{ 'laneAccessions' };
	
	my %accessions = %{ $$self{ 'laneAccessions' } };
	$vrtrack->transaction_start();
	foreach( keys( %accessions ) )
	{
		my $lane = VRTrack::Lane->new_by_name( $vrtrack, $_ );
		my $submission = $lane->add_submission($accessions{$_}[ 0 ]);
		$submission->acc( $accessions{$_}[ 1 ] );
	}
	$vrtrack->transaction_commit();
}

#private function to validate which of the lanes provided are eligible for submission
sub _validateLanes
{
	my $self = shift;
	my $vrtrack = shift;

	#read in the lane names
	my @laneNames;
	open( my $fh, $$self{'lanes'} ) or $self->throw('Cant open file of lane names: '.$$self{'lanes'});
	while( <$fh> )
	{
		chomp;
		
		#small sanity check
		$self->throw("Invalid lane name: ".$_) unless $_ !~ /^\d+_\d+$/;
		
		push( @laneNames, $_ );
	}
	close( $fh );
	
	my $originally = @laneNames;
	
	for(my $i = 0; $i < @laneNames; $i ++ )
	{
		my $lane = VRTrack::Lane->new_by_name( $vrtrack, $_ );
		
		if( $lane->submission_id() && $lane->submission_id() > 0 )
		{
			delete( $laneNames[ $i ] );
			$i --;
		}
		elsif( $lane->qc_status() ne 'passed' )
		{
			$self->warn("Found lane not passed QC: $_\n");
		}
	}
	
	print scalar( @laneNames ).' out of '.$originally.' lanes are eligible for submission\n';
	
	$$self{'lanes'} = \@laneNames;
}

sub _readAccessionInformation
{
	my $self = shift;
	my $vrtrack = shift;
	
	#read in the lane names
	my %laneAccessions;
	open( my $fh, $$self{'accessions'} ) or $self->throw('Cant open file of accessions: '.$$self{'lanes'});
	while( <$fh> )
	{
		chomp;
		
		if( ! $_ =~ /^(\d+_\d+)\t([ERR|SRR]\d+)$]/ ){$self->throw("Invalid accession information line: $_");}
		
		$laneAccessions{ $1 } = $2;
	}
	close( $fh );
	
	my $originally = keys( %laneAccessions );
	foreach( keys( %laneAccessions ) )
	{
		my $lane = VRTrack::Lane->new_by_name( $vrtrack, $_ );
		
		if( $lane->submission_id() && $lane->submission_id() > 0 )
		{
			print "Lane already submitted: $_\n";
			delete( $laneAccessions{ $_ } );
		}
		elsif( $lane->qc_status() ne 'passed' )
		{
			$self->warn("Found lane not passed QC: $_");
		}
	}
	
	print scalar( keys( %laneAccessions ) )." lanes out of ".$originally." remain after checking\n";
}

#private function for finding the names of all the lanes eligible for submission
sub _gatherAllUnsubmittedLanes
{
	my ($self, $project, $vrtrack) = @_;
	
	my @projects;
	
	if( defined( $project ) )
	{
		my $p = VRTrack::Project->new_by_name($vrtrack, $project);
		$self->throw("Cant find project: $project") unless defined( $p );
		push( @projects, $p );
	}
	else
	{
		@projects = @{ $vrtrack->projects() };
	}
	
	my @laneNames;
	foreach( @projects )
	{
		my $samples = $_->samples();
		foreach( @{$samples} )
		{
			my $sample = $_;
			my $libraries = $sample->libraries();
			
			foreach( @{$libraries} )
			{
				my $library = $_;
				my $lanes = $library->lanes();
				
				foreach( @{$lanes} )
				{
					my $lane = $_;
					if( ! $lane->submission_id() && $lane->qc_status() eq 'passed' )
					{
						push( @laneNames, $lane->name() );
					}
				}
			}
		}
	}
	
	$$self{'lanes'} = \@laneNames;
}

#private function for gathering all the meta info that is required for the XML files
#check the study names
#check the library names
sub _gatherMetaInformation
{
	my ($self, $vrtrack ) = @_;
	my %experimentInfo;
	my %laneInfo;
	my %sampleInfo;
	my %studyInfo;
	my %subnames;
	
	my @lanes = @{ $$self{ 'lanes' } };
	foreach( @lanes )
	{
		my $lane = VRTrack::Lane->new_by_name( $vrtrack, $_ );
		my $library = $lane->library();
		my $sample = $library->sample();
		my $project = $sample->project();
		my $study_acc = $project->acc() or $self->throw('No accession for project: $project');
		my $projectName = $project->acc;
		
		my $sample_name = $sample->individual->name;
        my $sample_acc = $sample->acc;
        my $species = $sample->individual->species->name;
		my $libraryName = $library->name();
		$libraryName =~ s/[_ ]/-/g;
		
		if( ! $studyInfo{ $study_acc } )
		{
			$studyInfo{ $study_acc } = $projectName;
			$subnames{ $study_acc } = lc(sprintf("%s-%s-%s", $projectName, lc($$self{'center'}) ,dateStr(time)));
		}
		
		my ($file, $md5) = get_lane_srf($lane->name);
		my $cycles = $lane->read_len;
		
		my $experiment_id = sprintf("%s-%s-%s-%s", $libraryName,$sample_name,$cycles,$$self{'machine_code'} );
		
		$sampleInfo{$study_acc}{$sample_name}{species} = $species;
		$laneInfo{$study_acc}{$file} = {'filetype'  => 'srf',
			'md5'       => $md5,
			'experiment'=> $experiment_id,
			'lane'      => $lane,
		};
		
		unless( defined $experimentInfo{$study_acc}{$experiment_id} )
		{
			$experimentInfo{$study_acc}{$experiment_id} = 
			{
				'sample'    => $sample_name,
				'species'   => $species,
				'library'   => $libraryName,
				'cycles'    => $cycles,
				'machine'   => $$$self{'machine_code'},
			};
		}
	}
	
	$$self{ 'experimentInfo' } = \%experimentInfo;
	$$self{ 'laneInfo' } = \%laneInfo;
	$$self{ 'sampleInfo' } = \%sampleInfo;
	$$self{ 'studyInfo' } = \%studyInfo;
	$$self{ 'subnames' } = \%subnames;
}

#find the srf file for a lane
sub _get_lane_srf 
{
    my $lane = shift;
    my $srf = "$lane.srf";
    my $md5 = `mpsa_download -m -f $srf`;
    unless ($md5)
	{
        $srf =~ s/_/_s_/;
        $md5 = `mpsa_download -m -f $srf`;
    }

    unless ($md5)
	{
        die "Can't get md5 for $lane\n";
    }
    ($md5) = split /\s+/,$md5;
    return ($srf,$md5);
}

sub _dateStr 
{
    my ($t) = @_;
    my ($sec,$min,$hour,$mday,$mon,$year) = localtime($t);
    my $date = sprintf("%4d%.2d%.2d", $year+1900,$mon+1,$mday);
    return $date;
}

sub _isoDate 
{
    my ($t) = @_;
    my ($sec,$min,$hour,$mday,$mon,$year) = localtime($t);
    my $date = sprintf("%4d-%2d-%2dT%2d:%2d:%2dZ",
                      $year+1900,$mon+1,$mday,$hour,$min,$sec);
    $date =~ tr/ /0/;
    return $date;
}
1;
