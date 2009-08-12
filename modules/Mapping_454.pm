package Mapping_454;
use strict;
use Carp;
use File::Basename;
use File::Spec;
use Cwd;

use G1KUtilities;

=pod
	A function that takes a mapping root directory and a list of projects
	Goes through the hierarchy and checks for unmapped 454 lanes, and starts the mapping
=cut
my $G1K = $ENV{ 'G1K' };
my $SSAHA2 = $G1K."/bin/ssaha2";
my $SAMTOOLS = $G1K."/bin/samtools";

my $MALE_REF_FAI = $G1K.'/ref/human_b36_male.fa.fai';
my $FEMALE_REF_FAI = $G1K.'/ref/human_b36_female.fa.fai';

my $LOCAL_CACHE = '/localcache';
my $REMOTE_REF_DIR= $G1K.'/ref';

sub map_454
{
	croak "Usage: map_454 mappingRootDir projects_list lsf_queue male_ref_fa female_ref_fa [index_file]" unless @_ != 5 || @_ != 6;
	
	my $mRoot = shift;
	my $projs = shift;
	my $lsf_queue = shift;
	
	my $male_ref_fa = shift;
	my $female_ref_fa = shift;
	
	my $indexF = '';
	if( @_ == 1 )
	{
		$indexF = shift;
		croak "Cant find index file: $indexF\n" unless -f $indexF;
	}
	
	croak "Cant find mapping directory: $mRoot\n" unless -d $mRoot;
	croak "Cant find index file: $indexF\n" unless -f $indexF;
	
	my %projects;
	open( P, $projs ) or die $!;
	while( <P> )
	{
		chomp;
		$_ =~ tr/ /_/;
		$projects{ $_} = 1;
	}
	close( P );
	
	chdir( $mRoot );
	
	foreach my $project (`ls`) 
	{
		chomp( $project );
		
		$project =~ tr/ /_/;
		print "Checking $project\n";
		if( -d $project && $projects{ $project } )
		{
			chdir( $project );
			print "In ".getcwd."\n";
			
			foreach my $individual (`ls`) 
			{
				chomp( $individual );
				if( -d $individual )
				{
					chdir( $individual );
					print "In ".getcwd."\n";
					
					my $gender = G1KUtilities::path2Gender('/lustre/sf4/1kgenomes/ref/genders.txt');
					
					foreach my $technology (`ls`)
					{
						chomp( $technology );
						if( -d $technology && $technology eq '454' )
						{
							chdir( $technology );
							print "In ".getcwd."\n";
							
							foreach my $library (`ls`) 
							{
								chomp( $library );
								if( -d $library )
								{
									chdir( $library );
									print "In ".getcwd."\n";
									
									my $libID = int( rand( 1000000 ) );
									
									my $libraryLanesCompleted = 1;
									
									foreach my $lane (`ls`) 
									{
										chomp( $lane );
										
										if( -d $lane )
										{
											chdir( $lane );
											print "In ".getcwd."\n";
											
											my $laneID = int( rand( 1000000 ) );
											
											my $paired = 0;
											my $lane_read0 = '';
											my $lane_read1 = '';
											my $lane_read2 = '';
										
											if( -f "meta.info" || -l "meta.info" )
											{
												$lane_read0=`grep read0 meta.info | head -1 | awk -F: '{print \$2}'`;
												$lane_read1=`grep read1 meta.info | head -1 | awk -F: '{print \$2}'`;
												$lane_read2=`grep read2 meta.info | head -1 | awk -F: '{print \$2}'`;
												
												chomp( $lane_read0 );
												chomp( $lane_read1 );
												chomp( $lane_read2 );
											}
											else
											{
												print "ERROR: Cant find meta.info file in ".getcwd()."\n";
												chdir( '..' );
												next;
											}
											
											#create a sym link to the ref fa
											if( $gender eq 'male' || $gender eq 'unknown' )
											{
												symLinkReferenceFasta( $male_ref_fa );
											}
											else
											{
												symLinkReferenceFasta( $female_ref_fa );
											}
											
											my $currentDir = getcwd();
											foreach my $fastq ($lane_read0, $lane_read1, $lane_read2)
											{
												if( -f $fastq.".touch" )
												{
													print "Fastq being mapped: $fastq\n";
												}
												elsif( $fastq )
												{
													my $fileID = int( rand( 1000000 ) );
													my @s = split( /\./, $fastq );
													my $cigarName = $s[ 0 ].'.'.$s[ 1 ].'.cigar.gz';
													if( -f $cigarName && `zcat $cigarName | tail -50 | grep "SSAHA2 finished" | wc -l` == 1 )
													{
														print "Fastq already mapped: $fastq\n";
													}
													else
													{
														system( "echo 'touched' > $fastq.touch" );
														
														my $preExec = qq/perl -w -e "use Mapping_454;Mapping_454::getLocalCopyRefSSAHA( \\"$gender\\" );"/;
														
														my $directory = int( rand( 1000000 ) );
														
														if( $gender eq 'male' || $gender eq 'unknown' )
														{
															my $cmd = qq{bsub -J map.$libID.$laneID.$fileID -R "select[type==X86_64 && mem > 6000] rusage[mem=6000:tmp=21000]" -M6000000 -q $lsf_queue -o $fastq.map.o -e $fastq.map.e -E '$preExec' "mkdir /tmp/$directory; gunzip -c $fastq > /tmp/$directory/reads.fastq; $SSAHA2 -454 -output cigar -diff 10 -save $LOCAL_CACHE/human_b36_male /tmp/$directory/reads.fastq | /nfs/team81/tk2/code/vert_reseq/user/tk2/miscScripts/filterCigarStreamTop10.pl | gzip -c > /tmp/$directory/$cigarName;cp /tmp/$directory/$cigarName $currentDir; rm -rf /tmp/$directory $fastq.touch"};
															#print $cmd."\n";
															system( $cmd );
														}
														elsif( $gender eq 'female' )
														{
															my $cmd = qq{bsub -J map.$libID.$laneID.$fileID -R "select[type==X86_64 && mem > 6000] rusage[mem=6000:tmp=21000]" -M6000000 -q $lsf_queue -o $fastq.map.o -e $fastq.map.e -E '$preExec' "mkdir /tmp/$directory; gunzip -c $fastq > /tmp/$directory/reads.fastq; $SSAHA2 -454 -output cigar -diff 10 -save $LOCAL_CACHE/human_b36_female.Y /tmp/$directory/reads.fastq | /nfs/team81/tk2/code/vert_reseq/user/tk2/miscScripts/filterCigarStreamTop10.pl | gzip -c > /tmp/$directory/$cigarName;cp /tmp/$directory/$cigarName $currentDir; rm -rf /tmp/$directory $fastq.touch"};
															#print $cmd."\n";
															system( $cmd );
														}
													}
												}
											}
											
											#setoff the BAM conversion after the mapping is done (if required)
											if( length( $indexF ) > 0 )
											{
												laneToBAM( getcwd(), $libID, $laneID, $lane_read0, $lane_read1, $lane_read2, basename( getcwd() ), $gender, $indexF, $lsf_queue );
											}
											
											chdir( ".." );
										}
									}
									chdir( ".." );
								}
							}
							chdir( ".." );
						}
					}
					chdir( ".." );
				}
			}
			chdir( ".." );
		}
	}
}

sub symLinkReferenceFasta
{
	croak "Usage: symLinkReference fasta" unless @_ == 1;
	
	my $fasta = shift;
	croak "Cant find reference fasta: $fasta\n" unless -f $fasta;
	
	if( ! -l 'ref.fa' )
	{
		symlink( $fasta, 'ref.fa' );
	}
}

sub getLocalCopyRefSSAHA
{
	croak "Usage: getLocalCopyRefSSAHA gender" unless @_ == 1;
	
	my $gender = shift;
	
	if( $gender eq 'male' || $gender eq 'unknown' )
	{
		checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_male.base", "$REMOTE_REF_DIR/human_b36_male.base.md5" );
		checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_male.body", "$REMOTE_REF_DIR/human_b36_male.body.md5" );
		checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_male.head", "$REMOTE_REF_DIR/human_b36_male.head.md5" );
		checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_male.name", "$REMOTE_REF_DIR/human_b36_male.name.md5" );
		checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_male.size", "$REMOTE_REF_DIR/human_b36_male.size.md5" );
	}
	elsif( $gender eq 'female' )
	{
		checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_female.Y.base", "$REMOTE_REF_DIR/human_b36_female.Y.base.md5" );
		checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_female.Y.body", "$REMOTE_REF_DIR/human_b36_female.Y.body.md5" );
		checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_female.Y.head", "$REMOTE_REF_DIR/human_b36_female.Y.head.md5" );
		checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_female.Y.name", "$REMOTE_REF_DIR/human_b36_female.Y.name.md5" );
		checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_female.Y.size", "$REMOTE_REF_DIR/human_b36_female.Y.size.md5" );
	}
}

sub checkDownloadLocalFile
{
	croak "Usage: downloadLocalFile remote_file remote_md5" unless @_ == 2;
	
	my $remoteF = shift;
	my $remote_md5 = shift;
	
	croak "Cant find remote file: $remoteF\n" unless -f $remoteF;
	croak "Cant find remote md5 file: $remote_md5\n" unless -f $remote_md5;
	
	my $remoteFilename = basename( $remoteF );
	if( ! -f $LOCAL_CACHE.'/'.$remoteFilename || `md5sum $LOCAL_CACHE/$remoteFilename | awk '{print \$1}'` eq `awk '{print \$1}' $remote_md5` )
	{
		#download
		system( "cp $remoteF $LOCAL_CACHE" );
		
		my $numRetries = 0;
		while( `md5sum $LOCAL_CACHE/$remoteFilename | awk '{print \$1}'` ne `awk '{print \$1}' $remote_md5` && $numRetries < 10 )
		{
			print "ERROR: Failed to download local copy of $remoteFilename Trying again....\n";
			system( "cp $remoteF $LOCAL_CACHE" );
			$numRetries ++;
		}
		
		croak "ERROR: Failed to download local copy of $remoteFilename after 10 retries\n" unless $numRetries < 10;
	}
}

#assume current directory is a lane
sub laneToBAM
{
	croak "Usage: laneToBAM cwd libID laneID read0 read1 read2 laneAccession gender index_file [lsf_queue]\n" unless @_ == 10;
	my $cwd = shift;
	my $libID = shift;
	my $laneID = shift;
	my $lane_read0 = shift;
	my $lane_read1 = shift;
	my $lane_read2 = shift;
	my $accession = shift;
	my $gender = shift;
	my $indexF = shift;
	my $lsf_queue = shift;
	$lsf_queue ||= 'normal';
	
	croak "Cant find lane directory: $cwd\n" unless -d $cwd;
	croak "Cant find index file: $indexF\n" unless -f $indexF;
	#unlink( "rmdup.bam" )  unless ! -f "rmdup.bam";
	chdir( $cwd );
	if( -f "rmdup.bam" && -s "rmdup.bam" > 2000 )
	{
		print "Skipping - found existing rmdup.bam file\n";
		return;
	}
	
	my $library = basename(dirname( $cwd )); #library is the directory name above lane
	chomp( $library );
	
	if( length( $library ) == 0 )
	{
		print "Cant find lane accession in index file - skipping BAM conversion\n";
		return;
	}
	
	my $individual = `grep $accession $indexF | head -1 | awk -F"\t" '{print \$10}'`;
	chomp( $individual );
	
	my $runName = `grep $accession $indexF | head -1 | awk -F"\t" '{print \$16}'`;
	chomp( $runName );
	
	my $platform = `grep $accession $indexF | head -1 | awk -F"\t" '{print \$13}'`;
	chomp( $platform );

	my $centre = `grep $accession $indexF | head -1 | awk -F"\t" '{print \$6}'`;
	chomp ($centre);

	if( length( $individual ) == 0 || length( $runName ) == 0 )
	{
		print "Cant get index file entries: $accession\n";
		return;
	}
	
	my $insertSize = 0;
	if( $lane_read1 && $lane_read2 )
	{
		$insertSize = `grep $accession $indexF | head -1 | awk -F"\t" '{print \$18}'`;
		chomp( $insertSize );
	}
	
	my $jobCondition = '';
	#check if there are pending mapping jobs for the lane
	if( `bjobs -w | grep "map.$libID.$laneID." | wc -l` > 0 )
	{
		$jobCondition = ' -w "done( map.$libID.$laneID.* )"';
	}
	
	system( "echo 'touched' > lane.touch" );
	
	if( $lane_read0 )
	{
		$lane_read0 =~ /^(.*)\.fastq\.gz$/;
		my $cigarName0 = $1.'.cigar.gz';
		
		open( S, ">unpaired.bam.sh" ) or die $!;
		print S qq/if [ -s $cigarName0 ] && [ `zcat $cigarName0 | tail -50 | grep "SSAHA2 finished" | wc -l` -eq 1 ]
		then
			perl -w -e \"use VertRes::Utils::Cigar; \$c = VertRes::Utils::Cigar->new(); \$c->cigar_to_sam(['$lane_read0'], ['$cigarName0'], undef, '$accession', 'unpaired.sam'); system('gzip unpaired.sam');\"
		fi/;
		close( S );
		
		#my $cmd = "bsub -J sam.$libID.$laneID.0 -q $lsf_queue -o bam.o -e bam.e $jobCondition perl -w -e \"use SamTools;SamTools::ssaha2samUnpaired( '$lane_read0', '$cigarName', '$accession', 'unpaired.sam.gz');\"";
		my $cmd = qq/bsub -J sam.$libID.$laneID.0 -q $lsf_queue -o bam.o -e bam.e $jobCondition "sh unpaired.bam.sh;rm unpaired.bam.sh"/;
		#print $cmd."\n";
		system( $cmd );
	}
	
	if( $lane_read1 && $lane_read2 )
	{
		$lane_read1 =~ /^(.*)\.fastq\.gz$/;
		my $cigarName1 = $1.'.cigar.gz';
		
		$lane_read2 =~ /^(.*)\.fastq\.gz$/;
		my $cigarName2 = $1.'.cigar.gz';
		
		open( S, ">paired.bam.sh" ) or die $!;
		print S qq/if [ -s $cigarName1 ] && [ -s $cigarName2 ] && [ `zcat $cigarName1 | tail -20 | grep "SSAHA2 finished" | wc -l` -eq 1 ] && [ `zcat $cigarName2 | tail -50 | grep "SSAHA2 finished" | wc -l` -eq 1 ]
		then
			perl -w -e \"use VertRes::Utils::Cigar; \$c = VertRes::Utils::Cigar->new(); \$c->cigar_to_sam(['$lane_read1', '$lane_read2'], ['$cigarName1', '$cigarName2'], '$insertSize', '$accession', 'paired.sam'); system('gzip paired.sam');\"
		fi/;
		close( S );
		
		#my $cmd = "bsub -J sam.$libID.$laneID.1 -q $lsf_queue -o bam.o -e bam.e $jobCondition perl -w -e \"use SamTools;SamTools::ssaha2samPaired( '$lane_read1', '$lane_read1.cigar.gz', '$lane_read2','$read2Cigar', '$insertSize', '$accession', 'paired.sam.gz');\"";
		my $cmd = qq/bsub -J sam.$libID.$laneID.1 -q $lsf_queue -o bam.o -e bam.e $jobCondition "sh paired.bam.sh;rm paired.bam.sh"/;
		#print $cmd."\n";
		system( $cmd );
	}
	
	#write the SAM header
	if( -f "raw.sam" )
	{
		unlink( "raw.sam" );
	}
	
	open( H, ">raw.sam" ) or die "Cannot create raw.sam: $!\n";
	print H "\@HD\tVN:1.0\tSO:coordinate\n";
	
# 2009-04-03 jws: removed this, as samtools adds its own set of SQ headers, so we had two.
#	#write the SQ headers
#	if( $gender eq 'male' || $gender eq 'unknown' )
#	{
#		open( FAI, $MALE_REF_FAI ) or die $!;
#		while( <FAI> )
#		{
#			chomp( $_ );
#			my @s = split( /\s+/ , $_ );
#			print H "\@SQ\tSN:$s[0]\tAS:NCBI36\tLN:$s[1]\n";
#		}
#		close( FAI );
#	}
#	else
#	{
#		open( FAI, $FEMALE_REF_FAI ) or die $!;
#		while( <FAI> )
#		{
#			chomp( $_ );
#			my @s = split( /\s+/ , $_ );
#			print H "\@SQ\tSN:$s[0]\tAS:NCBI36\tLN:$s[1]\n";
#		}
#		close( FAI );
#	}
	
	print H "\@RG\tID:$accession\tPU:$runName\tLB:$library\tSM:$individual";
	if( $insertSize =~ /\d+/ ){
	    print H "\tPI:$insertSize";
	}
	if( $centre ){
	    print H "\tCN:$centre";
	}
	if( $platform ){
	    print H "\tPL:$platform";
	}

	print H "\n";
	close( H );
	
	open( S, ">makeBam.sh" ) or die $!;
	print S qq/if [ $lane_read0 ] && [ ! -s unpaired.sam.gz ]
	then
		exit 1
	fi/;
	
	if( $lane_read1 && $lane_read2 )
	{
		print S qq/if [ $lane_read1 ] && [ $lane_read2 ] && [ ! -s paired.sam.gz ]
		then
			exit 1
		fi
		/;
	}
	
	print S "\n";
	print S qq/zcat *.sam.gz >> raw.sam;rm *.sam.gz;/;
	
	my $REF_FAI;
	if( $gender eq 'male' || $gender eq 'unknown' ){
	    $REF_FAI = $MALE_REF_FAI;
	}
	else {
	    $REF_FAI = $FEMALE_REF_FAI;
	}

	print S	qq/$SAMTOOLS import $REF_FAI raw.sam - | $SAMTOOLS sort - raw.sorted; $SAMTOOLS rmdup raw.sorted.bam rmdup.bam; $SAMTOOLS flagstat rmdup.bam > rmdup.bam.flagstat ; $SAMTOOLS index rmdup.bam; rm raw.sorted.*;rm raw.sam/;
	close( S );
	
	my $cmd = qq/bsub -J bam.$libID.$laneID -q $lsf_queue -o bam.o -e bam.e -w "done(sam.$libID.$laneID.*)" "sh makeBam.sh;rm makeBam.sh;rm lane.touch"/;
	#print $cmd."\n";
	system( $cmd );
}

1;
