package Mapping_454_ssaha;
use strict;
use Carp;
use File::Basename;
use File::Spec;
use Cwd;

use G1KUtilities;
use Mapping;

my $G1K = '/nfs/sf8/G1K'; #$ENV{ 'G1K' };
my $MOUSE_REF_FA = '/nfs/sf7/MOUSE/ref/NCBIM37_um.fa';
my $HUMAN_FEMALE_REF_FA = $G1K.'/ref/human_b36_female.fa';
my $HUMAN_MALE_REF_FA = $G1K.'/ref/human_b36_male.fa';

my $HUMAN_FEMALE_REF_FAI = $G1K.'/ref/human_b36_female.fa.fai';
my $HUMAN_MALE_REF_FAI = $G1K.'/ref/human_b36_male.fa.fai';

my $REMOTE_REF_DIR= $G1K.'/ref';

my $SSAHA2 = $G1K."/bin/ssaha2";
my $SAMTOOLS = $G1K."/bin/samtools";

my $filterCigarStreamTop10 = 'perl /nfs/users/nfs_t/tk2/code/tk2/miscScripts/filterCigarStreamTop10.pl';

=pod

=head1 NAME

Mapping_454_ssaha

=head1 SYNOPSIS

=head1 REQUIRES

Perl5.8.8

=head1 DESCRIPTION

=head1 METHODS

=head2 mappingInProgress

	Arg [1]    : lane directory
	Returntype : 0/1
=cut
sub mappingInProgress
{
	croak "Usage: mappingInProgress lane_dir" unless @_ == 1;
	
	my $laneDir = shift;
	croak "Cant find lane directory: $laneDir\n" unless -d $laneDir;
	
	chdir( $laneDir );
	
	my $paired = 0;
	my $lane_read0 = '';
	my $lane_read1 = '';
	my $lane_read2 = '';
	my $expectedInsert = '';
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
	
	#check if the ssaha jobs completed
	foreach my $fastq ($lane_read0, $lane_read1, $lane_read2)
	{
		if( $fastq && -f $fastq.'.touch' )
		{
			return 1;
		}
	}
	return 0;
}

=head2 isLaneMapped

	Arg [1]    : lane directory
	Returntype : 0/1
=cut
sub isLaneMapped
{
	croak "Usage: isLaneMapped lane_dir" unless @_ == 1;
	
	my $laneDir = shift;
	croak "Cant find lane directory: $laneDir\n" unless -d $laneDir;
	
	chdir( $laneDir );
	
	my $paired = 0;
	my $lane_read0 = '';
	my $lane_read1 = '';
	my $lane_read2 = '';
	my $expectedInsert = '';
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
	
	#check if the ssaha jobs completed
	foreach my $fastq ($lane_read0, $lane_read1, $lane_read2)
	{
		if( $fastq )
		{
			my @s = split( /\./, $fastq );
			my $cigarName = $s[ 0 ].'.'.$s[ 1 ].'.cigar.gz';
			if( ! -f $cigarName )
			{
				return 0;
			}
			elsif( -f $cigarName.".mapstat" || -l $cigarName.".mapstat" )
			{
				my $hasEntry = `grep "num_reads_mapped" $cigarName.mapstat | wc -l`;chomp( $hasEntry );
				return 0 unless $hasEntry == 1;
			}
		}
	}
	
	return 1;
}

=head2 mapLane

	Arg [1]    : lane directory
	Arg [2]    : lane directory
	Arg [3]    : index file
	Returntype : none
=cut
sub mapLane
{
	croak "Usage: mapLane lane_dir lsf_queue index_file" unless @_ == 3;
	
	my $laneDir = shift;
	my $lsf_queue = shift;
	my $indexF = shift;
	
	croak "Cant find lane directory: $laneDir\n" unless -d $laneDir;
	
	chdir( $laneDir );
	
	croak "Cant find meta.info file: $laneDir\n" unless -f 'meta.info';
	
	#if the directory is being operated on at the moment - skip it
	if( -f "lane.touch" )
	{
		print "Lane is currently being mapped - skipping\n";
		return;
	}
	
	my $paired = 0;
	my $lane_read0 = '';
	my $lane_read1 = '';
	my $lane_read2 = '';
	my $expectedInsert = '';
	if( -f "meta.info" || -l "meta.info" )
	{
		$lane_read0=`grep read0 meta.info | head -1 | awk -F: '{print \$2}'`;
		$lane_read1=`grep read1 meta.info | head -1 | awk -F: '{print \$2}'`;
		$lane_read2=`grep read2 meta.info | head -1 | awk -F: '{print \$2}'`;
		$expectedInsert=`grep insert meta.info | head -1 | awk -F: '{print \$2}'`;
		chomp( $lane_read0 );
		chomp( $lane_read1 );
		chomp( $lane_read2 );
		chomp( $expectedInsert );
	}
	else
	{
		print "ERROR: Cant find meta.info file in ".getcwd()."\n";
		chdir( '..' );
		next;
	}
	
	my $laneID = basename( $laneDir );
	
	#if the directory is being operated on at the moment - skip it
	foreach my $fastq ($lane_read0, $lane_read1, $lane_read2)
	{
		if( $fastq )
		{
			if( -f "$fastq.touch" )
			{
				print "Lane is currently being mapped - skipping\n";
				return;
			}
		}
	}
	
	#check if the ssaha jobs completed
	my $completed = 1;
	foreach my $fastq ($lane_read0, $lane_read1, $lane_read2)
	{
		if( $fastq )
		{
			my @s = split( /\./, $fastq );
			my $cigarName = $s[ 0 ].'.'.$s[ 1 ].'.cigar.gz';
			if( ! -f $cigarName )
			{
				$completed = 0;
				last;
			}

			if( -f $cigarName )
			{
				my $done = `zcat $cigarName | tail -50 | grep "SSAHA2 finished" | wc -l`;
				chomp( $done );
				if( $done == 0 )
				{
					$completed = 0;
					last;
				}
			}
		}
	}
	
	if( $completed == 1 )
	{
		print "Fastq files are already mapped - making BAM file now\n";
		
		#only do the bam conversion when there is an index file (i.e. doesnt apply to mouse)
		if( -f $indexF )
		{
			#laneToBAM( getcwd(), $laneID, $lane_read0, $lane_read1, $lane_read2, basename( getcwd() ), $indexF, $lsf_queue );
		}
		return;
	}
	
	#create a sym link to the ref bfa
	symLinkReference();
	
	my $gender = G1KUtilities::path2Gender();
	
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
			
			#my $jobReqs = '-R "select[type==X86_64 && mem > 6500] rusage[mem=6500:tmp=12000]" -M6500000';
			my $jobReqs = '-R "select[type==X86_64] rusage[tmp=35000]"';
			
			if( -f $cigarName )
			{
				my $done = `zcat $cigarName | tail -50 | grep "SSAHA2 finished" | wc -l`;
				chomp( $done );
				if( $done >= 1 )
				{
					print "Fastq already mapped: $fastq\n";
					next;
				}
				else
				{
					#$jobReqs = '-R "select[type==X86_64 && mem > 10000] rusage[mem=10000:tmp=12000]" -M10000000';
					$jobReqs = '-R "select[type==X86_64] rusage[tmp=35000]"';
				}
			}
			
			system( "echo 'touched' > $fastq.touch" );
													
			my $preExec = qq/perl -w -e "use Mapping_454_ssaha;Mapping_454_ssaha::getLocalCopyReference( \\"$gender\\" );"/;
			
			my $directory = int( rand( 1000000 ) );
			my $currentDir = getcwd();
			
			unlink( "$fastq.map.o" ) unless ! -f "$fastq.map.o";
			unlink( "$fastq.map.e" ) unless ! -f "$fastq.map.e";
			
			if( $gender eq 'male' || $gender eq 'unknown' )
			{
				#my $cmd = qq{bsub -J map.454.$laneID.$fileID $jobReqs -q normal -o $fastq.map.o -e $fastq.map.e -E '$preExec' 'mkdir /tmp/$directory; gunzip -c $currentDir/$fastq > /tmp/$directory/reads.fastq; perl -w -e "use AssemblyTools;AssemblyTools::filterOutShortReads( \\"/tmp/$directory/reads.fastq\\", \\"30\\", \\"/tmp/$directory/t.fastq\\" );"; mv /tmp/$directory/t.fastq /tmp/$directory/reads.fastq; $SSAHA2 -454 -output cigar -diff 10 -save $Mapping::LOCAL_CACHE_DIR/human_b36_male /tmp/$directory/reads.fastq | /nfs/team81/tk2/code/vert_reseq/user/tk2/miscScripts/filterCigarStreamTop10.pl | gzip -c > /tmp/$directory/$cigarName;cp /tmp/$directory/$cigarName $currentDir; rm -rf /tmp/$directory $fastq.touch};
				
				my $cmd = qq{bsub -J map.454.$laneID.$fileID $jobReqs -q normal -o $fastq.map.o -e $fastq.map.e -E '$preExec' 'mkdir /tmp/$directory; gunzip -c $currentDir/$fastq > /tmp/$directory/reads.fastq; perl -w -e "use AssemblyTools;AssemblyTools::filterOutShortReads( \\"/tmp/$directory/reads.fastq\\", \\"30\\", \\"/tmp/$directory/t.fastq\\" );"; mv /tmp/$directory/t.fastq /tmp/$directory/reads.fastq; $SSAHA2 -disk 1 -454 -output cigar -diff 10 -save $Mapping::LOCAL_CACHE_DIR/human_b36_male /tmp/$directory/reads.fastq | $filterCigarStreamTop10 | gzip -c > /tmp/$directory/$cigarName;cp /tmp/$directory/$cigarName $currentDir; rm -rf /tmp/$directory $fastq.touch};
				
				$cmd .= qq{; perl -w -e "use Mapping_454_ssaha;Mapping_454_ssaha::cigarStat( \\"$currentDir/$cigarName\\", \\"$currentDir/$cigarName.mapstat\\");"'};
				#print $cmd."\n";
				system( $cmd );
			}
			elsif( $gender eq 'female' )
			{
				#my $cmd = qq{bsub -J map.454.$laneID.$fileID $jobReqs -q normal -o $fastq.map.o -e $fastq.map.e -E '$preExec' 'mkdir /tmp/$directory; gunzip -c $currentDir/$fastq > /tmp/$directory/reads.fastq; perl -w -e "use AssemblyTools;AssemblyTools::filterOutShortReads( \\"/tmp/$directory/reads.fastq\\", \\"30\\", \\"/tmp/$directory/t.fastq\\" );"; mv /tmp/$directory/t.fastq /tmp/$directory/reads.fastq; $SSAHA2 -454 -output cigar -diff 10 -save $Mapping::LOCAL_CACHE_DIR/human_b36_female.Y /tmp/$directory/reads.fastq | /nfs/team81/tk2/code/vert_reseq/user/tk2/miscScripts/filterCigarStreamTop10.pl | gzip -c > /tmp/$directory/$cigarName;cp /tmp/$directory/$cigarName $currentDir; rm -rf /tmp/$directory $fastq.touch};
				
				my $cmd = qq{bsub -J map.454.$laneID.$fileID $jobReqs -q normal -o $fastq.map.o -e $fastq.map.e -E '$preExec' 'mkdir /tmp/$directory; gunzip -c $currentDir/$fastq > /tmp/$directory/reads.fastq; perl -w -e "use AssemblyTools;AssemblyTools::filterOutShortReads( \\"/tmp/$directory/reads.fastq\\", \\"30\\", \\"/tmp/$directory/t.fastq\\" );"; mv /tmp/$directory/t.fastq /tmp/$directory/reads.fastq; $SSAHA2 -disk 1 -454 -output cigar -diff 10 -save $Mapping::LOCAL_CACHE_DIR/human_b36_female.Y /tmp/$directory/reads.fastq | $filterCigarStreamTop10 | gzip -c > /tmp/$directory/$cigarName;cp /tmp/$directory/$cigarName $currentDir; rm -rf /tmp/$directory $fastq.touch};
				
				$cmd .= qq{; perl -w -e "use Mapping_454_ssaha;Mapping_454_ssaha::cigarStat( \\"$currentDir/$cigarName\\", \\"$currentDir/$cigarName.mapstat\\");"'};
				#print $cmd."\n";
				system( $cmd );
			}
		}
	}
	
	#setoff the BAM conversion after the mapping is done (if required)
	if( length( $indexF ) > 0 )
	{
		#laneToBAM( getcwd(), $laneID, $lane_read0, $lane_read1, $lane_read2, basename( getcwd() ), $indexF, $lsf_queue );
	}
}

#assume current directory is a lane
sub laneToBAM
{
	croak "Usage: laneToBAM lane_dir lsf_queue index_file\n" unless @_ == 3;
	
	my $lane_dir = shift;
	my $lsf_queue = shift;
	my $indexF = shift;
	
	croak "Cant find index file: $indexF\n" unless -f $indexF;
	croak "Cant find lane directory: $lane_dir\n" unless -d $lane_dir;
	
	chdir( $lane_dir );
	
	croak "Cant find meta.info file: $lane_dir\n" unless -f 'meta.info';
	
	if( ! isLaneMapped( $lane_dir ) )
	{
		croak "Lane is not mapped yet!\n";
	}
	
	if( -f "lane.touch" )
	{
		print "BAM file is currently being made - skipping";
		return;
	}
	
	my $accession = basename( getcwd() );
	
	my @fastq = @{ HierarchyUtilities::getFastqInfo( getcwd() ) };
	my $lane_read0 = $fastq[ 0 ][ 0 ];
	my $lane_read1 = $fastq[ 1 ][ 0 ];
	my $lane_read2 = $fastq[ 2 ][ 0 ];
	
	#unlink( "rmdup.bam" )  unless ! -f "rmdup.bam";
	chdir( $lane_dir );
	if( -f "rmdup.bam" && -s "rmdup.bam" > 2000 )
	{
		print "Skipping - found existing rmdup.bam file\n";
		return;
	}
	
	my $library = basename(dirname( $lane_dir )); #library is the directory name above lane
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
	
	my $centre = `grep $accession $indexF | head -1 | awk -F"\t" '{print \$6}'`;
	chomp( $centre );
	
	my $platform = `grep $accession $indexF | head -1 | awk -F"\t" '{print \$13}'`;
	chomp( $platform );
	
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
	my $numJobs = `bjobs -w | grep "map.454.$accession." | wc -l`;
	if( $numJobs > 0 )
	{
		$jobCondition = ' -w "done( map.454.$accession.* )"';
	}
	
	system( "echo 'touched' > lane.touch" );
	
	unlink( "bam.o" ) unless ! -f 'bam.o';
	unlink( "bam.e" ) unless ! -f 'bam.e';
	
	if( -f $lane_read0 || -l $lane_read0 )
	{
		$lane_read0 =~ /^(.*)\.fastq\.gz$/;
		my $cigarName0 = $1.'.cigar.gz';
		
		open( S, ">unpaired.bam.sh" ) or die $!;
		print S qq/if [ -s $cigarName0 ]
		then
			perl -w -e "use VertRes::Utils::Cigar; \$c = VertRes::Utils::Cigar->new(); \$lines = \$c->cigar_to_sam(['$lane_read0'], ['$cigarName0'], undef, '$accession', 'unpaired.sam'); print \$lines; system('gzip unpaired.sam');" > unpaired.samstat
		fi/;
		close( S );
		
		#my $cmd = "bsub -J sam.$accession.0 -q normal -o bam.o -e bam.e $jobCondition perl -w -e \"use SamTools;SamTools::ssaha2samUnpaired( '$lane_read0', '$cigarName', '$accession', 'unpaired.sam.gz');\"";
		my $cmd = qq/bsub -J sam.454.$accession.0 -q $lsf_queue -o bam.o -e bam.e $jobCondition "sh unpaired.bam.sh;rm unpaired.bam.sh"/;
		#print $cmd."\n";
		system( $cmd );
	}
	
	if( ( -f $lane_read1 || -l $lane_read1 ) && ( -f $lane_read2 || -l $lane_read2 ) )
	{
		$lane_read1 =~ /^(.*)\.fastq\.gz$/;
		my $cigarName1 = $1.'.cigar.gz';
		
		$lane_read2 =~ /^(.*)\.fastq\.gz$/;
		my $cigarName2 = $1.'.cigar.gz';
		
		open( S, ">paired.bam.sh" ) or die $!;
		print S qq/if [ -s $cigarName1 ] && [ -s $cigarName2 ]
		then
			perl -w -e "use VertRes::Utils::Cigar; \$c = VertRes::Utils::Cigar->new(); \$lines = \$c->cigar_to_sam(['$lane_read1', '$lane_read2'], ['$cigarName1', '$cigarName2'], '$insertSize', '$accession', 'paired.sam'); print \$lines; system('gzip paired.sam');" > paired.samstat
		fi/;
		close( S );
		
		my $cmd = qq/bsub -J sam.454.$accession.1 -q $lsf_queue -o bam.o -e bam.e $jobCondition "sh paired.bam.sh;rm paired.bam.sh"/;
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
	
	#write the SQ headers
	my $gender = G1KUtilities::path2Gender();
	if( $gender eq 'male' || $gender eq 'unknown' )
	{
		open( FAI, $HUMAN_MALE_REF_FAI ) or die $!;
		while( <FAI> )
		{
			chomp( $_ );
			my @s = split( /\s+/ , $_ );
			print H "\@SQ\tSN:$s[0]\tAS:NCBI36\tLN:$s[1]\n";
		}
		close( FAI );
	}
	else
	{
		open( FAI, $HUMAN_FEMALE_REF_FAI ) or die $!;
		while( <FAI> )
		{
			chomp( $_ );
			my @s = split( /\s+/ , $_ );
			print H "\@SQ\tSN:$s[0]\tAS:NCBI36\tLN:$s[1]\n";
		}
		close( FAI );
	}
	
	print H "\@RG\tID:$accession\tPU:$runName\tLB:$library\tSM:$individual";
	if( $insertSize =~ /\d+/ )
	{
		print H "\tPI:$insertSize";
	}
	if( $centre )
	{
      print H "\tCN:$centre";
	}
	if( $platform )
	{
      print H "\tPL:$platform";
	}
	print H "\n";
	close( H );
	
	open( S, ">makeBam.sh" ) or die $!;
	print S qq/if [ $lane_read0 ] && [ ! -s unpaired.sam.gz ]
	then
		exit 1
	fi/;
	
	print S "\n";
	
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
	
	if( $gender eq 'male' || $gender eq 'unknown' )
	{
		print S	qq/$SAMTOOLS view -bS raw.sam | $SAMTOOLS sort - raw.sorted;/;
	}
	else
	{
		print S	qq/$SAMTOOLS view -bS raw.sam | $SAMTOOLS sort - raw.sorted;/;
	}
	
	#if the lane is single ended - then run the rmdupse instead
	if( $lane_read0 && ! $lane_read1 && ! $lane_read2 )
	{
		print S qq/$SAMTOOLS rmdupse raw.sorted.bam rmdup.bam;/;
	}
	else
	{
		print S qq/$SAMTOOLS rmdup raw.sorted.bam rmdup.bam;/;
	}
	
	print S qq/rm raw.sam;/;
	print S qq/$SAMTOOLS flagstat rmdup.bam > rmdup.flagstat;/;
	print S qq/$SAMTOOLS flagstat raw.sorted.bam > raw.sorted.flagstat;/;
	print S qq[ perl -w -e "use SamTools;SamTools::bam_stat( \\"rmdup.bam\\", \\"rmdup.bamstat\\");";];
	print S qq[ perl -w -e "use SamTools;SamTools::bam_stat( \\"raw.sorted.bam\\", \\"raw.sorted.bamstat\\");";];
	print S qq[ perl -w -e "use Mapping_454_ssaha;Mapping_454_ssaha::unmapped2Bam( \\"$lane_dir\\", \\"$indexF\\");"; ];
	close( S );
	
	my $cmd = qq/bsub -J bam.454.$accession -q $lsf_queue -o bam.o -e bam.e -w "done(sam.454.$accession.*)" "sh makeBam.sh;rm makeBam.sh;rm lane.touch"/;
	#print $cmd."\n";
	system( $cmd );
}

sub verifyBamFile
{
	croak "Usage: verifyBamFile lane_dir\n" unless @_ == 1;
	
	my $lane_dir = shift;
	
	croak "Cant find lane directory: $lane_dir\n" unless -d $lane_dir;
	
	chdir( $lane_dir );
	
	croak "Cant find meta.info file: $lane_dir\n" unless -f 'meta.info';
	
	if( -f "lane.touch" )
	{
		print "Skipping - bam in progress\n";
		return;
	}
	
	if( ! -f "raw.sorted.bam" || ! -f "raw.sorted.bamstat" )
	{
		print "Skipping - cant find bam file\n";
		return;
	}
	
	#sum the reads in the cigar mapstat files
	my @fastq = @{ HierarchyUtilities::getFastqInfo( getcwd() ) };
	my $lane_read0 = $fastq[ 0 ][ 0 ];
	my $lane_read1 = $fastq[ 1 ][ 0 ];
	my $lane_read2 = $fastq[ 2 ][ 0 ];
	
	$lane_read0 =~ /^(.*)\.fastq\.gz$/;
	my $cigarName0 = $1.'.cigar.gz';
	$lane_read1 =~ /^(.*)\.fastq\.gz$/;
	my $cigarName1 = $1.'.cigar.gz';
	$lane_read2 =~ /^(.*)\.fastq\.gz$/;
	my $cigarName2 = $1.'.cigar.gz';
	if( ( $cigarName0 && ( ! -s $cigarName0 && ! -l $cigarName0 && ! -s $cigarName0.".mapstat" && ! -l $cigarName0.".mapstat" ) ) || ( $cigarName1 && ( ! -s $cigarName1 && ! -l $cigarName1 && ! -s $cigarName1.".mapstat" && ! -l $cigarName1.".mapstat" ) ) || ( $cigarName2 && ( ! -l $cigarName2 && ! -l $cigarName2 && ! -s $cigarName2.".mapstat" && ! -l $cigarName2.".mapstat" ) ) )
	{
		print "WARNING: Cant find cigar files for lane: $lane_dir\n";
		return;
	}
	
	my $totalReadsMapped = 0;
	for my $mapstat (`ls *.mapstat`)
	{
		chomp( $mapstat );
		my $t = `grep "num_reads_mapped" $mapstat | awk -F":" '{print \$2}'`;chomp( $t );
		$totalReadsMapped += $t;
	}
	
	my $totalReadsSam = 0;
	my $totalReadsDiscarded = 0;
	for my $samstat (`ls *.samstat`)
	{
		chomp( $samstat );
		my $t = `grep "num_reads_written" $samstat | awk -F":" '{print \$2}'`;chomp( $t );
		$totalReadsSam += $t;
		$t = `grep "num_reads_discarded" $samstat | awk -F":" '{print \$2}'`;chomp( $t );
		$totalReadsDiscarded += $t;
	}
	
	my $totalReadsBam = `grep "num_reads" raw.sorted.bamstat | awk -F":" '{print \$2}'`;chomp( $totalReadsBam );
	
	if( $totalReadsSam + $totalReadsDiscarded != $totalReadsMapped || $totalReadsSam != $totalReadsBam )
	{
		print "ERROR: $totalReadsSam vs. $totalReadsMapped v.s. $totalReadsBam for lane $lane_dir\n";
	}
	else
	{
		print "GOOD $lane_dir\n";
	}
}

=head2 getLocalCopyReference

	Arg [1]	: gender (male/female)
=cut
sub getLocalCopyReference
{
	croak "Usage: getLocalCopyReference gender" unless @_ == 1;
	
	my $gender = shift;
	
	if( $gender eq 'male' || $gender eq 'unknown' )
	{
		Mapping::checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_male.base", "$REMOTE_REF_DIR/human_b36_male.base.md5" );
		Mapping::checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_male.body", "$REMOTE_REF_DIR/human_b36_male.body.md5" );
		Mapping::checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_male.head", "$REMOTE_REF_DIR/human_b36_male.head.md5" );
		Mapping::checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_male.name", "$REMOTE_REF_DIR/human_b36_male.name.md5" );
		Mapping::checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_male.size", "$REMOTE_REF_DIR/human_b36_male.size.md5" );
	}
	elsif( $gender eq 'female' )
	{
		Mapping::checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_female.Y.base", "$REMOTE_REF_DIR/human_b36_female.Y.base.md5" );
		Mapping::checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_female.Y.body", "$REMOTE_REF_DIR/human_b36_female.Y.body.md5" );
		Mapping::checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_female.Y.head", "$REMOTE_REF_DIR/human_b36_female.Y.head.md5" );
		Mapping::checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_female.Y.name", "$REMOTE_REF_DIR/human_b36_female.Y.name.md5" );
		Mapping::checkDownloadLocalFile( "$REMOTE_REF_DIR/human_b36_female.Y.size", "$REMOTE_REF_DIR/human_b36_female.Y.size.md5" );
	}
}

=head2 isLaneMapped

	Returntype : none
=cut
sub symLinkReference
{
	return unless ! -l 'ref.fa';
	
	if( getcwd() =~ /MOUSE/ )
	{
		symlink( $MOUSE_REF_FA, 'ref.fa' );
	}
	elsif( getcwd() =~ /G1K/ )
	{
		my $gender = G1KUtilities::path2Gender('/lustre/sf4/1kgenomes/ref/genders.txt');
		if( $gender eq 'male' || $gender eq 'unknown' )
		{
			symlink( $HUMAN_MALE_REF_FA, 'ref.fa' );
		}
		else
		{
			symlink( $HUMAN_FEMALE_REF_FA, 'ref.fa' );
		}
	}
}

sub getMappingStats
{
	croak "Usage: getMappingStats lane_dir" unless @_ == 1;
	
	my $dir = shift;
	
	croak "Cant find lane directory: $dir\n" unless -d $dir;
	
	chdir( $dir );
	
	#get the fastq names
	my @fastq = @{ HierarchyUtilities::getFastqInfo( $dir ) };
	
	my $rawMappedBases = 0;
	my $rawMappedReads = 0;
	my $rmdupMappedBases = 0;
	my $rmdupMappedReads = 0;
	for( my $i = 0; $i < @fastq; $i ++ )
	{
		next unless length( $fastq[ $i ][ 0 ] ) > 0;
		
		if( ! -f $fastq[ $i ][ 0 ] )
		{
			croak "Cant find fastq file: $fastq[ $i ][ 0 ]\n";
		}
		
		my @s = split( /\./, $fastq[ $i ][ 0 ] );
		my $cigarName = $s[ 0 ].'.'.$s[ 1 ].'.cigar.gz';
		if( ! -f $cigarName )
		{
			print "no cigar $fastq[ $i ][ 0 ]\n";
			return undef;
		}
		
		my $mapStat = $cigarName.'.mapstat';
		if( -f $mapStat )
		{
			my $t = `grep num_reads_mapped $mapStat | head -1 | awk -F":" '{print \$2}'`;chomp( $t );
			$rawMappedReads += $t;
			$t = `grep bases_mapped $mapStat | head -1 | awk -F":" '{print \$2}'`;chomp( $t );
			$rawMappedBases += $t;
		}
		else
		{
			print "no mapstat $mapStat\n";
			return undef;
		}
		
		if( -f "rmdup.bamstat" )
		{
			my $t = `grep 'num_reads' rmdup.bamstat | awk -F":" '{print \$2}'`;chomp( $t );
			$rmdupMappedReads += $t;
			$t = `grep 'num_reads' rmdup.bamstat | awk -F":" '{print \$2}'`;chomp( $t );
			$rmdupMappedBases += $t;
		}
		else
		{
			#return what info we have already
			my $t = [ $rawMappedReads, $rawMappedBases, 0, 0, 0, 0 ];
			return $t;
		}
	}
	
	my $t = [ $rawMappedReads, $rawMappedBases, 0, $rmdupMappedReads, $rmdupMappedBases, 0 ];
	return $t;
}

#get some stats on a cigar file
sub cigarStat
{
	croak "Usage: cigarStat cigar output" unless @_ == 2;
	
	my $cigar = shift;
	my $output = shift;
	
	croak "Cant find cigar: $cigar\n" unless -f $cigar;
	
	if( $cigar =~ /\.gz$/ )
	{
		open( CIGAR, "gunzip -c $cigar |" ) or die "Cant open cigar file: $!\n";
	}
	else
	{
		open( CIGAR, $cigar ) or die "Cant open cigar file: $!\n";
	}
	
	my %reads;
	my $basesMapped = 0;
	while( <CIGAR> )
	{
		chomp;
		
		next unless $_ =~ /^cigar::/;
		my @s = split( /\s+/, $_ );
		
		if( !defined( $reads{ $s[ 1 ] } ) )
		{
			$reads{ $s[ 1 ] } = 1;
			$basesMapped += abs( $s[ 3 ] - $s[ 2 ] );
		}
	}
	close( CIGAR );
	
	open( OUT, ">$output" ) or die "Cannot create output file: $output\n";
	print OUT "num_reads_mapped:".scalar( keys( %reads ) )."\n";
	print OUT "bases_mapped:$basesMapped\n";
	close( OUT );
}

=head2 unmapped2Bam

	Arg [1]    : lane directory
	Arg [2]    : index file (empty string if none)
	Example    : unmapped2Bam( '$G1K/MOUSE/MAPPING/Proj/SLX/Lib1/Lane1', 'normal', 'sequence.index');
	Description: Looks for the unmapped.gz file with the unmapped reads and converts them to BAM
	Returntype : none
=cut

sub unmapped2Bam
{
	croak "Usage: unmapped2Bam lane_dir index_file" unless @_ == 2;
	
	my $laneDir = shift;
	my $indexF = shift;
	
	croak "Cant find lane directory: $laneDir\n" unless -d $laneDir;
	
	chdir( $laneDir );
	
	croak "Cant find meta.info file: $laneDir\n" unless -f 'meta.info';
	
	if( ! isLaneMapped( $laneDir ) )
	{
		croak "Lane is not mapped yet! Cannot get unmapped reads\n";
	}
	
	my $accession = basename( getcwd() );
	
	my $library = basename(dirname( $laneDir )); #library is the directory name above lane
	chomp( $library );
	
	if( length( $library ) == 0 )
	{
		print "Cant find lane accession in index file - skipping unmapped BAM conversion\n";
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
	
	my $insertSize = `grep $accession $indexF | head -1 | awk -F"\t" '{print \$18}'`;
	chomp( $insertSize );
	
	my $paired = `grep $accession $indexF | grep "PAIRED" | head -1 | wc -l`;
	chomp( $paired );
	
	if( length( $individual ) == 0 || length( $runName ) == 0 )
	{
		print "Cant get index file entries: $accession\n";
		return;
	}
	
	#get the list of mapped reads from the cigar files
	my @fileInfo = @{ HierarchyUtilities::getFastqInfo( $laneDir ) };
	
	$fileInfo[ 0 ][ 0 ] =~ /^(.*)\.fastq\.gz$/;
	if( -s $1.'.cigar.gz' > 1000 )
	{
		my $cmd = qq[zcat $1.cigar.gz | grep "^cigar::" | awk '{print \$2}' > mapped0.fofn];
		#print $cmd."\n";
		system( $cmd );
	}
	
	$fileInfo[ 1 ][ 0 ] =~ /^(.*)\.fastq\.gz$/;
	if( -s $1.".cigar.gz" > 1000 )
	{
		my $cmd = qq[zcat $1.cigar.gz | grep "^cigar::" | awk '{print \$2}' > mapped1.fofn];
		#print $cmd."\n";
		system( $cmd );
	}
	
	$fileInfo[ 2 ][ 0 ] =~ /^(.*)\.fastq\.gz$/;
	if( -s $1.".cigar.gz" > 1000 )
	{
		my $cmd = qq[zcat $1.cigar.gz | grep "^cigar::" | awk '{print \$2}' > mapped2.fofn];
		#print $cmd."\n";
		system( $cmd );
	}
	
	my %mappedReads0;
	open( my $mappedFh, "mapped0.fofn" ) or die "cant open mapped0.fofn: $!\n";
	while( <$mappedFh> )
	{
		chomp;
		$mappedReads0{ $_ } = 1;
	}
	close( $mappedFh );
	
	my %mappedReads1;
	open( $mappedFh, "mapped1.fofn" ) or die "cant open mapped1.fofn: $!\n";
	while( <$mappedFh> )
	{
		chomp;
		$mappedReads1{ $_ } = 1;
	}
	close( $mappedFh );
	
	my %mappedReads2;
	open( $mappedFh, "mapped2.fofn" ) or die "cant open mapped2.fofn: $!\n";
	while( <$mappedFh> )
	{
		chomp;
		$mappedReads2{ $_ } = 1;
	}
	close( $mappedFh );
	
	unlink( "mapped0.fofn" );
	unlink( "mapped1.fofn" );
	unlink( "mapped2.fofn" );
	
	my %readPairsIncluded;
	#read each fastq file, pull out the unmapped reads, and write to sam file
	open( my $samFh, ">unmapped.raw.sam" ) or die "Cant create unmapped.raw.sam file: $!\n";
	
	#print the SAM header
	print $samFh "\@HD\tVN:1.0\tSO:coordinate\n";
  
	print $samFh "\@RG\tID:$accession\tPU:$runName\tLB:$library\tSM:$individual";
	if( $insertSize =~ /\d+/ )
	{
		print $samFh "\tPI:$insertSize";
	}
	if( $centre )
	{
		print $samFh "\tCN:$centre";
	}
	if( $platform )
	{
		print $samFh "\tPL:$platform\n";
	}
	
	if( -f $fileInfo[ 0 ][ 0 ] || -l $fileInfo[ 0 ][ 0 ] )
	{
		open( my $read0, "gunzip -c $fileInfo[ 0 ][ 0 ] |" ) or die "cant open fastq file: $!\n";
		while( <$read0> )
		{
			chomp;
			my $rname = substr( $_, 1 );
			my $seq = <$read0>;
			chomp( $seq );
			my $qname = <$read0>;
			chomp( $qname );
			my $quals = <$read0>;
			chomp( $quals );
			
			if( ! defined( $mappedReads0{ ( split( /\s+/, $rname ) )[ 0 ] } ) )
			{
				my $flags = SamTools::determineUnmappedFlag( $paired, -1 );
				print $samFh ( split( /\s+/, $rname ) )[ 0 ]."\t$flags\t*\t0\t0\t*\t*\t0\t0\t$seq\t$quals\tRG:Z:$accession\n";
			}
		}
		close( $read0 );
	}
	
	if( ( -f $fileInfo[ 1 ][ 0 ] || -l $fileInfo[ 1 ][ 0 ] ) && ( -f $fileInfo[ 2 ][ 0 ] || -l $fileInfo[ 2 ][ 0 ] ) )
	{
		open( my $read1, "gunzip -c $fileInfo[ 1 ][ 0 ] |" ) or die "cant open fastq file: $!\n";
		while( <$read1> )
		{
			chomp;
			my $rname = substr( $_, 1 );
			my $seq = <$read1>;
			chomp( $seq );
			my $qname = <$read1>;
			chomp( $qname );
			my $quals = <$read1>;
			chomp( $quals );
			
			if( ! defined( $mappedReads1{ ( split( /\s+/, $rname ) )[ 0 ] } ) )
			{
				my $flags = SamTools::determineUnmappedFlag( $paired, -1 );
				print $samFh ( split( /\s+/, $rname ) )[ 0 ]."\t$flags\t*\t0\t0\t*\t*\t0\t0\t$seq\t$quals\tRG:Z:$accession\n";
			}
		}
		close( $read1 );
		
		open( my $read2, "gunzip -c $fileInfo[ 2 ][ 0 ] |" ) or die "cant open fastq file: $!\n";
		while( <$read2> )
		{
			chomp;
			my $rname = substr( $_, 1 );
			my $seq = <$read2>;
			chomp( $seq );
			my $qname = <$read2>;
			chomp( $qname );
			my $quals = <$read2>;
			chomp( $quals );
			
			if( ! defined( $mappedReads2{ ( split( /\s+/, $rname ) )[ 0 ] } ) )
			{
				my $flags = SamTools::determineUnmappedFlag( $paired, -1 );
				print $samFh ( split( /\s+/, $rname ) )[ 0 ]."\t$flags\t*\t0\t0\t*\t*\t0\t0\t$seq\t$quals\tRG:Z:$accession\n";
			}
		}
		close( $read2 );
	}
	close( $samFh );
	
	#convert to BAM
	my $cmd = qq[ $SAMTOOLS view -bS raw.unmapped.sam > raw.unmapped.bam;$SAMTOOLS flagstat raw.unmapped.bam > raw.unmapped.bam.flagstat];
	system( $cmd );
}

1;
