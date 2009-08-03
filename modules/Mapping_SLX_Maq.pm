package Mapping_SLX_Maq;
use strict;
use Carp;
use File::Basename;
use File::Spec;
use Cwd;

use G1KUtilities;
use AssemblyTools;
use SamTools;

my $MAQ = "maq";
my $SAMTOOLS = "samtools";
my $BASES_PER_CHUNK = 57600000;
my $MAQ2SAM = "maq2sam-long";

#read length ranges and what to set the -e parameter in maq to
my %eParameter =
(
		37 	=> 	80, 	#2 * 40
		63	=> 	120,	#3 * 40
		92 	=> 	160, 	#4 * 40
		123	=>	200, 	#5 * 40
		156	=>	240 	#6 * 40
);

=pod

=head1 NAME

Mapping_SLX_Maq

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
	
	if( -f 'lane.touch' )
	{
		return 1;
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
	
	if( -f 'raw.map' && -s 'raw.map' > 2000 && -s 'raw.map.mapstat' && ! -d 'split' )
	{
		return 1;
	}
	return 0;
}

=head2 unmapped2Bam

	Arg [1]    : lane directory
	Arg [2]    : index file (empty string if none)
	Example    : unmapped2Bam( '$G1K/MOUSE/MAPPING/Proj/SLX/Lib1/Lane1', 'normal', 'sequence.index');
	Description: Looks for the unmapped.gz file with the unmapped reads and converts them to BAM. Note this function may take a while to run as it calls samtools import - therefore it might be advisable to call the function via LSF
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
	
	#make a sam file first
	open( my $umFh, "zcat unmapped.gz |" ) or die "Cant open unmapped.gz file: $!\n";
	
	#write the sam header
	open( my $sam, ">raw.unmapped.sam" ) or die "Cannot create raw.sam: $!\n";
	print $sam "\@HD\tVN:1.0\n";
	
	print $sam "\@RG\tID:$accession\tPU:$runName\tLB:$library\tSM:$individual";
	if( $insertSize =~ /\d+/ )
	{
		print $sam "\tPI:$insertSize";
	}
	if( $centre )
	{
		print $sam "\tCN:$centre";
	}
	if( $platform )
	{
      print $sam "\tPL:$platform";
	}
	
	print $sam "\n";
	
	my %readsIncluded;
	while( <$umFh> )
	{
		chomp;
		my @s = split( /\s+/, $_ );
		
		#work out the SAM flags
		my $flags;
		if( defined( $readsIncluded{ $s[ 0 ] } ) )
		{
			$flags = SamTools::determineUnmappedFlag( $paired, 0 );
		}
		else
		{
			$flags = SamTools::determineUnmappedFlag( $paired, 1 );
			$readsIncluded{ $s[ 0 ] } = 1;
		}
		
		print $sam qq/$s[ 0 ]\t$flags\t*\t0\t0\t*\t*\t0\t0\t$s[ 2 ]\t$s[ 3 ]\tRG:Z:$accession\n/;
	}
	close( $umFh );
	close( $sam );
	
	#convert to BAM
	my $cmd = qq[ $SAMTOOLS view -bS raw.unmapped.sam > raw.unmapped.bam; $SAMTOOLS flagstat raw.unmapped.bam > raw.unmapped.bam.flagstat ];
	system( $cmd );
	$cmd = qq[ rm raw.unmapped.sam ];
	system( $cmd );
}

=head2 mapLane

	Arg [1]    : lane directory
	Arg [2]    : LSF queue for jobs
	Arg [3]    : index file (empty string if none)
	Example    : map( '$G1K/MOUSE/MAPPING/Proj/SLX/Lib1/Lane1', 'normal', 'sequence.index');
	Description: Maps the lane specified
	Returntype : none
=cut

sub mapLane
{
	croak "Usage: mapLane lsf_queue lanes_fofn index_file" unless @_ == 3;
	
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
	
	my $laneID = basename( $laneDir ).int(rand(10000000000));
	
	#check if the lane is already mapped
	if( -f "raw.map" && -s "raw.map" > 2000 )
	{
		print "Lane is already mapped\n";
		
		#only do the bam conversion when there is an index file (i.e. doesnt apply to mouse)
		if( -f $indexF )
		{
			#laneToBAM( getcwd(), $laneID, $lane_read0, $lane_read1, $lane_read2, basename( getcwd() ), $indexF, $lsf_queue );
		}
		return;
	}
	
	if( ! -l 'ref.bfa' && ! -f 'ref.bfa' )
	{
		#create a sym link to the ref bfa
		croak "Error: Cant find ref.bfa file!";
	}
    
	my $mapping = 0;
	my $mergeNoMapping = 0;
    my $maqE = "";
	
	if( $lane_read1 && $lane_read2 )
	{
		if( $lane_read1 && ! -f "$lane_read1.fastqcheck" && ! -l "$lane_read1.fastqcheck" )
		{
			print "WARNING: Cant find fastqcheck for $lane_read1\n";
			chdir( ".." );
			next;
		}
		
		if( $lane_read1 && ! -f "$lane_read2.fastqcheck" && ! -l "$lane_read2.fastqcheck" )
		{
			print "WARNING: Cant find fastqcheck for $lane_read2\n";
			chdir( ".." );
			next;
		}
		
		#get the read lengths and adjust the maq -e paramter if necessary
		my $readLen1 = `cat $lane_read1.fastqcheck | head -1 | awk '{print \$8}'`;chomp( $readLen1 );
		my $readLen2 = `cat $lane_read1.fastqcheck | head -1 | awk '{print \$8}'`;chomp( $readLen2 );
		
		#work out what -e parameter to use for maq
		foreach( sort {$a <=> $b}( keys( %eParameter ) ) )
		{
			if( $readLen1 < $_ && $readLen2 < $_ )
			{
				$maqE = "-e $eParameter{ $_ }";
				last;
			}
		}
	}
	elsif( $lane_read0 )
	{
		if( ! -f "$lane_read0.fastqcheck" && ! -l "$lane_read0.fastqcheck" )
		{
			print "WARNING: Cant find fastqcheck for $lane_read0\n";
			chdir( ".." );
			next;
		}
		
		#get the read lengths and adjust the maq -e paramter if necessary
		my $readLen0 = `cat $lane_read0.fastqcheck | head -1 | awk '{print \$8}'`;chomp( $readLen0 );
		
		#work out what -e parameter to use for maq
		foreach( sort {$a <=> $b}( keys( %eParameter ) ) )
		{
			if( $readLen0 < $_ )
			{
				$maqE = "-e $eParameter{ $_ }";
				last;
			}
		}
	}
	else
	{
		print "WARNING: Found no reads in meta.info file\n";
		chdir( ".." );
		next;
	}
	
	#check if the lane fastq splits have been made and files are not empty (i.e. the splitting job finished ok)
	if( -d "split" )
	{
		foreach my $s (`ls split/*.fastq`)
		{
			chomp $s;
			if( ! -f $s || -z $s ) #each split fastq file is non-zero size
			{
				system( "rm -rf split" ); #delete the split directory
				last;
			}
		}
	}
    
	my $lane_dir = getcwd();
	
	#if no mapping has been done
	if( ( -f "raw.map" && -s "raw.map" < 2000 ) || ( ! -f "raw.map" && ! -d "split" ) )
	{
		system( "touch lane.touch" );
        
		mkdir( "split" );
		
		#split the reads fastq
		if( $lane_read1 && $lane_read2 )
		{
			my $cmd = qq[bsub -J split.$laneID -q $lsf_queue -o $lane_dir/split/split.o -e $lane_dir/split/split.e perl -w -e "use AssemblyTools;AssemblyTools::splitPairedFastq( '$lane_read1', '$lane_read2', 'split', '$BASES_PER_CHUNK', '$lane_dir/split' );"];
			#print "$cmd\n";
			system( $cmd );
		}
		elsif( $lane_read0 )
		{
			my $cmd = qq[bsub -J split.$laneID -q $lsf_queue -o $lane_dir/split/split.o -e $lane_dir/split/split.e perl -w -e "use AssemblyTools;AssemblyTools::splitUnpairedFastq( '$lane_read0', 'split', '$BASES_PER_CHUNK', '$lane_dir/split' );"];
			#print "$cmd\n";
			system( $cmd );
		}
		else
		{
			print "Cant determine if lane is paired or unpaired\n";
		}
        
		chdir( "split" );
		
		if( $lane_read1 && $lane_read2 )
		{
			#determine how many splits there will be
			my $totalLengthR1=`cat ../$lane_read1.fastqcheck | head -1 | awk '{print \$3}'`;chomp( $totalLengthR1 );
			my $totalLengthR2=`cat ../$lane_read2.fastqcheck | head -1 | awk '{print \$3}'`;chomp( $totalLengthR2 );
			
			my $totalLength = $totalLengthR1+$totalLengthR2;
			my $numSplits= int($totalLength/$BASES_PER_CHUNK) + 1;
			print "Num splits: $numSplits\n";
            
			#check if the clip points exist - else dont clip
			my $clipPoint1=-1;
			my $clipPoint2=-1;
			if( length( `grep "clip1" ../meta.info | wc -l`) > 1 && length( `grep "clip2" ../meta.info | wc -l`) > 1 )
			{
				$clipPoint1=`grep "clip1" ../meta.info | awk -F ":" '{print \$2}' | head -1`;chomp( $clipPoint1 );if( length( $clipPoint1 ) == 0 ){$clipPoint1=-1;}
				$clipPoint2=`grep "clip2" ../meta.info | awk -F ":" '{print \$2}' | head -1`;chomp( $clipPoint2 );if( length( $clipPoint2 ) == 0 ){$clipPoint2=-1;}
			}

			for( my $s=0;$s<$numSplits;$s++ )
			{
				my $read1="split$s"."_1.fastq";
				my $read2="split$s"."_2.fastq";
				my $prefix = (split( /_/, $read1 ) )[ 0 ];
				
				#map the split
				my $cmdPrefix = qq[bsub -w "done(split.$laneID)" -J map.slx.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -u $prefix.unmapped $maqE -a 1000];
				
				my $cmdMiddle = '';
				if( $clipPoint1 != -1 )
				{
					$cmdMiddle .= qq[ -1 $clipPoint1 ];
				}
				
				if( $clipPoint2 > -1 )
				{
					$cmdMiddle .= qq[ -2 $clipPoint2 ];
				}
				
				if( $expectedInsert && $expectedInsert > 1000 )
				{
					my $ins = $expectedInsert * 2;
					$cmdMiddle .= qq[ -A $ins ];
				}
				
				my $cmdSuffix = qq[$prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
				
				my $cmd = qq[$cmdPrefix $cmdMiddle $cmdSuffix];
				system( $cmd );
                $mapping=1;
			}
		}
		else #unpaired lane
		{
			#determine how many splits there will be
			my $totalLength=`cat ../$lane_read0.fastqcheck | head -1 | awk '{print \$3}'`;chomp( $totalLength );
			my $numSplits= int($totalLength/$BASES_PER_CHUNK) + 1;
			print "Num splits: $numSplits\n";
			
            my $clipPoint0=-1;
            if( length( `grep "clip0" ../meta.info | wc -l` ) > 1 )
            {
            	$clipPoint0=`grep clip0 ../meta.info | awk -F ":" '{print \$2}' | head -1`;if( length( $clipPoint0 ) == 0 ){$clipPoint0=-1;}
                chomp( $clipPoint0 );
			}
			
			for(my $s=0;$s<$numSplits;$s++)
            {
				my $read0="split$s.fastq";
				my $prefix = (split( /\./, $read0 ) )[ 0 ];
				
				#map the split
				my $cmdPrefix = qq[bsub -w "done(split.$laneID)" -J map.slx.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read0 $read0.bfq; $MAQ match -u $prefix.unmapped $maqE];
				my $cmdMiddle = '';
				if( $clipPoint0 > -1 )
				{
					$cmdMiddle = qq[-1 $clipPoint0];
				}
				
				my $cmdSuffix = qq[$prefix.raw.map ../ref.bfa $read0.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read0.bfq"];
				my $cmd = qq[$cmdPrefix $cmdMiddle $cmdSuffix];
				system( $cmd );
                $mapping=1;
            }
		}
		
		chdir( ".." );
    }
	elsif( ( -f "raw.map" && -s "raw.map" < 2000 && -d "split" ) || ( ! -f "raw.map" && -d "split" ) )
    {
		system( "touch lane.touch" );

		#the lane splits are partially mapped
		chdir( "split" );
		if( $lane_read1 && $lane_read2 )
		{
			foreach my $read1 (`ls split*_1.fastq`)
			{
				chomp $read1;
				if( ! -f $read1 )
				{
					last;
				}
				
				my $prefix = (split( /_/, $read1 ) )[ 0 ];
				my $read2=$prefix."_2.fastq";
				
				#if already mapped - skip
				if( -f "$prefix.raw.map" && -s "$prefix.raw.map" > 2000 )
				{
					print "Split $prefix already mapped\n";
					next;
				}
				
				#check if the clip points exist - else dont clip
				my $clipPoint1=-1;
				my $clipPoint2=-1;
				if( length( `grep "clip1" ../meta.info | wc -l` ) > 1 )
				{
					$clipPoint1=`grep "clip1" ../meta.info | awk -F ":" '{print \$2}' | head -1`;chomp( $clipPoint1 );if( length( $clipPoint1 ) == 0 ){$clipPoint1=-1;}
					$clipPoint2=`grep "clip2" ../meta.info | awk -F ":" '{print \$2}' | head -1`;chomp( $clipPoint2 );if( length( $clipPoint2 ) == 0 ){$clipPoint2=-1;}
				}
                
				my $cmdPrefix = qq[bsub -J map.slx.$laneID.$read2 -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -u $prefix.unmapped $maqE -a 2000];
				
				my $cmdMiddle = '';
				if( $clipPoint1 != -1 )
				{
					$cmdMiddle .= qq[ -1 $clipPoint1 ];
				}
				
				if( $clipPoint2 > -1 )
				{
					$cmdMiddle .= qq[ -2 $clipPoint2 ];
				}
				
				if( $expectedInsert && $expectedInsert > 1000 )
				{
					my $ins = $expectedInsert * 2;
					$cmdMiddle .= qq[ -A $ins ];
				}
				
				my $cmdSuffix = qq[$prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
				
				system( qq[$cmdPrefix $cmdMiddle $cmdSuffix] );
				$mapping=1;
			}
			$mergeNoMapping = 1
		}
		elsif( $lane_read0 )
		{
				my $s = 0;
				foreach my $read0 (`ls split*.fastq`)
				{
					chomp $read0;
					if( ! -f $read0 )
					{
						last;
					}
					
					my $prefix = (split( /\./, $read0 ) )[ 0 ];
					my $clipPoint0=-1;
					if( length( `grep "clip0" ../meta.info | wc -l` ) > 1 )
					{
						$clipPoint0=`grep clip0 ../meta.info | awk -F ":" '{print \$2}' | head -1`;chomp( $clipPoint0 );if( length( $clipPoint0 ) == 0 ){$clipPoint0=-1;}
					}
					
					if( ! -f "$prefix.raw.map" || ( -f "$prefix.raw.map" && -s "$prefix.raw.map" < 2000 ) )
					{
						my $cmdPrefix = qq[bsub -J map.slx.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read0 $read0.bfq; $MAQ match];
						my $cmdMiddle = '';
						if( $clipPoint0 > -1 )
						{
							$cmdMiddle = qq[-1 $clipPoint0];
						}
						my $cmdSuffix = qq[-u $prefix.unmapped $maqE $prefix.raw.map ../ref.bfa $read0.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read0.bfq"];
                        system( qq[$cmdPrefix $cmdMiddle $cmdSuffix] );
						$mapping=1;
					}
				}
				$mergeNoMapping = 1;
		}
		chdir( ".." );
	}
	else
	{
		croak "Unknown state";
	}
	
	#check all the map files are there
	if( $mapping == 1 || $mergeNoMapping == 1 )
	{
				#write out a little script to merge the lane splits
				open( S, ">check_and_merge.sh" ) or die $!;
				if ( $lane_read1 && $lane_read2 )
				{
					print S "totalMapJobs=`ls split/*_1.fastq | wc -l`\n";
				}
				else
				{
					print S "totalMapJobs=`ls split/*.fastq | wc -l`\n";
				}
                        
				print S qq~
					completed=0
					
                        for map in split/*.raw.map
                        do
                          size=`ls -l \$map | awk '{print \$5}'`
                          if \[ \$size -gt 2000 \]
                          then
                            ((completed=completed+1))
                          fi
                        done
                        
                        #if all the jobs are completed - then merge the maps and delete the split directory
                        if \[ \$completed -eq \$totalMapJobs \]
                        then
                          cat split/*.unmapped | gzip -c > unmapped.gz
                          $MAQ mapmerge raw.map split/*.raw.map
                          $MAQ mapstat raw.map > raw.map.mapstat
                          $MAQ rmdup rmdup.map raw.map
                          $MAQ mapstat rmdup.map > rmdup.map.mapstat
                          rm -rf split
                        else
							echo "Number of splits mapped not equal to number of jobs: \$completed \$totalMapJobs"
						fi~;
                        close( S );
                        
						sleep( 5 );
						my $numMapJobs = `bjobs -w | grep "map.slx.$laneID." | wc -l`;chomp( $numMapJobs );
						my $jobCondition = '';
						if( $numMapJobs > 0 )
						{
							$jobCondition = '-w "done(map.slx.'.$laneID.'.*)"';
						}
						
                        my $cmd = qq[bsub -J merge.slx.$laneID -q $lsf_queue -o $lane_dir/merge.o -e $lane_dir/merge.e $jobCondition "sh check_and_merge.sh; rm check_and_merge.sh; rm lane.touch"];
                        #print "$cmd\n";
                        system( $cmd );
	}
	
	if( -f $indexF )
	{
		#laneToBAM( getcwd(), $laneID, $lane_read0, $lane_read1, $lane_read2, basename( getcwd() ), $indexF, $lsf_queue );
	}
}

=head2 mapLane

	Arg [1]    : lane directory
	Arg [2]    : LSF queue for jobs
	Arg [3]    : sequence index file
	Arg [4]    : male fai file
	Arg [5]    : female fai file
	Example    : laneToBAM( '$G1K/MOUSE/MAPPING/Proj/SLX/Lib1/Lane1', 'normal', 'sequence.index', 'male.fai', 'female.fai' );
	Description: Converts the raw.map file to Bam
	Returntype : none
=cut

sub laneToBAM
{
    croak "Usage: laneToBAM lane_dir lsf_queue index_file male_fai female_fai\n" unless @_ == 5;
	
	my $lane_dir = shift;
	my $lsf_queue = shift;
	my $indexF = shift;
	my $male_fai = shift;
	my $female_fai = shift;
	
	croak "Cant find index file: $indexF\n" unless -f $indexF;
	croak "Cant find lane directory: $lane_dir\n" unless -d $lane_dir;
	croak "Cant find male fai file: $!\n" unless -f $male_fai;
	croak "Cant find female fai file: $!\n" unless -f $female_fai;
	
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
	
	if( -s "raw.sorted.bam" )
	{
		print "Skipping - found existing raw.sorted.bam file\n";
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
	#check if there are pending merging job for the lane
	if( `bjobs -w | grep "merge.slx.$accession" | wc -l` > 0 )
	{
		$jobCondition = ' -w "done( merge.slx.'.$accession.'.* )"';
	}
	elsif( ! -f "raw.map" )
	{
		print "No jobs running and no raw.map file\n";
		return;
	}
	
	my $tmp_directory = "/tmp/".int( rand( 1000000) );
	
	#write the SAM header
	open( H, ">raw.sam" ) or die "Cannot create raw.sam: $!\n";
	print H "\@HD\tVN:1.0\tSO:coordinate\n";
	
	my $REF_FAI;
	my $gender = G1KUtilities::path2Gender();
	if( $gender eq 'male' || $gender eq 'unknown' )
	{
		$REF_FAI = $male_fai;
	}
	else
	{
		$REF_FAI = $female_fai;
	}
	
	open( my $faiH, $REF_FAI ) or die "cant open ref fai file: $REF_FAI $!\n";
	while( <$faiH> )
  {
	  chomp;
	  my @s = split( /\s+/, $_ );
	  print H "\@SQ\tSN:$s[ 0 ]\tAS:NCBI36\tLN:$s[ 1 ]\n";
  }
  close( $faiH );
  
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
  system( "touch lane.touch" );
  
  unlink( "bam.o" ) unless ! -f "bam.o";
  unlink( "bam.e" ) unless ! -f "bam.e";
  
  #convert to bam
  my $cmd = qq[bsub -R "select[type==X86_64 && mem > 6000] rusage[mem=6000:tmp=15000]" -M6000000 -J bam.slx.$accession -q $lsf_queue $jobCondition -o bam.o -e bam.e 'mkdir $tmp_directory;$MAQ2SAM raw.map $accession > $tmp_directory/tmp.sam; cat raw.sam $tmp_directory/tmp.sam > $tmp_directory/raw.sam; rm $tmp_directory/tmp.sam;rm raw.sam;];
  
  $cmd .= qq[$SAMTOOLS view -bS $tmp_directory/raw.sam | $SAMTOOLS sort -n - $tmp_directory/raw.name_sort; samtools fixmate $tmp_directory/raw.name_sort.bam - | $SAMTOOLS sort - $tmp_directory/raw.sorted;];
  
  if( $lane_read0 && ! $lane_read1 && ! $lane_read2 )
  {
	  $cmd .= qq[ $SAMTOOLS rmdupse $tmp_directory/raw.sorted.bam $lane_dir/rmdup.bam; ];
  }
  else
  {
	  $cmd .= qq[ $SAMTOOLS rmdup $tmp_directory/raw.sorted.bam $lane_dir/rmdup.bam; ];
  }
  
  $cmd .= qq[ $SAMTOOLS flagstat $tmp_directory/raw.sorted.bam > $lane_dir/raw.sorted.flagstat ;];
  $cmd .= qq[ mv $tmp_directory/raw.sorted.bam $lane_dir/raw.sorted.bam; ];
  $cmd .= qq[ perl -w -e "use SamTools;SamTools::bam_stat( \\"$lane_dir/raw.sorted.bam\\", \\"$lane_dir/raw.sorted.bamstat\\");";];
  
  $cmd .= qq[ $SAMTOOLS index $lane_dir/rmdup.bam;];
  $cmd .= qq[ $SAMTOOLS flagstat $lane_dir/rmdup.bam > $lane_dir/rmdup.flagstat;];
  $cmd .= qq[ perl -w -e "use SamTools;SamTools::bam_stat( \\"$lane_dir/rmdup.bam\\", \\"$lane_dir/rmdup.bamstat\\");";];
  if( -f $lane_dir."/unmapped.gz" && -s $lane_dir."/unmapped.gz" > 100 )
  {
	  $cmd .= qq[ perl -w -e "use Mapping_SLX_Maq;Mapping_SLX_Maq::unmapped2Bam( \\"$lane_dir\\", \\"$indexF\\");"; ];
  }
  $cmd .= qq[ rm -rf $tmp_directory ; rm lane.touch'];
  
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
		print "Skipping - bam in progress: $lane_dir\n";
		return;
	}
	
	if( ! -f "raw.sorted.bam" || ! -f "raw.sorted.bamstat" )
	{
		print "Skipping - cant find bam file: $lane_dir\n";
		return;
	}
	
	if( ! -f "raw.map.mapstat" && ! -l "raw.map.mapstat" )
	{
		print "WARNING: Cant find mapstat file: $lane_dir\n";
		return;
	}
	
	my $raw_map = `grep "Total number of reads" raw.map.mapstat | awk '{print \$6}'`;chomp($raw_map);
	my $raw_bam = `grep "num_reads" raw.sorted.bamstat | awk -F":" '{print \$2}'`;chomp( $raw_bam );
	
	if( $raw_map > $raw_bam * 1.05 )
	{
		print "ERROR: $raw_map vs. $raw_bam for lane $lane_dir\n";
	}
	else
	{
		print "GOOD $lane_dir\n";
	}
	
}

sub removeIntermediateMappingFiles
{
	croak "Usage: removeMappingFiles lane_dir" unless @_ == 1;
	
	my $laneDir = shift;
	croak "Cant find lane directory: $laneDir\n" unless -d $laneDir;
	
	if( isLaneMapped( $laneDir ) )
	{
		print "Lane is fully mapped - skipping";
		return;
	}
	
	chdir( $laneDir );
	
	unlink 'lane.touch' unless ! -f 'lane.touch';
	
	system( "rm -rf split*" );
	
	unlink "map.o"  unless ! -f "map.o";
	
	unlink "map.e"  unless ! -f "map.e";
	
	unlink "merge.o"  unless ! -f "merge.o";
	
	unlink "merge.e"  unless ! -f "merge.e";
}

sub calculateClipPointsLane
{
	croak "Usage: calculateClipPointsLane lane_dir" unless @_ == 1;
	
	my $laneDir = shift;
	croak "Cant find lane directory: $laneDir\n" unless -d $laneDir;
	
	chdir( $laneDir );
	
	croak "Cant find meta.info file in $laneDir\n" unless -f 'meta.info';
	
	my $alreadyClipped = `grep "clip" meta.info | wc -l`;chomp( $alreadyClipped );
	
	if( $alreadyClipped > 0 || -f 'raw.map' )
	{
		return;
	}
	
	print "Calculating clip points for lane: $laneDir\n";
	
	my $fastq = HierarchyUtilities::getFastqNamesLane( $laneDir );
	my $counter = 0;
	my $laneID = basename( getcwd() );
	foreach( @{ $fastq } )
	{print $_."\n";
		if( -f $_ )
		{
			my $cmd = qq[bsub -J clip.$laneID.$counter -q normal -o $laneDir/import.o -e $laneDir/import.e perl -w -e "use Mapping_SLX_Maq;Mapping_SLX_Maq::calculateClipPointsFastq( '$laneDir/$_', '$laneDir/meta.info', '$counter' );"];
			#print $cmd."\n";
			system( $cmd );
		}
		$counter ++;
	}
}

sub calculateClipPointsFastq
{
	croak "Usage: calculateClipPointsFastq fastq meta_file fastq_counter" unless @_ == 3;
	
	my $fastq = shift;
	my $meta = shift;
	my $counter = shift;
	
	croak "Cant find fastq file $fastq\n" unless -f $fastq;
	croak "Cant find meta.info file\n" unless -f $meta;
	croak "Fastq counter must be 0,1,2\n" unless $counter >-1 && $counter < 3;
	
	my $point = AssemblyTools::determineClipPointMaq( "$fastq" );
	system( "echo 'clip$counter:$point' >> $meta" );
}

sub getMappingStats
{
	croak "Usage: getMappingStats lane_dir" unless @_ == 1;
	
	my $dir = shift;
	
	croak "Cant find lane directory: $dir\n" unless -d $dir;
	
	chdir( $dir );
	
	my $rawMappedReads = 0;
	my $rawMappedBases = 0;
	my $rmdupMappedReads = 0;
	my $rmdupMappedBases = 0;
	my $errorRate = 0;
	my $rawNumReads = 0;
	my $mappedNumReads = 0;
	my $rawPairedNumReads = 0;
	if( -s "raw.map.mapstat" )
	{
		my $isLastIteration = 0;
		my $lastFileName = '';
=pod		
		my @outputs = (`ls -l split*.maq.out.gz | sort -k 4 -r | awk '{print \$9}'`);

		for( my $i = 0; $i < @outputs; $i ++ )
		{
			chomp( $outputs[ $i ] );
			
			my $line=`zcat $outputs[ $i ] | tail -5 | grep '(total, isPE, mapped, paired)'`;
			my $zeroReads=`zcat $outputs[ $i ] | grep "reads loaded"`;
			
			if( $zeroReads !~ /\s+0\*2\s+/ )
			{
				if( length( $line ) > 0 )
				{
					$line =~ /.*=\s+\((\d+),\s+(\d+),\s+(\d+),\s+(\d+)\)/;
					
					$rawNumReads += $1;
					$mappedNumReads += $3;
					$rawPairedNumReads += $4;
				}
				else
				{
					print "1BAD_LANE_SPLIT: $outputs[ $i ] ".getcwd."\n";
					return undef;
				}
			}
		}
=cut		
		#verify the number of reads agrees with the mapstat file
		$rawMappedReads = `head raw.map.mapstat | grep 'Total number of reads' | awk -F":" '{print \$2}'`; $rawMappedReads =~ s/\s//g;
		$rawMappedBases = `head raw.map.mapstat | grep 'Sum of read length' | awk -F":" '{print \$2}'`; $rawMappedBases =~ s/\s//g;
		$errorRate = `head raw.map.mapstat | grep 'Error rate' | awk -F":" '{print \$2}'`; $errorRate =~ s/\s//g;
		if( -s "rmdup.map.mapstat" )
		{
			$rmdupMappedReads = `head rmdup.map.mapstat | grep 'Total number of reads' | awk -F":" '{print \$2}'`; $rmdupMappedReads =~ s/\s//g;
			$rmdupMappedBases = `head rmdup.map.mapstat | grep 'Sum of read length' | awk -F":" '{print \$2}'`; $rmdupMappedBases =~ s/\s//g;
		}
		else
		{
			print "No rmdup.map.mapstat: ".getcwd."\n";
			return undef;
		}
=pod		
		if( $rawMappedReads != $mappedNumReads )
		{
			print "Number of reads in mapstat not equal to split outputs: $rawMappedReads vs. $mappedNumReads\n";
			print "BAD_LANE_MAPSTAT: ".getcwd."\n";
			return undef;
		}
=cut		
		if( -f "raw.sorted.bam" && -s "raw.sorted.flagstat" )
		{
			my $numReadsInBam = `head -1 raw.sorted.flagstat | awk '{print \$1}'`;
			
			if( $rawMappedReads - $numReadsInBam > ( $rawMappedReads * 0.05 ) )
			{
				return undef;
			}
			
			$rawPairedNumReads = `cat raw.sorted.flagstat | grep 'properly paired' | awk '{print \$1}'`; $rawPairedNumReads =~ s/\s//g;
		}
		else
		{
			print "INCOMPLETE: No raw.sorted.bam file\n";
			return undef;
		}
	}
	else
	{
		print "INCOMPLETE: Can't find mapstat file\n";
		return undef;
	}
	
	my $t = [ $rawMappedReads, $rawMappedBases, $rawPairedNumReads, $rmdupMappedReads, $rmdupMappedBases, $errorRate ];
	return $t;
}

1;
