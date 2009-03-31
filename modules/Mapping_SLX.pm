package Mapping_SLX;
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
my $MAQ = $G1K."/bin/maq";
my $SAMTOOLS = $G1K."/bin/samtools";
my $BASES_PER_CHUNK = 57600000;
my $MAQ2SAM = $G1K."/bin/maq2sam-long";

my $FEMALE_REF_FAI = $G1K.'/ref/human_b36_female.fa.fai';
my $MALE_REF_FAI = $G1K.'/ref/human_b36_male.fa.fai';

sub map_SLX
{
	croak "Usage: map_SLX mappingRootDir projects_list lsf_queue male_ref_fa male_ref_bfa female_ref_fa female_ref_bfa [index_file]" unless @_ != 8 || @_ != 7;
	
	my $mRoot = shift;
	my $projs = shift;
	my $lsf_queue = shift;
	
	my $male_ref_fa = shift;
	my $male_ref_bfa = shift;
	my $female_ref_fa = shift;
	my $female_ref_bfa = shift;
	
	my $indexF = '';
	if( @_ == 1 )
	{
		$indexF = shift;
		croak "Cant find index file: $indexF\n" unless -f $indexF;
	}
	
	croak "Cant find mapping directory: $mRoot\n" unless -d $mRoot;
	croak "Invalid LSF queue: $lsf_queue\n" unless $lsf_queue eq 'normal' || $lsf_queue eq 'yesterday' || $lsf_queue eq 'basement';
	
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
						if( -d $technology && $technology eq 'SLX' )
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
									
									foreach my $lane (`ls`) 
									{
										chomp( $lane );
										
										if( -d $lane )
										{
											chdir( $lane );
											print "In ".getcwd."\n";
											
											my $laneID = int( rand( 1000000 ) );
											
											#if the directory is being operated on at the moment - skip it
											if( -f "lane.touch" )
											{
												chdir( ".." );
												print "Lane is currently being mapped - skipping\n";
												next;
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
											
											#check if the lane is already mapped
											if( -f "raw.map" && -s "raw.map" > 2000 )
											{
												print "Lane is already mapped\n";
												
												#only do the bam conversion when there is an index file (i.e. doesnt apply to mouse)
												if( length( $indexF ) > 0 )
												{
													laneToBAM( getcwd(), $libID, $laneID, $lane_read0, $lane_read1, $lane_read2, basename( getcwd() ), $gender, $indexF );
												}
												
												chdir( ".." );
												next;
											}
											chdir( ".." );
											next;
											
											#create a sym link to the ref bfa
											if( $gender eq 'male' || $gender eq 'unknown' )
											{
												symLinkReference( $male_ref_fa, $male_ref_bfa );
											}
											else
											{
												symLinkReference( $female_ref_fa, $female_ref_bfa );
											}
											
											my $mapping = 0;
											
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
											}	
											elsif( $lane_read0 )
											{
												if( ! -f "$lane_read0.fastqcheck" && ! -l "$lane_read0.fastqcheck" )
												{
													print "WARNING: Cant find fastqcheck for $lane_read0\n";
													chdir( ".." );
													next;
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
													if( ! -f $s || -z $s ) #each split fastq file is non-zero size
													{
														#system( "rm -rf split" ); #delete the split directory
														last;
													}
												}
											}
											
											my $lane_dir = getcwd();
					
											#if no mapping has been done
											if( ( -f "raw.map" && -s "raw.map" < 2000 ) || ( ! -f "raw.map" && ! -d "split" ) )
											{
												system( "echo 'touched' > lane.touch" );
					
												mkdir( "split" );
					
												#split the reads fastq
												if( $lane_read1 && $lane_read2 )
												{
													my $cmd = qq[bsub -J split.$libID.$laneID -q $lsf_queue -o $lane_dir/split/split.o -e $lane_dir/split/split.e perl -w -e "use AssemblyTools;AssemblyTools::splitPairedFastq( '$lane_read1', '$lane_read2', 'split', '$BASES_PER_CHUNK', '$lane_dir/split' );"];
													print "$cmd\n";
												}
												elsif( $lane_read0 )
												{
													my $cmd = qq[bsub -J split.$libID.$laneID -q $lsf_queue -o $lane_dir/split/split.o -e $lane_dir/split/split.e perl -w -e "use AssemblyTools;AssemblyTools::splitUnpairedFastq( '$lane_read0', 'split', '$BASES_PER_CHUNK', '$lane_dir/split' );"];
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
													if( length( `grep "clip1" ../meta.info | wc -l`) > 1 )
													{
														$clipPoint1=`grep "clip1" ../meta.info | awk -F ":" '{print \$2}' | head -1`;chomp( $clipPoint1 );if( length( $clipPoint1 ) == 0 ){$clipPoint1=-1;}
														$clipPoint2=`grep "clip2" ../meta.info | awk -F ":" '{print \$2}' | head -1`;chomp( $clipPoint2 );if( length( $clipPoint2 ) == 0 ){$clipPoint2=-1;}
													}
													
													for( my $s=0;$s<$numSplits;$s++ )
													{
														my $read1="split$s"."_1.fastq";
														my $read2="split$s"."_2.fastq";
														my $prefix=`echo $read1 | sed 's/\(.*\)_1.fastq/\1/'`;
														chomp( $prefix );
														
														#map the split
														if( $clipPoint1 == -1 && $clipPoint2 == -1 )
														{
															if( $expectedInsert && $expectedInsert > 1000 )
															{
																my $insert = $expectedInsert*2;
																my $cmd = qq[bsub -w "done(split.$libID.$laneID)" -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -u $prefix.unmapped -a 2000 -A $insert $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																print "$cmd\n";
																#system( $cmd );
															}
															else
															{
																my $cmd = qq[bsub -w "done(split.$libID.$laneID)" -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -u $prefix.unmapped -a 2000 $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																print "$cmd\n";
																#system( $cmd );
															}
														}
														elsif( $clipPoint1 == -1 && $clipPoint2 > -1 )
														{
															if( $expectedInsert && $expectedInsert > 1000 )
															{
																my $insert = $expectedInsert*2;
																my $cmd = qq[bsub -w "done(split.$libID.$laneID)" -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -2 $clipPoint2 -u $prefix.unmapped -a 2000 -A $insert $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																print "$cmd\n";
																#system( $cmd );
															}
															else
															{
																my $cmd = qq[bsub -w "done(split.$libID.$laneID)" -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -2 $clipPoint2 -u $prefix.unmapped -a 2000 $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																print "$cmd\n";
																#system( $cmd );
															}
														}
														elsif( $clipPoint1 > -1 && $clipPoint2 == -1 )
														{
															if( $expectedInsert && $expectedInsert > 1000 )
															{
																my $insert = $expectedInsert*2;
																my $cmd = qq[bsub -w "done(split.$libID.$laneID)" -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -1 $clipPoint1  -u $prefix.unmapped -a 2000 -A $insert $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																print "$cmd\n";
																#system( $cmd );
															}
															else
															{
																my $cmd = qq[bsub -w "done(split.$libID.$laneID)" -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -1 $clipPoint1  -u $prefix.unmapped -a 2000 $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																print "$cmd\n";
																#system( $cmd );
															}
														}
														else
														{
															if( $expectedInsert && $expectedInsert > 1000 )
															{
																my $insert = $expectedInsert*2;
																my $cmd = qq[bsub -w "done(split.$libID.$laneID)" -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -1 $clipPoint1 -2 $clipPoint2 -u $prefix.unmapped -a 2000 -A $insert $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																print "$cmd\n";
																#system( $cmd );
															}
															else
															{
																my $cmd = qq[bsub -w "done(split.$libID.$laneID)" -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -1 $clipPoint1 -2 $clipPoint2 -u $prefix.unmapped -a 2000 $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																print "$cmd\n";
																#system( $cmd );
															}
														}
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
														my $prefix=`echo $read0 | sed 's/\(.*\).fastq/\1/'`;
														chomp( $prefix );
														
														#map the split
														if( $clipPoint0 == -1 )
														{
															my $cmd = qq[bsub -w "done(split.$libID.$laneID)" -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read0 $read0.bfq; $MAQ match -u $prefix.unmapped $prefix.raw.map ../ref.bfa $read0.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read0.bfq"];
															print "$cmd\n";
															#system( $cmd );
														}
														else
														{
															my $cmd = qq[bsub -w "done(split.$libID.$laneID)" -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read0 $read0.bfq; $MAQ match -1 $clipPoint0 -u $prefix.unmapped $prefix.raw.map ../ref.bfa $read0.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read0.bfq"];
															print "$cmd\n";
															#system( $cmd );
														}
														$mapping=1;
													}
												}
											}
											elsif( ! -f "raw.map" && -d "split" )
											{
												#system( "echo 'touched' > lane.touch" );
												
												#the lane splits are partially mapped
												chdir( "split" );
												if( $lane_read1 && $lane_read2 )
												{
													foreach my $read1 (`ls split*_1.fastq`)
													{
														if( ! -f $read1 )
														{
															last;
														}

														my $prefix=`echo $read1 | sed 's/\(.*\)_1.fastq/\1/'`;
														chomp( $prefix );
														my $read2=$prefix."_2.fastq";
														
														#check if the clip points exist - else dont clip
														my $clipPoint1=-1;
														my $clipPoint2=-1;
														if( length( `grep "clip1" ../meta.info | wc -l` ) > 1 )
														{
															$clipPoint1=`grep "clip1" ../meta.info | awk -F ":" '{print \$2}' | head -1`;chomp( $clipPoint1 );if( length( $clipPoint1 ) == 0 ){$clipPoint1=-1;}
															$clipPoint2=`grep "clip2" ../meta.info | awk -F ":" '{print \$2}' | head -1`;chomp( $clipPoint2 );if( length( $clipPoint2 ) == 0 ){$clipPoint2=-1;}
														}
														
														#map the split
														if( ( -f "raw.map" && -s "raw.map" < 2000 ) || ( ! -f "raw.map" && ! -d "split" ) )
														{
															if( $clipPoint1 == -1 && $clipPoint2 == -1 ) 
															{
																if( $expectedInsert && $expectedInsert > 1000 )
																{
																	my $insert = $expectedInsert*2;
																	my $cmd = qq[bsub -J map.$libID.$laneID.$read2 -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -u $prefix.unmapped -a 2000 -A $insert $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																	print "$cmd\n";
																	#system( $cmd );
																}
																else
																{
																	my $cmd = qq[bsub -J map.$libID.$laneID.$read2 -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -u $prefix.unmapped -a 2000 $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																	print "$cmd\n";
																	#system( $cmd );
																}
															}
															elsif( $clipPoint1 == -1 && $clipPoint2 > -1 )
															{
																if( $expectedInsert && $expectedInsert > 1000 )
																{
																	my $insert = $expectedInsert*2;
																	my $cmd = qq[bsub -J map.$libID.$laneID.$read2 -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -2 $clipPoint2 -u $prefix.unmapped -a 2000 -A $insert $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																	print "$cmd\n";
																	#system( $cmd );
																}
																else
																{
																	my $cmd = qq[bsub -J map.$libID.$laneID.$read2 -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -2 $clipPoint2 -u $prefix.unmapped -a 2000 $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																	print "$cmd\n";
																	#system( $cmd );
																}
															}
															elsif( $clipPoint1 > -1 && $clipPoint2 == -1 )
															{
																if( $expectedInsert && $expectedInsert > 1000 )
																{
																	my $insert = $expectedInsert*2;
																	my $cmd = qq[bsub -J map.$libID.$laneID.$read2 -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -1 $clipPoint1  -u $prefix.unmapped -a 2000 -A $insert $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																	print "$cmd\n";
																	#system( $cmd );
																}
																else
																{
																	my $cmd = qq[bsub -J map.$libID.$laneID.$read2 -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -1 $clipPoint1  -u $prefix.unmapped -a 2000 $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																	print "$cmd\n";
																	#system( $cmd );
																}
															}
															else
															{
																if( $expectedInsert && $expectedInsert > 1000 )
																{
																	my $insert = $expectedInsert*2;
																	my $cmd = qq[bsub -J map.$libID.$laneID.$read2 -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -1 $clipPoint1 -2 $clipPoint2 -u $prefix.unmapped -a 2000 -A $insert $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																	print "$cmd\n";
																	#system( $cmd );
																}
																else
																{
																	my $cmd = qq[bsub -J map.$libID.$laneID.$read2 -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read1 $read1.bfq; $MAQ fastq2bfq $read2 $read2.bfq; $MAQ match -1 $clipPoint1 -2 $clipPoint2 -u $prefix.unmapped -a 2000 $prefix.raw.map ../ref.bfa $read1.bfq $read2.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read1.bfq $read2.bfq"];
																	print "$cmd\n";
																	#system( $cmd );
																}
															}
															$mapping=1;
														}
													}
												}	
												elsif( $lane_read0 )
												{
													my $s = 0;
													foreach my $read0 (`ls split*.fastq`)
													{
														if( ! -f $read0 )
														{
															last;
														}
														
														my $prefix=`echo $read0 | sed 's/\(.*\).fastq/\1/'`;chomp( $prefix );
														my $clipPoint0=-1;
														
														if( length( `grep "clip0" ../meta.info | wc -l` ) > 1 )
														{
															$clipPoint0=`grep clip0 ../meta.info | awk -F ":" '{print \$2}' | head -1`;chomp( $clipPoint0 );if( length( $clipPoint0 ) == 0 ){$clipPoint0=-1;}
														}
														
														if( ! -f "$prefix.raw.map" || ( -f "$prefix.raw.map" && -s "$prefix.raw.map" < 2000 ) )
														{
															#map the split
															if( $clipPoint0 == -1 )
															{
																my $cmd = qq[bsub -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read0 $read0.bfq; $MAQ match -u $prefix.unmapped $prefix.raw.map ../ref.bfa $read0.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read0.bfq"];
																print "$cmd\n";
																#system( $cmd );
															}
															else
															{
																my $cmd = qq[bsub -J map.$libID.$laneID.$s -q $lsf_queue -o ../map.o -e ../map.e "rm ../$prefix.maq.out.gz; $MAQ fastq2bfq $read0 $read0.bfq; $MAQ match -1 $clipPoint0 -u $prefix.unmapped $prefix.raw.map ../ref.bfa $read0.bfq 2>&1 | gzip -c > ../$prefix.maq.out.gz; rm $read0.bfq"];
																print "$cmd\n";
																#system( $cmd );
															}
														}
														$mapping=1;
													}
												}
											}
											
											chdir( ".." );
											#check all the map files are there
											if( $mapping == 1 )
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
												
												print S qq[												
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
													rm -rf split
												fi];
												close( S );
												
												my $cmd = qq[bsub -J merge.$libID.$laneID -q $lsf_queue -o $lane_dir/merge.o -e $lane_dir/merge.e -w "done(map.$libID.$laneID.*)" "sh check_and_merge.sh; rm check_and_merge.sh; rm lane.touch"];
												print "$cmd\n";
												#system( $cmd );
											}
											
											chdir( ".." );
											if( length( $indexF ) > 0 )
											{
												laneToBAM( getcwd(), $libID, $laneID, $lane_read0, $lane_read1, $lane_read2, basename( getcwd() ), $gender, $indexF );
											}
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

sub laneToBAM
{
	croak "Usage: laneToBAM cwd libID laneID read0 read1 read2 laneAccession gender index_file\n" unless @_ == 9;
	my $cwd = shift;
	my $libID = shift;
	my $laneID = shift;
	my $lane_read0 = shift;
	my $lane_read1 = shift;
	my $lane_read2 = shift;
	my $accession = shift;
	my $gender = shift;
	my $indexF = shift;
	
	croak "Cant find lane directory: $cwd\n" unless -d $cwd;
	croak "Cant find index file: $indexF\n" unless -f $indexF;
	#unlink( "rmdup.bam" )  unless ! -f "rmdup.bam";
	chdir( $cwd );
	if( -s "rmdup.bam" )
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
	elsif( ! -f "raw.map" )
	{
		print "No jobs running and no raw.map file\n";
		return;
	}
	
	my $tmp_directory = "/tmp/".int( rand( 1000000) );
	
	#write the SAM header
	open( H, ">raw.sam" ) or die "Cannot create raw.sam: $!\n";
	print H "\@HD\tVN:1.0\n";
	
	#write the SQ headers
	if( $gender eq 'male' || $gender eq 'unknown' )
	{
		open( FAI, $MALE_REF_FAI ) or die $!;
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
		open( FAI, $FEMALE_REF_FAI ) or die $!;
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
		print H "\tPI:$insertSize\n";
	}
	else
	{
		print H "\n";
	}
	close( H );
	
	system( "echo 'touched' > lane.touch" );
	
	unlink( "bam.o" ) unless ! -f "bam.o";
	unlink( "bam.e" ) unless ! -f "bam.e";
	
	#convert to bam
	my $cmd = qq[bsub -R "select[type==X86_64] rusage[tmp=15000]" -J bam.$libID.$laneID -q normal $jobCondition -o bam.o -e bam.e "mkdir $tmp_directory;maq2sam-long raw.map $library > $tmp_directory/tmp.sam; cat raw.sam $tmp_directory/tmp.sam > $tmp_directory/raw.sam; rm $tmp_directory/tmp.sam;rm raw.sam];
	if( $gender eq 'male' || $gender eq 'unknown' )
	{
		$cmd .= qq[; $SAMTOOLS import $MALE_REF_FAI $tmp_directory/raw.sam - | $SAMTOOLS sort - $tmp_directory/raw.sorted; $SAMTOOLS rmdup $tmp_directory/raw.sorted.bam $cwd/rmdup.bam; rm -rf $tmp_directory;rm lane.touch"]
	}
	else
	{
		$cmd .= qq[; $SAMTOOLS import $FEMALE_REF_FAI $tmp_directory/raw.sam - | $SAMTOOLS sort - $tmp_directory/raw.sorted; $SAMTOOLS rmdup $tmp_directory/raw.sorted.bam $cwd/rmdup.bam; rm -rf $tmp_directory;rm lane.touch"];
	}
	
	#print $cmd."\n";
	system( $cmd );
}

sub symLinkReference
{
	croak "Usage: symLinkReference ref_fa ref_bfa" unless @_ == 2;
	
	my $ref_fa = shift;
	my $ref_bfa = shift;
	
	if( ! -l 'ref.fa' )
	{
		symlink( $ref_fa, 'ref.fa' );
		symlink( $ref_bfa, 'ref.bfa' );
	}
}

1;
