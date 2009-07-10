package Recalibration;
use strict;
use Carp;
use Cwd;
use POSIX qw( floor );
use File::Basename;
use File::Copy;
use Data::Dumper;

use AssemblyTools;

#a perl module of many commonly used for recalibrating g1k sequencing data

sub cigar2MatchMismatchTable
{
	croak "Usage: cigar2MatchMismatch reads_fastq reference_fasta cigar dbSNP_pos output_file\n" unless @_ == 5;
	my $fastq = shift;
	my $reference = shift;
	my $cigar = shift;
	my $ignorePosFile = shift;
	my $outputFile = shift;
	
	croak "Cannot find file: $fastq\n" unless -f $fastq;
	croak "Cannot find file: $reference\n" unless -f $reference;
	croak "Cannot find file: $cigar\n" unless -f $cigar;
	croak "Cannot find file: $ignorePosFile\n" unless -f $ignorePosFile;
	
	my %ignorePos;
	open( IN, "$ignorePosFile" ) or die "Cannot open $ignorePosFile\n";
	while( <IN> )
	{
		chomp;
		#assume 1st 2 columns are chr and pos
		my @s = split( /\s+/, $_ );
		$ignorePos{ $s[ 0 ].'_'.$s[ 1 ] } = 1;
	}
	close( IN );
	
	print "Read in ignore positions\n";
	
	my %reads_hash = %{ AssemblyTools::createFastqHash( $fastq ) };
	
	print "Read in reads\n";
	
	my $ref = AssemblyTools::hash_ssaha_cigar_output( $cigar );
	print "Hashed cigar\n";
	my %cigarHash = %{ $ref };
	
	#remove multi hit reads
	foreach( keys( %cigarHash ) )
	{
		my @hits = @{ $cigarHash{ $_ } };
		if( @hits > 1 )
		{
			delete( $cigarHash{ $_ } );
		}
	}
	
	my $currentChrSeq = '';
	my $currentChrName = '';
	my %mismatch;
	open( REF, "$reference" ) or die "Cannot open reference file\n";
	while( <REF> )
	{
		chomp;
		if( $_ =~ /^>/ )
		{
			if( length( $currentChrName ) > 0 )
			{
				$ref = updateMismatchTable( \%cigarHash, \%reads_hash, \$currentChrSeq, \$currentChrName, \%mismatch, \%ignorePos );
				%mismatch = %{ $ref };
				$currentChrName = $_;
			}
			
			$currentChrSeq = '';
			my $t = ( split( /\s+/, $_ ) )[ 0 ];
			$currentChrName = substr( $t, 1 );
		}
		else
		{
			$currentChrSeq .= $_;
		}
	}
	close( REF );
	
	$ref = updateMismatchTable( \%cigarHash, \%reads_hash, \$currentChrSeq, \$currentChrName, \%mismatch, \%ignorePos );
	%mismatch = %{ $ref };
	
	open( OUT, ">$outputFile" ) or die "Cannot create output file: $outputFile in ".getcwd()."\n";
	print OUT "Q\tMatch\tMismatch\n";
	foreach( sort( keys( %mismatch ) ) )
	{
		my @t = @{ $mismatch{ $_ } };
		print OUT "$_\t$t[ 0 ]\t$t[ 1 ]\n";
	}
	close( OUT );
}

sub updateMismatchTable
{
	croak "Usage: updateMismatchTable cigar_ref_hash reads_ref_hash seq_hash chr_name_ref mismatch_ref ignorePos_ref\n" unless @_ == 6;
	my $cigar_ref = shift;
	my $reads_ref = shift;
	my $seq_ref = shift;
	my $chr_name_ref = shift;
	my $mismatch_ref = shift;
	my $ignore_ref = shift;
	
	my %cigar = %{ $cigar_ref };
	my %reads_hash = %{ $reads_ref };
	my $seq = ${ $seq_ref };
	my $chr_name = ${ $chr_name_ref };
	my %mismatch = %{ $mismatch_ref };
	my %ignorePos = %{ $ignore_ref };
	
	foreach( keys( %cigar ) )
	{
		my @hits = @{ $cigar{ $_ } };
		
		#assumes there's only 1 hit per read so just take first entry
		my @splits = split( /\s+/, $hits[ 0 ] );
		if( $chr_name eq $splits[ 5 ] )
		{
			croak "Cannot find read in fastq file: $splits[ 1 ] in ".getcwd()."\n" unless defined( $reads_hash{ $splits[ 1 ] } );
			
			#use this read alignment
			my $contig_index = $splits[ 6 ] - 1; #ssaha starts counting from 1
			my $read_index = $splits[ 2 ] - 1;
			if( $splits[ 4 ] eq '+' ) #if forward match
			{
				for( my $i = 10; $i < @splits; $i = $i + 2 )
				{
					if( $splits[ $i ] eq 'M' )
					{
						my $match_count = 0; #exact matching positions counter
						while( $match_count < $splits[ $i + 1 ] )
						{
							if( ! defined( $ignorePos{ $splits[ 5 ].'_'.( $contig_index + 1 ) } ) )
							{
								my @readSplits = split( /\n/, $reads_hash{ $splits[ 1 ] } );
								my $read_base = uc( substr( $readSplits[ 1 ], $read_index, 1 ) );
								my $qual = (unpack("C", substr( $readSplits[ 3 ], $read_index, 1 ) ) - 33);
								
								my $contig_base = substr( $seq, $contig_index, 1 );
								
								if( ! defined $mismatch{ $qual } )
								{
									$mismatch{ $qual } = [ 0, 0 ];
								}
								
								if( $contig_base eq $read_base )
								{
									$mismatch{ $qual }[ 0 ] ++;
								}
								else
								{
									$mismatch{ $qual }[ 1 ] ++;
									$mismatch{ $qual }[ 0 ] ++;
								}
							}
							$contig_index ++;
							$read_index ++;
							$match_count ++;
						}
					}
					elsif( $splits[ $i ] eq 'I' )
					{
						$read_index += $splits[ $i + 1 ];
					}
					elsif( $splits[ $i ] eq 'D' )
					{
						$contig_index += $splits[ $i + 1 ];
					}
				}
			}
			else #reverse alignment
			{
				for( my $i = 10; $i < @splits; $i = $i + 2 )
				{
					if( $splits[ $i ] eq 'M' )
					{
						my $match_count = 0; #exact matching positions counter
						while( $match_count < $splits[ $i + 1 ] )
						{
							if( ! defined( $ignorePos{ $splits[ 5 ].'_'.( $contig_index + 1 ) } ) )
							{
								my @readSplits = split( /\n/, $reads_hash{ $splits[ 1 ] } );
								my $read_base = uc( complement( substr( $readSplits[ 1 ], $read_index, 1 ) ) );
								my $qual = (unpack("C", substr( $readSplits[ 3 ], $read_index, 1 ) ) - 33);
								
								my $contig_base = substr( $seq, $contig_index, 1 );
								
								if( ! defined $mismatch{ $qual } )
								{
									$mismatch{ $qual } = [ 0, 0 ];
								}
								
								if( $contig_base eq $read_base )
								{
									$mismatch{ $qual }[ 0 ] ++;
								}
								else
								{
									$mismatch{ $qual }[ 1 ] ++;
									$mismatch{ $qual }[ 0 ] ++;
								}
							}
							$contig_index ++;
							$read_index --;
							$match_count ++;
						}
					}
					elsif( $splits[ $i ] eq 'I' )
					{
						$read_index -= $splits[ $i + 1 ];
					}
					elsif( $splits[ $i ] eq 'D' )
					{
						$contig_index += $splits[ $i + 1 ];
					}
				}
			}
			delete( $cigar{ $_ } );
		}
	}
	
	return \%mismatch;
}

#assumes you are giving it a cigar file with 1 hit per read!
sub cigar2MatchMismatchTableSingleChr
{
	croak "Usage: cigar2MatchMismatch reads_fastq reference_fasta cigar dbSNP_pos output_file\n" unless @_ == 5;
	my $fastq = shift;
	my $reference = shift;
	my $cigar = shift;
	my $ignorePosFile = shift;
	my $outputFile = shift;
	
	croak "Cannot find file: $fastq\n" unless -f $fastq;
	croak "Cannot find file: $reference\n" unless -f $reference;
	croak "Cannot find file: $cigar\n" unless -f $cigar;
	croak "Cannot find file: $ignorePosFile\n" unless -f $ignorePosFile;
	
	my %ignorePos;
	open( IN, "$ignorePosFile" ) or die "Cannot open $ignorePosFile\n";
	while( <IN> )
	{
		chomp;
		#assume 1st 2 columns are chr and pos
		my @s = split( /\s+/, $_ );
		$ignorePos{ $s[ 0 ].'_'.$s[ 1 ] } = 1;
	}
	close( IN );
	
	print "Read in ignore positions\n";
	
	#read in reference
	my %reference = %{ AssemblyTools::fastaHash( $reference ) };
	
	print 'REF: '.keys( %reference );
	
	print "Read in reference\n";
	
	my %reads_hash = %{ AssemblyTools::createFastqHash( $fastq ) };
	
	print "Read in reads\n";
	
	my %matchMismatch; #contains pointers to arrays of match/mismatch counts
	my $cigarCount = 0;
	my $validHits = 0;
	
	if( $cigar =~ /\.gz$/ )
	{
		open( CIGAR, "gunzip -c $cigar |" ) or die "Cannot open gzipped cigar file\n";
	}
	else
	{
		open( CIGAR, $cigar ) or die "cannot open cigar file\n";
	}
	while( <CIGAR> )
	{
		chomp;
		if( $_ =~ /^cigar::/ )
		{
			$cigarCount ++;
			my @splits = split( /\s+/, $_ );
			if( defined $reference{ $splits[ 5 ] } )
			{
				$validHits ++;
				croak "Cannot find read in fastq file: $splits[ 1 ] in ".getcwd()."\n" unless defined( $reads_hash{ $splits[ 1 ] } );
				
				#use this read alignment
				
				my $contig_index = $splits[ 6 ] - 1; #ssaha starts counting from 1
				my $read_index = $splits[ 2 ] - 1;
				if( $splits[ 4 ] eq '+' ) #if forward match
				{
					for( my $i = 10; $i < @splits; $i = $i + 2 )
					{
						if( $splits[ $i ] eq 'M' )
						{
							my $match_count = 0; #exact matching positions counter
							while( $match_count < $splits[ $i + 1 ] )
							{
								if( ! defined( $ignorePos{ $splits[ 5 ].'_'.( $contig_index + 1 ) } ) )
								{
									my @readSplits = split( /\n/, $reads_hash{ $splits[ 1 ] } );
									my $read_base = uc( substr( $readSplits[ 1 ], $read_index, 1 ) );
									my $qual = (unpack("C", substr( $readSplits[ 3 ], $read_index, 1 ) ) - 33);
									
									my $contig_base = substr( $reference{ $splits[ 5 ] }, $contig_index, 1 );
									
									if( ! defined $matchMismatch{ $qual } )
									{
										$matchMismatch{ $qual } = [ 0, 0 ];
									}
									
									if( $contig_base eq $read_base )
									{
										#print "MA1 $readSplits[0] $read_index $read_base $contig_index $contig_base\n";
										$matchMismatch{ $qual }[ 0 ] ++;
									}
									else
									{
										#print "MM1 $readSplits[0] $read_index $read_base $contig_index $contig_base\n";
										$matchMismatch{ $qual }[ 1 ] ++;
										$matchMismatch{ $qual }[ 0 ] ++;
									}
								}
								$contig_index ++;
								$read_index ++;
								$match_count ++;
							}
						}
						elsif( $splits[ $i ] eq 'I' )
						{
							$read_index += $splits[ $i + 1 ];
						}
						elsif( $splits[ $i ] eq 'D' )
						{
							$contig_index += $splits[ $i + 1 ];
						}
					}
				}
				else #reverse alignment
				{
					for( my $i = 10; $i < @splits; $i = $i + 2 )
					{
						if( $splits[ $i ] eq 'M' )
						{
							my $match_count = 0; #exact matching positions counter
							while( $match_count < $splits[ $i + 1 ] )
							{
								if( ! defined( $ignorePos{ $splits[ 5 ].'_'.( $contig_index + 1 ) } ) )
								{
									my @readSplits = split( /\n/, $reads_hash{ $splits[ 1 ] } );
									my $read_base = uc( complement( substr( $readSplits[ 1 ], $read_index, 1 ) ) );
									my $qual = (unpack("C", substr( $readSplits[ 3 ], $read_index, 1 ) ) - 33);
									
									my $contig_base = substr( $reference{ $splits[ 5 ] }, $contig_index, 1 );
									
									if( ! defined $matchMismatch{ $qual } )
									{
										$matchMismatch{ $qual } = [ 0, 0 ];
									}
									
									if( $contig_base eq $read_base )
									{
										#print "MA2 $readSplits[0] $read_index $read_base $contig_index $contig_base\n";
										$matchMismatch{ $qual }[ 0 ] ++;
									}
									else
									{
										#print "MM2 $readSplits[0] $read_index $read_base $contig_index $contig_base\n";
										$matchMismatch{ $qual }[ 1 ] ++;
										$matchMismatch{ $qual }[ 0 ] ++;
									}
								}
								$contig_index ++;
								$read_index --;
								$match_count ++;
							}
						}
						elsif( $splits[ $i ] eq 'I' )
						{
							$read_index -= $splits[ $i + 1 ];
						}
						elsif( $splits[ $i ] eq 'D' )
						{
							$contig_index += $splits[ $i + 1 ];
						}
					}
				}
			}
		}
	}
	close( CIGAR );
	
	print "Found total of $cigarCount in cigar file\n";
	print "Found total of $validHits hits to reference in cigar file\n";
	print "Got ".scalar( keys( %matchMismatch ) )." entries in match/mismatch table\n";
	
	open( OUT, ">$outputFile" ) or die "Cannot create output file: $outputFile in ".getcwd()."\n";
	print OUT "Q\tMatch\tMismatch\n";
	foreach( sort( keys( %matchMismatch ) ) )
	{
		my @t = @{ $matchMismatch{ $_ } };
		print OUT "$_\t$t[ 0 ]\t$t[ 1 ]\n";
	}
	close( OUT );
}

sub mergeMatchMismatchTables
{
	croak "Usage: mergeMatchMismatchTables outputTableFile directory tables_prefix" unless @_ == 3;
	my $outputFile = shift;
	my $directory = shift;
	my $prefix = shift;
	
	my %mergedMatchMismatch;
	opendir( DIR, $directory ) or die "Cannot open directory\n";
	while( (my $filename = readdir( DIR ) ) )
	{
		if( $filename =~ /^$prefix/ )
		{
			open( TABLE, "$filename" ) or die "Cannot open match/mismatch table: $filename\n";
			<TABLE>;
			while( <TABLE> )
			{
				chomp;
				my @s = split( /\s+/, $_ );
				if( ! defined( $mergedMatchMismatch{ $s[ 0 ] } ) )
				{
					$mergedMatchMismatch{ $s[ 0 ] } = [ $s[ 1 ], $s[ 2 ] ];
				}
				else
				{
					$mergedMatchMismatch{ $s[ 0 ] }[ 0 ] += $s[ 1 ];
					$mergedMatchMismatch{ $s[ 0 ] }[ 1 ] += $s[ 2 ];
				}
			}
			close( TABLE );
		}
	}
	close( DIR );
	
	open( OUT, ">$outputFile" ) or die "Cannot create output file: $outputFile in ".getcwd()."\n";
	foreach( sort( keys( %mergedMatchMismatch ) ) )
	{
		my @t = @{ $mergedMatchMismatch{ $_ } };
		print OUT "$_ $t[ 0 ] $t[ 1 ]\n";
	}
	close( OUT );
}

#a function to sample a subset of bases from a solid lane
#works for paired or unpaired lanes
#produces *.sample.fastq files
sub sampleSubsetBasesLaneSolid
{
	croak "Usage: sampleSubsetBases_solid outputDirectory sampleSizeBp fastq1 output1 [fastq2 output2]" unless @_ == 4 || @_ == 6;
	my $outDirectory = shift;
	my $sampleSize = shift;
	my $fastq1 = shift;
	my $output1 = shift;
	
	croak "Cant find fastq1 file: $fastq1 in ".getcwd()."\n" unless -f $fastq1;
	
	my $fastq2 = '';
	my $output2 = '';
	if( @_ == 2 )
	{
		$fastq2 = shift;
		$output2 = shift;
		
		croak "Cant find fastq2 file: $fastq2 in ".getcwd()."\n" unless -f $fastq2;
	}
	
	my $numReads = ( $fastq1 =~ /\.fastq\.gz$/ ? `zcat $fastq1 | wc -l` : `cat $fastq1 | wc -l` ) / 4;
	my $readLength = $fastq1 =~ /\.fastq\.gz$/ ? `zcat $fastq1 | head -2 | tail -1 | wc -m` : `head -2 $fastq1 | tail -1 | wc -m`;
	
	if( length( $fastq2 ) > 0 ){$readLength *= 2;}
	
	my $numReadsToSample = int( $sampleSize / $readLength ) + 1;
	my $readChunkSize = int( $numReads / $numReadsToSample );
	
	print "Total Reads: $numReads\n";
	print "Read Length: $readLength\n";
	
	print "Chunk Size: $readChunkSize\n";
	
	if( $readChunkSize > $numReads )
	{
		#not enough reads => whole file is the sample
		copy( $fastq1, $output1 );
		if( length( $fastq2 ) > 0 )
		{
			copy( $fastq2, $output2 );
		}
	}
	
	my %currentReads1;
	my %currentReads2;
	if( $fastq1 =~ /\.gz$/ )
	{
		open( FASTQ1, "gunzip -c $fastq1 |" ) or die "Cannot open gzipped fastq: $fastq1\n";
	}
	else
	{
		open( FASTQ1, "$fastq1" ) or die "Cannot open fastq: $fastq1\n";
	}
	
	if( length( $fastq2 ) > 0 )
	{
		if( $fastq2 =~ /\.gz$/ )
		{
			open( FASTQ2, "gunzip -c $fastq2 |" ) or die "Cannot open gzipped fastq: $fastq2\n";
		}
		else
		{
			open( FASTQ2, "$fastq2" ) or die "Cannot open fastq: $fastq2\n";
		}
	}
	
	open( OUT1, ">$output1" ) or die "Cannot create output file: $output1\n";
	
	if( length( $fastq2 ) > 0 )
	{
		open( OUT2, ">$output2" ) or die "Cannot create output file: $output2\n";
	}
	
	my $numReadsInChunk = 0;
	my $totalReadsSampled = 0;
	my $totalNReads = 0;
	my $totalReadCounter = 0;
	while( <FASTQ1> )
	{
		chomp;
		
		#print $totalReadsSampled.' '.$sampleSize.' '.$numReadsInChunk.' '.$readChunkSize."\n";
		last unless $totalReadsSampled < $numReadsToSample + 1;
		
		my $rn1 = $_;
		my $seq = <FASTQ1>;
		chomp( $seq );
		
		#dont include reads with all N's
		if( $seq =~ /^T\.+$/ || $seq =~ /^G\.+$/ )
		{
			<FASTQ1>;
			<FASTQ1>;
			
			if( length( $fastq2 ) > 0 )
			{
				<FASTQ2>;
				<FASTQ2>;
				<FASTQ2>;
				<FASTQ2>;
			}
			$totalNReads ++;
			next;
		}
				
		my $qn = <FASTQ1>;
		chomp( $qn );
		my $quals = <FASTQ1>;
		chomp( $quals );
				
		$currentReads1{ substr( $rn1, 1 ) } = '@'.$totalReadCounter."/1\n".$seq."\n".substr( $qn, 0, 1)."\n".$quals;
		
		if( length( $fastq2 ) > 0 )
		{
			my $rn2 = <FASTQ2>;
			chomp( $rn2 );
			$seq = <FASTQ2>;
			
			#dont include reads with all N's
			if( $seq =~ /^T\.+$/ || $seq =~ /^G\.+$/ )
			{
				<FASTQ2>;
				<FASTQ2>;
				
				delete $currentReads1{ substr( $rn1, 1 ) };
				
				$totalNReads ++;
				next;
			}
			
			chomp( $seq );
			
			$qn = <FASTQ2>;
			chomp( $qn );
			$quals = <FASTQ2>;
			chomp( $quals );
			
			$currentReads2{ substr( $rn2, 1 ) } = '@'.$totalReadCounter."/2\n".$seq."\n".substr( $qn, 0, 1)."\n".$quals;
			$totalReadCounter ++;
		}
		
		$numReadsInChunk ++;
		
		if( $numReadsInChunk >= $readChunkSize )
		{
			#pick a random read
			my $randomKey = (keys( %currentReads1 ))[ int( rand( keys( %currentReads1 ) ) ) ];
			
			my $read1 = $currentReads1{ $randomKey };
			print OUT1 $read1."\n";
			
			if( length( $fastq2 ) > 0 )
			{
				my $read2 = $currentReads2{ $randomKey };
				print OUT2 $read2."\n";
			}
			
			$numReadsInChunk = 0;
			%currentReads1 = ();
			%currentReads2 = ();
			$totalReadsSampled ++;
		}
	}
	
	if( keys( %currentReads1 ) > 0 )
	{
	#pick a random read
			my $randomKey = (keys( %currentReads1 ))[ int( rand( keys( %currentReads1 ) ) ) ];
			
			my $read1 = $currentReads1{ $randomKey };
			print OUT1 $read1."\n";
			
			if( length( $fastq2 ) > 0 )
			{
				my $read2 = $currentReads2{ $randomKey };
				print OUT2 $read2."\n";
			}
			
			$numReadsInChunk = 0;
			%currentReads1 = ();
			%currentReads2 = ();
			$totalReadsSampled ++;
	}
	
	print "Total Sampled: $totalReadsSampled\n";
	
	close( FASTQ1 );
	close( OUT1 );
	if( length( $fastq2 ) > 0 )
	{
		close( FASTQ2 );
		close( OUT2 );
	}
	
	print "NReads in $fastq1: $totalNReads\n";
}

#a function to sample a subset of bases from a solexa lane
#works for paired or unpaired lanes
#produces *.sample.fastq files
my $MONOMER_CUTOFF = 0.7;
sub sampleSubsetBasesLaneSolexa
{
	croak "Usage: sampleSubsetBasesLaneSolexa outputDirectory sampleSizeBp fastq1 output1 [fastq2 output2]" unless @_ == 4 || @_ == 6;
	my $outDirectory = shift;
	my $sampleSize = shift;
	my $fastq1 = shift;
	my $output1 = shift;
	
	croak "Cant find fastq1 file: $fastq1 in ".getcwd()."\n" unless -f $fastq1;
	
	my $fastq2 = '';
	my $output2 = '';
	if( @_ == 2 )
	{
		$fastq2 = shift;
		$output2 = shift;
		
		croak "Cant find fastq2 file: $fastq2 in ".getcwd()."\n" unless -f $fastq2;
	}
	
	my $numReads = ( $fastq1 =~ /\.fastq\.gz$/ ? `zcat $fastq1 | wc -l` : `cat $fastq1 | wc -l` ) / 4;
	my $readLength1 = $fastq1 =~ /\.fastq\.gz$/ ? `zcat $fastq1 | head -2 | tail -1 | wc -m` : `head -2 $fastq1 | tail -1 | wc -m`;
	print "Read1 Length: $readLength1 bp\n";
	
	my $readLength2 = 0;
	if( length( $fastq2 ) > 0 )
	{
		$readLength2 = $fastq2 =~ /\.fastq\.gz$/ ? `zcat $fastq2 | head -2 | tail -1 | wc -m` : `head -2 $fastq2 | tail -1 | wc -m`;
		print "Read2 Length: $readLength2 bp\n";
	}
	
	#see if can work out clip points to remove any trailing Ns
	my $clip1Point = -1;
	my $clip2Point = -1;
	
	$clip1Point = AssemblyTools::determineClipPointMaq( $fastq1 );
	if( length( $fastq2 ) > 0 )
	{
		$clip2Point = AssemblyTools::determineClipPointMaq( $fastq2 );
	}
	
	my $readLength = $readLength1 + $readLength2;
	
	my $numReadsToSample = int( $sampleSize / $readLength ) + 1;
	my $readChunkSize = int( $numReads / $numReadsToSample ) - 1;
	
	print "Total Reads: $numReads\n";
	print "Paired Read Length: $readLength\n";
	
	print "Chunk Size: $readChunkSize\n";
	
	if( $readChunkSize > $numReads )
	{
		print "Not enough reads for sample!!!";
		exit;
	}
	
	my %currentReads1;
	my %currentReads2;
	if( $fastq1 =~ /\.gz$/ )
	{
		open( FASTQ1, "gunzip -c $fastq1 |" ) or die "Cannot open gzipped fastq: $fastq1\n";
	}
	else
	{
		open( FASTQ1, "$fastq1" ) or die "Cannot open fastq: $fastq1\n";
	}
	
	if( length( $fastq2 ) > 0 )
	{
		if( $fastq2 =~ /\.gz$/ )
		{
			open( FASTQ2, "gunzip -c $fastq2 |" ) or die "Cannot open gzipped fastq: $fastq2\n";
		}
		else
		{
			open( FASTQ2, "$fastq2" ) or die "Cannot open fastq: $fastq2\n";
		}
	}
	
	open( OUT1, ">$output1" ) or die "Cannot create output file: $output1\n";
	
	if( length( $fastq2 ) > 0 )
	{
		open( OUT2, ">$output2" ) or die "Cannot create output file: $output2\n";
	}
	
	my $numReadsInChunk = 0;
	my $totalReadsSampled = 0;
	my $totalNReads1 = 0;
	my $totalNReads2 = 0;
	my $totalReadCounter = 0;
	my $totalBasesSampled = 0;
	my $fileIterations = 0;
	while()
	{
		if( $totalBasesSampled > $sampleSize || $fileIterations > 10 )
		{
			last;
		}
		
		my $rn1 = <FASTQ1>;
		
		if( ! defined( $rn1 ) )
		{
			print "Sampled $totalBasesSampled bases. Target: $sampleSize bp\n";
			print "Resetting to start of file....\n";
			
			#reset to the start of the file again
			if( $fastq1 =~ /\.gz$/ )
			{
				open( FASTQ1, "gunzip -c $fastq1 |" ) or die "Cannot open gzipped fastq: $fastq1\n";
			}
			else
			{
				open( FASTQ1, "$fastq1" ) or die "Cannot open fastq: $fastq1\n";
			}
			
			if( length( $fastq2 ) > 0 )
			{
				if( $fastq2 =~ /\.gz$/ )
				{
					open( FASTQ2, "gunzip -c $fastq2 |" ) or die "Cannot open gzipped fastq: $fastq2\n";
				}
				else
				{
					open( FASTQ2, "$fastq2" ) or die "Cannot open fastq: $fastq2\n";
				}
			}
			
			$numReadsInChunk = 0;
			%currentReads1 = ();
			%currentReads2 = ();
			
			$fileIterations ++; #count the number of times iterated over the file trying to sample enough bases
			
			$rn1 = <FASTQ1>;
		}
		
		chomp( $rn1 );
		my $seq1 = <FASTQ1>;
		chomp( $seq1 );
		
		if( $clip1Point > -1 )
		{
			$seq1 = substr( $seq1, 0, $clip1Point );
		}
		
		#count the frequency of each base in the read
		$seq1 = uc( $seq1 );
		my $countA = ($seq1=~tr/A//);
		my $countC = ($seq1=~tr/C//);
		my $countG = ($seq1=~tr/G//);
		my $countT = ($seq1=~tr/T//);
		my $countN = ($seq1=~tr/N//);
		
		my $monomerCutoff = length( $seq1 ) * $MONOMER_CUTOFF;
		if( $countA > $monomerCutoff || $countC > $monomerCutoff || $countG > $monomerCutoff || $countT > $monomerCutoff || $countN > 2 )
		{
			<FASTQ1>;
			<FASTQ1>;
			
			if( length( $fastq2 ) > 0 )
			{
				<FASTQ2>;
				<FASTQ2>;
				<FASTQ2>;
				<FASTQ2>;
			}
			$totalNReads1 ++;
			$numReadsInChunk ++;
			next;
		}
				
		my $qn = <FASTQ1>;
		chomp( $qn );
		my $quals = <FASTQ1>;
		chomp( $quals );
		if( $clip1Point > -1 )
		{
			$quals = substr( $quals, 0, $clip1Point );
		}
		
		$currentReads1{ $numReadsInChunk } = '@'.$totalReadCounter."/1\n".$seq1."\n".substr( $qn, 0, 1)."\n".$quals;
		
		my $seq2 = '';
		if( length( $fastq2 ) > 0 )
		{
			my $rn2 = <FASTQ2>;
			chomp( $rn2 );
			$seq2 = <FASTQ2>;
			chomp( $seq2 );
			
			if( $clip2Point > -1 )
			{
				$seq2 = substr( $seq2, 0, $clip2Point );
			}
		
			#count the frequency of each base in the read
			$seq2 = uc( $seq2 );
			$countA = ($seq2=~tr/A//);
			$countC = ($seq2=~tr/C//);
			$countG = ($seq2=~tr/G//);
			$countT = ($seq2=~tr/T//);
			$countN = ($seq2=~tr/N//);
			
			my $monomerCutoff = length( $seq2 ) * $MONOMER_CUTOFF;
			if( $countA > $monomerCutoff || $countC > $monomerCutoff || $countG > $monomerCutoff || $countT > $monomerCutoff || $countN > 2 )
			{
				<FASTQ2>;
				<FASTQ2>;
				
				delete $currentReads1{ $numReadsInChunk };
				
				$totalNReads2 ++;
				$numReadsInChunk ++;
				next;
			}
			
			$qn = <FASTQ2>;
			chomp( $qn );
			$quals = <FASTQ2>;
			chomp( $quals );
			
			if( $clip2Point > -1 )
			{
				$quals = substr( $quals, 0, $clip2Point );
			}
			
			$currentReads2{ $numReadsInChunk } = '@'.$totalReadCounter."/2\n".$seq2."\n".substr( $qn, 0, 1)."\n".$quals;
		}
		
		$totalReadCounter ++;
		$numReadsInChunk ++;
		
		if( $numReadsInChunk >= $readChunkSize )
		{
			#pick a random read
			my $randomIndex = int( rand( keys( %currentReads1 ) ) );
			my $randomKey = (sort( keys( %currentReads1 ) ) )[ $randomIndex ];
			
			my $read1 = $currentReads1{ $randomKey };
			print OUT1 $read1."\n";
			
			if( length( $fastq2 ) > 0 )
			{
				croak "Number of reads in sample chunk 1 and 2 not equal: ".scalar(keys( %currentReads1 )).' '.scalar(keys( %currentReads2 ))."\n" unless scalar(keys( %currentReads1 )) == scalar(keys( %currentReads2 ));
				
				my $read2 = $currentReads2{ (sort( keys( %currentReads2 ) ) )[ $randomIndex ] };
				print OUT2 $read2."\n";
			}
			
			$numReadsInChunk = 0;
			%currentReads1 = ();
			%currentReads2 = ();
			$totalReadsSampled ++;
			$totalBasesSampled += length( $seq1 ) + length( $seq2 );
		}
	}
	
	close( FASTQ1 );
	close( OUT1 );
	
	print "Total reads Sampled: $totalReadsSampled\n";
	print "Total bases Sampled: $totalBasesSampled vs. $sampleSize\n";
	print "NReads in read1: $totalNReads1\n";
	
	if( length( $fastq2 ) > 0 )
	{
		print "NReads in read2: $totalNReads2\n";
		close( FASTQ2 );
		close( OUT2 );
	}
	
	#check that we sampled enough bases - otherwise make empty files to indicate this
	if( $totalBasesSampled < $sampleSize )
	{
		print "Failed to sample enough bases: $totalBasesSampled bp\n";
		unlink( $output1 );
		open( T, ">$output1" );close( T );
		
		if( length( $fastq2 ) > 0 )
		{
			unlink( $output2 );
			open( T, ">$output2" );close( T );
		}
	}
}

#a function to sample a subset of bases from a 454 fastq file (variable read lengths)
sub sampleSubsetBasesFastq
{
	croak "Usage: sampleSubsetBasesFastq fastq outputDirectory output_file sampleSizeBp" unless @_ == 4;
	my $fastq = shift;
	my $outDirectory = shift;
	my $outputFile = shift;
	my $sampleSize = shift;
	my ( $volume, $directories, $file ) = File::Spec->splitpath( $fastq );
	
	if( $fastq =~ /\.gz$/ )
	{
		open( FASTQ, "gunzip -c $fastq |" ) or die "Cannot open gzipped fastq\n";
	}
	else
	{
		open( FASTQ, "$fastq" ) or die "Cannot open fastq\n";
	}
	
	#count the number of bases in the file
	my $totalBases = 0;
	my %readNames;
	while( <FASTQ> )
	{
		chomp;
		my $name = $_;
		$readNames{ $name } = 1;
		
		my $seq = <FASTQ>;
		chomp( $seq );
		$totalBases += length( $seq );
		
		my $qname = <FASTQ>;
		chomp( $qname );
		my $quals = <FASTQ>;
		chomp( $quals );
	}
	close( FASTQ );
	
	print "$fastq has $totalBases bp\n";
	
	if( $totalBases < $sampleSize )
	{
		#the whole file is the sample!
		system( "cp $fastq $outDirectory/$outputFile" );
	}
	else
	{
		my $blockSizeBp = int( ( $totalBases / $sampleSize ) + .5);
		my %readsBlock;
		my $totalSampled = 0;
		my $totalBlockBp = 0;
		
		if( $fastq =~ /\.gz$/ )
		{
			open( FASTQ, "gunzip -c $fastq |" ) or die "Cannot open gzipped fastq\n";
		}
		else
		{
			open( FASTQ, "$fastq" ) or die "Cannot open fastq\n";
		}
		
		open( OUT, ">$outDirectory/$outputFile" ) or die "Cannot create output file $outDirectory/$outputFile\n";
		
		while( <FASTQ> )
		{
			chomp;
			my $name = $_;
			my $seq = <FASTQ>;
			chomp( $seq );
			my $qname = <FASTQ>;
			chomp( $qname );
			my $quals = <FASTQ>;
			chomp( $quals );
			
			$totalBlockBp += length( $seq );
			
			$readsBlock{ $name."\n".$seq."\n".$qname."\n".$quals } = 1;
			
			if( $totalBlockBp > $blockSizeBp )
			{
				#randomly sample a read the block
				my $read = (keys( %readsBlock ))[ int( rand( keys( %readsBlock ) ) ) ];
				print OUT $read."\n";
				
				$totalSampled += length( ( split( /\n/, $read ) )[ 1 ] );
				
				$totalBlockBp = 0;
				%readsBlock = ();
			}
			last unless $totalSampled < $sampleSize;
		}
		close( FASTQ );
		close( OUT );
		
		print "Total Sampled: $totalSampled bp\n";
	}
}

sub applyMatchMismatchTable
{
	croak "Usage: applyMatchMismatchTable originalFastq matchMismatch outputRecalibratedFastq outputSummaryTable" unless @_ == 4;
	my $original = shift;
	my $table = shift;
	my $outputFastq = shift;
	my $outputSummaryTable = shift;
	
	my %matchMismatchTable;
	open( IN, $table ) or die "cannot open match/mismatch table: $table\n";
	<IN>; #header line
	while( <IN> )
	{
		chomp;
		my @s = split( /\s+/, $_ );
		$matchMismatchTable{ $s[ 0 ] } = [ $s[ 1 ], $s[ 2 ] ];
	}
	close( IN );
	
	my %recalibrations;
	open( SUM, ">$outputSummaryTable" ) or die "Cannot create summary table\n";
	print SUM "Qval\tMatch\tmismatch\tNewQval\n";
	foreach( sort( keys( %matchMismatchTable ) ) )
	{
		my @entries = @{ $matchMismatchTable{ $_ } };
		
		if( $entries[ 1 ] > 0 )
		{
			$recalibrations{ $_ } = estimateQuality( ($entries[ 0 ] + $entries[ 1 ]), $entries[ 1 ] );
			print SUM "$_\t$entries[ 0 ]\t$entries[ 1 ]\t".$recalibrations{ $_ }."\n";
		}
	}
	close( SUM );
	
	#now apply the recalibration to the reads
	if( $original =~ /\.gz$/ )
	{
		open( IN, "gunzip -c $original |" ) or die "Cannot open gzipped fastq\n";
		open( OUT, "| gzip -c >$outputFastq" ) or die "Cannot create recalibrated fastq";
	}
	else
	{
		open( IN, "$original" ) or die "Cannot open fastq\n";
		open( OUT, ">$outputFastq" ) or die "Cannot create recalibrated fastq";
	}
	while( <IN> )
	{
		chomp;
		print OUT $_."\n";
		my $t = <IN>;
		print OUT $t;
		$t = <IN>;
		print OUT $t;
		$t = <IN>;
		chomp( $t );
		my @s = split( //, $t );
		foreach( @s )
		{
			my $qual = (unpack("C", $_ ) - 33);
			if( defined $recalibrations{ $qual } )
			{
				print OUT chr( $recalibrations{ $qual } + 33 );
			}
			else
			{
				my $nearest = $qual - 1; #find next lower recalibrated value
				while( ! defined $recalibrations{ $nearest } && $nearest > -2 )
				{
					$nearest --;
				}
				
				print OUT $nearest < 0 ? 0 : chr( $recalibrations{ $nearest } + 33 );
			}
		}
		print OUT "\n";
	}
	close( OUT );
	close( IN );
}

# directly from makeQualitiesBayesian.pl
# Arguments: n bases, of which m are (we assume) errors.
# Prior distribution is uniform over qualities 0 ... $QMAX.
my $QMAX = 40;
sub estimateQuality {
    my ($n,$m) = @_;
    # If we have no data or bad data, posterior = prior so we return its mean.
    return $QMAX/2 if ($n == 0 || $m > $n);
    # $lch is log(n Choose m).
    my $lch = 0;
    foreach my $i (0 .. $m-1) {
	$lch += log(($n-$i)/($m-$i));
    }
    my ($top,$bot);
    # For each quality in the prior range...
    foreach my $q (0 .. $QMAX) {
	# $p is error probability derived from $q. We
	# set an upper limit of 0.5 (very low qualities are
	# not very meaningful anyway).
	my $p = min(0.5,exp(-log(10)*0.1*$q));
	# Log prob (log likelihood) of m given n and p
	# is log(p^m * (1-p)^(n-m)), times n Choose m...
	my $lpm = $m*log($p)+($n-$m)*log(1-$p);
	my $lp = $lpm+$lch;
	# Posterior prob is exp($lp) (times the prior which is
	# constant and can be ignored)
	my $post = exp($lp);
	# Work towards estimating the posterior mean...
	$top += $post*$q;
	$bot += $post;
    }
    # If we didn't gather any data, return middle of range.
    return $QMAX/2 if ($bot == 0);
    # Find nearest integer to posterior mean $top/$bot.
    my $q = int(0.5+$top/$bot);
    return $q;
}

sub min
{
	my $a = shift;
	my $b = shift;
	return $a < $b ? $a : $b;
}

sub log10 
{
	my $n = shift;
	return log($n)/log(10);
}

=pod
A function to apply the per base position bayesian qual map
=cut
sub applyPosQualMap
{
	croak "Usage: applyPosQualMap qualmap destDir originalFastq1 recalFastq1 [originalFastq2 recalFastq2 originalUnpaired recalUnpaired]" unless @_ == 4 || @_ == 6 || @_ == 8;
	my $qualmap = shift;
	my $dest = shift;
	my $original1 = shift;
	my $newFastq1 = shift;
	croak "Cant find fastq1" unless -f $original1;
	
	my $original2 = '';
	my $newFastq2 = '';
	if( @_ >= 2 )
	{
		$original2 = shift;
		$newFastq2 = shift;
		croak "Cant find fastq2" unless -f $original2;
	}
	
	my $originalUnpaired = '';
	my $recalUnpaired = '';
	if( @_ == 2 )
	{
		$originalUnpaired = shift;
		$recalUnpaired = shift;
		croak "Cant find originalUnpaired" unless -f $originalUnpaired;
	}
	
	croak "Destination not a directory!\n" unless -d $dest;
	my %recal1;
	my %recal2;
	
	open( QUALMAP, $qualmap ) or die "Cannot open qualmap: $!\n";
	while( <QUALMAP> ) 
	{
		chomp;
		my @s = split( /\s+/, $_ );
		if( $s[ 0 ] == 0 )
		{
			$recal1{ $s[ 1 ].'_'.$s[ 2 ] } = $s[ 3 ];
		}
		if( $s[ 0 ] == 1 )
		{
			$recal1{ $s[ 1 ].'_'.$s[ 2 ] } = $s[ 3 ];
		}
		elsif( $s[ 0 ] == 2 )
		{
			$recal2{ $s[ 1 ].'_'.$s[ 2 ] } = $s[ 3 ];
		}
	}
	close( QUALMAP );
	
	if( $original1 =~ /\.gz$/ )
	{
		open( ORIG1, "gunzip -c $original1 |" ) or die "Cannot open gzipped fastq\n";
		open( NEW1, "| gzip -c > $dest/$newFastq1" ) or die "Cannot create paired1 fastq file\n";
		
		if( length( $original2 ) > 0 )
		{
			open( ORIG2, "gunzip -c $original2 |" ) or die "Cannot open gzipped fastq\n";
			open( NEW2, "| gzip -c > $dest/$newFastq2" ) or die "Cannot create paired2 fastq file\n";
		}
		
		if( length( $originalUnpaired ) > 0 )
		{
			open( UP, "gunzip -c $originalUnpaired |" ) or die "Cannot open up gzipped fastq\n";
			open( NEWUP, "| gzip -c > $dest/$recalUnpaired" ) or die "Cannot create up fastq file\n";
		}
	}
	else
	{
		open( ORIG1, $original1 ) or die "Cannot open fastq file1\n";
		open( NEW1, ">$dest/$newFastq1" ) or die "Cannot create paired1 fastq file\n";
		
		if( length( $original2 ) > 0 )
		{
			open( ORIG2, $original2 ) or die "Cannot open fastq file2\n";
			open( NEW2, ">$dest/$newFastq2" ) or die "Cannot create paired2 fastq file\n";
		}
		
		if( length( $originalUnpaired ) > 0 )
		{
			open( UP, "$originalUnpaired |" ) or die "Cannot open up gzipped fastq\n";
			open( NEWUP, ">$dest/$recalUnpaired" ) or die "Cannot create up fastq file\n";
		}
	}

	while( <ORIG1> ) 
	{
		print NEW1 $_;
		my $t = <ORIG1>;
		print NEW1 $t;
		$t = <ORIG1>;
		print NEW1 $t;
		$t = <ORIG1>;
		chomp( $t );
		my @s = split( //, $t );
		my $pos = 1;
		foreach( @s )
		{
			if( ! defined $recal1{ $pos.'_'.(unpack("C", $_ ) - 33) } )
			{
				#take the next bases recal value
				while( ! defined $recal1{ $pos.'_'.(unpack("C", $_ ) - 33) } && $pos > -1)
				{
					$pos --;
				}
				
				while( $pos == 0 && ! defined $recal1{ $pos.'_'.(unpack("C", $_ ) - 33) } && $pos < length( $t ) )
				{
					$pos ++;
				}
				
				if( $pos == 0 || ! defined( $recal1{ $pos.'_'.(unpack("C", $_ ) - 33) } ) )
				{
					print NEW1 $_;
				}
				else
				{
					print NEW1 chr( $recal1{ $pos.'_'.(unpack("C", $_ ) - 33) } + 33 );
				}
			}
			else
			{
				print NEW1 chr( $recal1{ $pos.'_'.(unpack("C", $_ ) - 33) } + 33 );
			}
			$pos ++;
		}
		print NEW1 "\n";
	}
	close( ORIG1 );
	close( NEW1 );

	if( length( $original2 ) > 0 )
	{
		while( <ORIG2> ) 
		{
			print NEW2 $_;
			my $t = <ORIG2>;
			print NEW2 $t;
			$t = <ORIG2>;
			print NEW2 $t;
			$t = <ORIG2>;
			chomp( $t );
			my @s = split( //, $t );
			my $pos = 1;
			foreach( @s )
			{
				if( ! defined $recal2{ $pos.'_'.(unpack("C", $_ ) - 33) } )
				{
					#take the next bases recal value
					while( ! defined $recal2{ $pos.'_'.(unpack("C", $_ ) - 33) } && $pos > -1)
					{
						$pos --;
					}
					
					while( $pos == 0 && ! defined $recal2{ $pos.'_'.(unpack("C", $_ ) - 33) } && $pos < length( $t ) )
					{
						$pos ++;
					}
					
					if( $pos == 0 || ! defined( $recal2{ $pos.'_'.(unpack("C", $_ ) - 33) } ) )
					{
						print NEW2 $_;
					}
					else
					{
						print NEW2 chr( $recal2{ $pos.'_'.(unpack("C", $_ ) - 33) } + 33 );
					}
				}
				else
				{
					print NEW2 chr( $recal2{ $pos.'_'.(unpack("C", $_ ) - 33) } + 33 );
				}
				$pos ++;
			}
			print NEW2 "\n";
		}
		close( ORIG2 );
		close( NEW2 );
	}
	
	if( length( $originalUnpaired ) > 0 )
	{
		#HACK - USE THE READ1 RECAL TABLE TO RECALIBRATE.....
		while( <UP> )
		{
			print NEWUP $_;
			my $t = <UP>;
			print NEWUP $t;
			$t = <UP>;
			print NEWUP $t;
			$t = <UP>;
			chomp( $t );
			my @s = split( //, $t );
			my $pos = 1;
			foreach( @s )
			{
				if( ! defined $recal1{ $pos.'_'.(unpack("C", $_ ) - 33) } )
				{
					print "! def $pos _".(unpack("C", $_ ) - 33)."\n";
					exit;
				}
				else
				{
					print NEWUP chr( $recal1{ $pos.'_'.(unpack("C", $_ ) - 33) } + 33 );
				}
				$pos ++;
			}
			print NEWUP "\n";
		}
		close( UP );
		close( NEWUP );
	}
}

=pod
A function to apply the per base position bayesian qual map
=cut
sub applyPosQualMapFastq
{
	croak "Usage: applyPosQualMapFastq qualmap destDir originalFastq1 recalFastq1 qualMapRead[0,1,2]" unless @_ == 5;
	my $qualmap = shift;
	my $dest = shift;
	my $original = shift;
	my $newFastq = shift;
	my $qualMapRead = shift; #which read label (1st column) to use for recal
	
	croak "Cant find fastq1" unless -f $original;
	
	croak "0,1,2 are the valid values for qual map read label to use!" unless $qualMapRead == 0 || $qualMapRead == 1 || $qualMapRead == 2;
	
	croak "Destination not a directory!\n" unless -d $dest;
	my %recal;
	
	open( QUALMAP, $qualmap ) or die "Cannot open qualmap: $!\n";
	while( <QUALMAP> ) 
	{
		chomp;
		my @s = split( /\s+/, $_ );
		if( $s[ 0 ] == $qualMapRead )
		{
			$recal{ $s[ 1 ].'_'.$s[ 2 ] } = $s[ 3 ];
		}
	}
	close( QUALMAP );
	
	if( $original =~ /\.gz$/ )
	{
		open( ORIG1, "gunzip -c $original |" ) or die "Cannot open gzipped fastq\n";
		open( NEW1, "| gzip -c > $dest/$newFastq" ) or die "Cannot create paired fastq file\n";
	}
	else
	{
		open( ORIG1, $original ) or die "Cannot open fastq file1\n";
		open( NEW1, ">$dest/$newFastq" ) or die "Cannot create paired fastq file\n";
	}
	
	while( <ORIG1> ) 
	{
		print NEW1 $_;
		my $t = <ORIG1>;
		print NEW1 $t;
		$t = <ORIG1>;
		print NEW1 $t;
		$t = <ORIG1>;
		chomp( $t );
		my @s = split( //, $t );
		my $pos = 1;
		foreach( @s )
		{
			if( ! defined $recal{ $pos.'_'.(unpack("C", $_ ) - 33) } )
			{
				print "! def $pos _".(unpack("C", $_ ) - 33)."\n";
				exit;
			}
			else
			{
				print NEW1 chr( $recal{ $pos.'_'.(unpack("C", $_ ) - 33) } + 33 );
			}
			$pos ++;
		}
		print NEW1 "\n";
	}
	close( ORIG1 );
	close( NEW1 );
}

=pod
A function to make append to a summary of all the 454 recalibrations per lane
Call this function per lane and it will append to the outputFile parameter.
The outputSummary will be in csv format
=cut
sub summariseMultipleRecalTables
{
	croak "Usage: summariseMultipleRecalTables laneRecalTable outputFile laneName" unless @_ == 3;
	my $table = shift;
	my $outputSummary = shift;
	my $laneName = shift;
	
	if( -f $outputSummary )
	{
		open( OUT, ">>$outputSummary" ) or die "Cannot create summary file\n";
	}
	else
	{
		open( OUT, ">$outputSummary" ) or die "Cannot create summary file\n";
		
		print OUT "Original Q Value,";
		for( my $i = 0; $i < 41; $i ++ )
		{
			print OUT "$i,";
		}
		print OUT "\n";
	}
	
	open( IN, $table ) or die "Cannot open input recal table\n";
	<IN>;
	my %recalVals;
	while( <IN> )
	{
		chomp;
		my @s = split( '\s+', $_ );
		$recalVals{ $s[ 0 ] } = $s[ 3 ];
	}
	close( IN );
	
	print OUT $laneName.",";
	for( my $i = 0; $i < 41; $i ++ )
	{
		print OUT defined( $recalVals{ $i } ) ? $recalVals{ $i }."," : ',';
	}
	print OUT "\n";
	close( OUT );
}

=pod
Summarise multiple Solexa/Solid recalibration tables where data is per pos (not run)
Average the recalibrated values over the run
=cut
sub summariseMultipleRecalTablesSolid
{
	croak "Usage: summariseMultipleRecalTables laneRecalTable outputFile laneName" unless @_ == 3;
	my $table = shift;
	my $outputSummary = shift; #append to this file
	my $laneName = shift;
	
	if( -f $outputSummary )
	{
		open( OUT, ">>$outputSummary" ) or die "Cannot create summary file\n";
	}
	else
	{
		open( OUT, ">$outputSummary" ) or die "Cannot create summary file\n";
		
		print OUT "Run,Read Length,";
		
		for( my $i = 0; $i < 41; $i ++ )
		{
			print OUT "$i,";
		}
		print OUT "\n";
	}
	
	my %newQvals; #oldVal -> [newAvgVal count]
	open( IN, $table ) or die "cannot open table file:$table";
	while( <IN> )
	{
		chomp;
		my @s = split( /\s+/, $_ );
		
		if( ! defined( $newQvals{ $s[ 2 ] } ) )
		{
			$newQvals{ $s[ 2 ] } = [ $s[ 3 ], 1 ];
		}
		else
		{
			$newQvals{ $s[ 2 ] }[ 0 ] = ( ( $newQvals{ $s[ 2 ] }[ 0 ] * $newQvals{ $s[ 2 ] }[ 1 ] ) + $s[ 3 ] ) / ( $newQvals{ $s[ 2 ] }[ 1 ] + 1 );
			$newQvals{ $s[ 2 ] }[ 1 ] ++;
		}
	}
	close( IN );
	
	my $readFile = `ls ../*.recal.fastq.gz | head -1`;
	my $seq = '';
	$readFile = length( $readFile ) > 1 ? substr( $readFile, 0, length( $readFile ) - 1 ) : '';
	if( length( $readFile ) > 1 && -f $readFile )
	{
		$seq = `zcat $readFile | head -2 | tail -1`;
	}
	print "SEQ: $seq\nLength: ".length( $seq )."\n";
	print OUT $laneName.','.length( $seq ).",";

	for( my $i = 0; $i < 41; $i ++ )
	{
		print OUT defined( $newQvals{ $i } ) ? floor( $newQvals{ $i }[ 0 ] )."," : ',';
	}
	print OUT "\n";
	close( OUT );
}

#examines the two fastq files and appends the avg qual per position for each fastq to the output files
#used to compare (via R) the avg qualities per base position before and after recalibration
sub compareAvqQualities
{
	croak "Usage: compareAvqQualities before.fastq after.fastq before_output_file after_output_file" unless @_ == 4;
	my $before = shift;
	my $after = shift;
	my $boutput = shift;
	my $aoutput = shift;
	
	my $readLen = 0;
	if( $before =~ /\.gz$/ )
	{
		open( BEFORE, "gunzip -c $before |" ) or die "Cannot open gzipped fastq\n";
		$readLen = `zcat $before | head -2 | tail -1`;
		$readLen = length( $readLen );
	}
	else
	{
		open( BEFORE, "$before" ) or die "Cannot open before fastq file\n";
		$readLen = `cat $before | head -2 | tail -1`;
		$readLen = length( $readLen );
	}
	
	if( $after =~ /\.gz$/ )
	{
		open( AFTER, "gunzip -c $after |" ) or die "Cannot open gzipped fastq\n";
	}
	else
	{
		open( AFTER, "$after" ) or die "cannot open after fastq file\n";
	}
	
	my @beforeQuals;
	my @afterQuals;
	for( my $i = 0; $i < $readLen; $i ++ )
	{
		push( @beforeQuals, 0 );
		push( @afterQuals, 0 );
	}
	
	my $readCount = 0;
	while( <BEFORE> )
	{
		chomp;
		my $rn1 = $_;
		my $seq1 = <BEFORE>;
		my $qn1 = <BEFORE>;
		my $qual1 = <BEFORE>;
		chomp( $qual1 );
		
		my @t = split( //, $qual1 );
		my $i = 0;
		foreach( @t )
		{
			$beforeQuals[ $i ] = ( ( $beforeQuals[ $i ] * $readCount ) + (unpack("C", $_ ) - 33) ) / ( $readCount + 1 );
			$i ++;
		}
		
		my $rn2 = <AFTER>;
		my $seq2 = <AFTER>;
		my $qn2 = <AFTER>;
		my $qual2 = <AFTER>;
		chomp( $qual2 );
		
		@t = split( //, $qual2 );
		$i = 0;
		foreach( @t )
		{
			$afterQuals[ $i ] = ( ( $afterQuals[ $i ] * $readCount ) + (unpack("C", $_ ) - 33) ) / ( $readCount + 1 );
			$i ++;
		}
		
		$readCount ++;
	}
	
	$before =~ /(.*)\.fastq\.gz/;
	my $prefix = $1;
	open( BOUT, ">>$boutput" ) or die "Cannot write to output file\n";
	print BOUT "$prefix: ";
	foreach( @beforeQuals )
	{
		print BOUT " $_,";
	}
	print BOUT "\n";
	close( BOUT );
	
	open( AOUT, ">>$aoutput" ) or die "Cannot write to output file\n";
	print AOUT "$prefix: ";
	foreach( @afterQuals )
	{
		print AOUT "$_,";
	}
	print AOUT "\n";
	close( AOUT );
}

=pod
Verify that the reads after recalibration are the same as before = apart from quality
=cut
sub verifyBeforeAfterSequence
{
	croak "Usage: verifyBeforeAfterSequence before.fastq after.fastq destination_dir" unless @_ == 3;
	my $before = shift;
	my $after = shift;
	my $dest = shift;
	
	croak "Cant find destination directory: $dest\n" unless -d $dest;
	
	my $totalBases = 0;
	if( $before =~ /\.gz$/ )
	{
		open( BEF, "gunzip -c $before |" ) or die "Cannot open gzipped fastq\n";
	}
	else
	{
		open( BEF, "$before" ) or die "Cannot open before fastq file\n";
	}
	
	if( $after =~ /\.gz$/ )
	{
		open( AFTER, "gunzip -c $after |" ) or die "Cannot open gzipped fastq\n";
	}
	else
	{
		open( AFTER, "$after" ) or die "cannot open after fastq file\n";
	}
	
	my $numSame = 0;
	my $numReads = 0;
	while( <BEF> )
	{
		my $b = $_;
		my $a = <AFTER>;
		
		if( ! defined $a )
		{
			print "rm ".getcwd()."/$after\n";
			croak "More reads in uncalibrated file: $before vs. $after in ".getcwd()."\n";
		}
		
		if( $b ne $a )
		{
			print "rm ".getcwd()."/$after\n";
			croak "E: Read names not same: $b vs. $a\n in ".getcwd()."\n";
		}
		
		$b = <BEF>;
		$a = <AFTER>;
		
		if( $b ne $a )
		{
			print "rm ".getcwd()."/$after\n";
			croak "E: Read sequences not same: $b vs. $a\n in ".getcwd()."\n";
		}
		
		$b = <BEF>;
		$a = <AFTER>;
		
		if( $b ne $a )
		{
			print "rm ".getcwd()."/$after\n";
			croak "E: Qual names not same: $b vs. $a\n in ".getcwd()."\n";
		}
		
		$b = <BEF>;
		$a = <AFTER>;
		
		$numSame ++ unless $b ne $a;
		$numReads ++;
	}
	my $t = <AFTER>;
	if( defined $t )
	{
		print "rm ".getcwd()."/$after\n";
		croak "More reads in recalibrated file: $before vs. $after in ".getcwd()."\n";
	}
	close( AFTER );
	close( BEF );
	
	if( $numReads == 0 )
	{
		print "rm ".getcwd()."/$after\n";
		croak "Empty recalibrated fastq file found: $before vs. $after in ".getcwd()."\n";
	}
	
	#move the recal file to the destination and delete the uncal file
	system( "mv $after $dest" );
	system( "rm -f $before" );
	
	print "$before Verified!";
}

=pod
function for copying the new recalibrated files back to the original place where the reads are
=cut
sub followLinkAndCopy
{
	croak "Usage: symLinkFile originalFile destFile"  unless @_ == 2;
	my $original = shift;
	my $fileToCopy = shift;
	
	croak "Not a sym link: $original\n" unless -l $original;
	
	my $directory = dirname( readlink( $original ) );
	
	copy( $fileToCopy, $directory ) or die "Failed to copy file $fileToCopy to $directory\n";
	
	my $fileName = basename( $fileToCopy );
	if( -s $directory.'/'.$fileName )
	{
		print "SUCCESS: $fileToCopy\n";
		#unlink( $fileToCopy );
	}
}

=pod

=head1 NAME

complement - get the complement of a base.

=head1 SYNOPSIS

complement 

=head1 ARGUMENTS

The function takes the following arguments:

B<> : a string of one character

The function returns the complement of the base - croaks if not a valid base.

=head1 DESCRIPTION

=head1 AUTHOR

Thomas Keane I<tk2@sanger.ac.uk>

=cut

sub complement
{
	if(  @_ != 1 )
	{
		print "Usage: complement base";
		exit;
	}
	
	my $base = shift;
	my $lc = lc( $base );
	if( $lc eq 'n' ){ return 'N'; }
	if( $lc eq 'a' ){ return 'T'; }
	if( $lc eq 't' ){ return 'A'; }
	if( $lc eq 'g' ){ return 'C'; }
	if( $lc eq 'c' ){ return 'G'; }
	croak "Cannot complement base: $base\n";
}

=pod

=head1 NAME

revCompDNA - returns the reverse complement of a DNA string

=head1 SYNOPSIS

revCompDNA string

=head1 ARGUMENTS

The function takes the following arguments:

B<string> : string to reverse complement
	
The function returns a string.

=head1 DESCRIPTION

Reverse complements a dna string argument. The string should only consist of uppercase A, C, G, T characters.

=head1 AUTHOR

Thomas Keane I<tk2@sanger.ac.uk>

=cut

sub revCompDNA
{
	croak "Usage: revCompDNA string\n" unless @_ == 1;
	my $seq = shift;
	$seq= uc( $seq );
	$seq=reverse( $seq );
	
	$seq =~ tr/ACGT/TGCA/;
	return $seq;
}

#a function to take a solid fastq file, convert to base space, and split it into 3 files (if necessary)
sub solid2SensibleFastqs
{
	croak "Usage: solid2SensibleFastqs fastq1 fastq2 desDir\n" unless @_ == 3;
	my $original1 = shift;
	my $original2 = shift;
	my $dest = shift;
	my $prefix = ( split( /_[1-2]\./, basename( $original1 ) ) )[ 0 ];
	
	if( ! -d $dest )
	{
		mkdir( $dest );
		
		croak "Cannot create destination directory\n" unless -d $dest;
	}
	
	my $num1 = 0;
	my $num2 = 0;
	my $numUnpaired = 0;
	if( $original1 =~ /\.gz$/ )
	{
		open( ORIG1, "gunzip -c $original1 |" ) or die "Cannot open gzipped fastq\n";
		open( ORIG2, "gunzip -c $original2 |" ) or die "Cannot open gzipped fastq\n";
		open( PAIRED1, "| gzip -c > $dest/$prefix"."_1.fastq.gz" ) or die "Cannot create paired1 fastq file\n";
		open( PAIRED2, "| gzip -c > $dest/$prefix"."_2.fastq.gz" ) or die "Cannot create paired2 fastq file\n";
		open( UNPAIRED, "| gzip -c > $dest/$prefix.fastq.gz" ) or die "Cannot create unpaired fastq file\n";
		
	}
	else
	{
		open( ORIG1, $original1 ) or die "Cannot open fastq file1\n";
		open( ORIG2, $original2 ) or die "Cannot open fastq file2\n";
		open( PAIRED1, ">$dest/$prefix"."_1.fastq" ) or die "Cannot create paired1 fastq file\n";
		open( PAIRED2, ">$dest/$prefix"."_2.fastq" ) or die "Cannot create paired2 fastq file\n";
		open( UNPAIRED, ">$dest/$prefix.fastq" ) or die "Cannot create unpaired fastq file\n";
	}
	
	while( <ORIG1> )
	{
		chomp;
		my $header1 = $_;
		my $header2 = <ORIG2>;
		chomp( $header2 );
		
		my $seq1 = <ORIG1>;
		chomp( $seq1 );
		#solid - so clip off the first 2 bases
		$seq1 = substr( $seq1, 2 );
		$seq1 =~ tr/0123./ACGTN/;
		
		my $seq2 = <ORIG2>;
		chomp( $seq2 );
		#solid - so clip off the first 2 bases
		$seq2 = substr( $seq2, 2 );
		$seq2 =~ tr/0123./ACGTN/;
		
		my $qheader1 = <ORIG1>;
		chomp( $qheader1 );
		my $qheader2 = <ORIG2>;
		chomp( $qheader2 );
		
		my $qual1 = <ORIG1>;
		chomp( $qual1 );
		#solid - so clip off the first 2 bases
		$qual1 = substr( $qual1, 2 );

		my $qual2 = <ORIG2>;
		chomp( $qual2 );
		#solid - so clip off the first 2 bases
		$qual2 = substr( $qual2, 2 );
		
		if( ( $seq1 =~ /^N+$/ ) )
		#if( ( $seq1 =~ /^T\.+$/ || $seq1 =~ /^G\.+$/ ) && ( $seq2 !~ /^T\.+$/ && $seq2 !~ /^G\.+$/ ) )
		{
			print UNPAIRED "$header2\n$seq2\n$qheader2\n$qual2\n";
			$numUnpaired ++;
		}
		elsif( ( $seq2 =~ /^N+$/ ) )
		#elsif( ( $seq1 !~ /^T\.+$/ && $seq1 !~ /^G\.+$/ ) && ( $seq2 =~ /^T\.+$/ || $seq2 =~ /^G\.+$/ ) )
		{
			print UNPAIRED "$header1\n$seq1\n$qheader1\n$qual1\n";
			$numUnpaired ++;
		}
		elsif( $seq2 !~ /^N+$/ && $seq1 !~ /^N+$/ )
		#elsif( ( $seq1 !~ /^T\.+$/ && $seq1 !~ /^G\.+$/ ) && ( $seq2 !~ /^T\.+$/ && $seq2 !~ /^G\.+$/ ) )
		{
			print PAIRED1 "$header1\n$seq1\n$qheader1\n$qual1\n";
			print PAIRED2 "$header2\n$seq2\n$qheader2\n$qual2\n";
			$num1 ++;
			$num2 ++;
		}
	}
	close( ORIG1 );
	close( PAIRED1 );
	close( PAIRED2 );
	close( ORIG2 );
	
	if( $num1 == 0 )
	{
		if( $original1 =~ /\.gz$/ )
		{
			unlink( "$dest/$prefix"."_1.fastq.gz" );
			unlink( "$dest/$prefix"."_2.fastq.gz" );
		}
		else
		{
			unlink( "$dest/$prefix"."_1.fastq" );
			unlink( "$dest/$prefix"."_2.fastq" );
		}
	}
	elsif( $numUnpaired == 0 )
	{
		if( $original1 =~ /\.gz$/ )
		{
			unlink( "$dest/$prefix.fastq.gz" );
		}
		else
		{
			unlink( "$dest/$prefix.fastq" );
		}
	}
}

sub verifyDirectorySolexaFastqSizes
{
	croak "Usage: verifyDirectorySolexaFastqSizes directory" unless @_ == 1;
	my $directory = shift;
	
	foreach my $fastq1 (`ls $directory/*_1.fastq.gz $directory/*_1.recal.fastq.gz`)
	{
		chomp( $fastq1 );
		
		$fastq1 =~ /(.*)_1(\..*fastq\.gz)$/;
		my $fastq2 = $1."_2$2";
		
		if( -f $fastq2 )
		{
			my $cmd = qq[ bsub -q normal -o $directory/verify.o -e $directory/verify.e perl -w -e "use Recalibration;Recalibration::verifySolexaFastqSizes( '$fastq1', '$fastq2' );" ];
			print $cmd."\n";
			system($cmd);
			#verifySolexaFastqSizes( $fastq1, $fastq2 );
		}
		else
		{
			print "No fastq2 for: $fastq1 $fastq2\n";
		}
	}
}

sub verifySolexaFastqSizes
{
	croak "Usage: verifySolexaFastqSizes 1.fastq 2.fastq" unless @_ == 2;
	my $fastq1 = shift;
	my $fastq2 = shift;
	
	my $numReads1 = 0;
	my $numReads2 = 0;
	if( $fastq1 =~ /\.gz$/ )
	{
		open( F1, "gunzip -c $fastq1 |" ) or die "Cannot open gzipped fastq1\n";
	}
	else
	{
		open( F1, "$fastq1" ) or die "Cannot open fastq1 file\n";
	}
	
	if( $fastq2 =~ /\.gz$/ )
	{
		open( F2, "gunzip -c $fastq2 |" ) or die "Cannot open gzipped fastq2\n";
	}
	else
	{
		open( F2, "$fastq2" ) or die "cannot open fastq2 file\n";
	}
	
	while( <F1> )
	{
		if( $_ =~ /^@.*/ )
		{
			$numReads1 ++;
			<F1>;<F1>;<F1>;
		}
		
		if( eof( F2 ) )
		{
			croak "ERROR: Fastq1 more reads: $fastq1 $fastq2\n";
			close( F1 );
			close( F2 );
			return;
		}
		else
		{
			my $t = <F2>;
			if( $_ =~ /^@.*/ )
			{
				$numReads2 ++;
				<F2>;<F2>;<F2>;
			}
		}
	}
	close( F1 );
	close( F2 );
	
	if( $numReads1 != $numReads2 )
	{
		print "ERROR: Number of reads different: $numReads1 vs. $numReads2 in $fastq1 $fastq2\n";
		return;
	}
	print "GOOD: $fastq1 $fastq2\n";
}

1;
