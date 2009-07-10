package SequencingQuality;
use strict;
use Carp;
use File::Spec;
use Cwd;
use File::Basename;

my $MAQ_CMD = "maq";
my $LSF_NORMAL_PREFIX = 'bsub -q normal -o out.o -e out.e';
my $LSF_LONG_PREFIX = 'bsub -q long -o out.o -e out.e';
my $R_LOCATION = '/software/R-2.7.1/bin/R';
my $SLX_LINKER = 'A*GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGA';
my $IMAGE_MAGIC = '/usr/bin/convert';

#a perl module of functions for assessing the quality of sequencing reads

sub gcContentHistogram
{
	croak "Usage: gcContentHistogram gc_hist_output1 fastq1 gc_bmp_output [gc_hist_output2 fastq2]" unless @_ == 3 || @_ == 5;
	
	my $gc_hist1 = shift;
	my $fastq1 = shift;
	my $gc_bmp = shift;
	my $fastq2 = '';
	my $gc_hist2 = '';
	if( @_ == 2 )
	{
		$gc_hist2 = shift;
		$fastq2 = shift;
	}
	
	my %gcContent1;
	my %gcContent2;
	
	my $avgContent1 = 0;
	my $totalReads1 = 0;
	my $avgContent2 = 0;
	my $totalReads2 = 0;
	
	if( $fastq1 =~ /.*\.gz$/ )
	{
		open( IN, "gunzip -c $fastq1 |" ) or die "Cannot open gzipped fastq1\n";
	}
	else
	{
		open( IN, $fastq1 ) or die "Cannot open fastq1 file\n";
	}
	
	while( <IN> )
	{
		chomp;
		my $name = substr( $_, 1 );
		my $seq = <IN>;
		chomp( $seq );
		
		<IN>;
		<IN>;
		
		my @seq1 = split( //, $seq );
		my $gc = 0;
		foreach( @seq1 )
		{
			if( $_ eq 'G' || $_ eq 'g' || $_ eq 'C' || $_ eq 'c' )
			{
				$gc ++;
			}
		}
		
		$gc = $gc / @seq1;
		
		if( defined $gcContent1{ $gc } )
		{
			$gcContent1{ $gc } ++;
		}
		else
		{
			$gcContent1{ $gc } = 1;
		}
		
		$avgContent1 = ( ( $avgContent1 * $totalReads1 ) + $gc ) / ( $totalReads1 + 1 );
		$totalReads1 ++;
	}
	close( IN );
	
	if( length( $fastq2 ) > 0 )
	{
		if( $fastq2 =~ /.*\.gz$/ )
		{
			open( IN, "gunzip -c $fastq2 |" ) or die "Cannot open gzipped fastq1\n";
		}
		else
		{
			open( IN, $fastq2 ) or die "Cannot open fastq2 file\n";
		}
		
		while( <IN> )
		{
			chomp;
			my $name = substr( $_, 1 );
			my $seq = <IN>;
			chomp( $seq );
			
			<IN>;
			<IN>;
			
			my @seq1 = split( //, $seq );
			my $gc = 0;
			foreach( @seq1 )
			{
				if( $_ eq 'G' || $_ eq 'g' || $_ eq 'C' || $_ eq 'c' )
				{
					$gc ++;
				}
			}
		
			$gc = $gc / @seq1;
		
			if( defined $gcContent2{ $gc } )
			{
				$gcContent2{ $gc } ++;
			}
			else
			{
				$gcContent2{ $gc } = 1;
			}
			
			$avgContent2 = ( ( $avgContent2 * $totalReads2 ) + $gc ) / ( $totalReads2 + 1 );
			$totalReads2 ++;
		}
		close( IN );
	}
	
	#also print out a GC histogram
	open( OUT, ">$gc_hist1" ) or die "Cannot create file\n";
	foreach( sort( keys( %gcContent1 ) ) )
	{
		print OUT $_.' '.$gcContent1{ $_ }."\n";
	}
	close( OUT );
	
	if( length( $fastq2 ) > 0 )
	{
		open( OUT, ">$gc_hist2" ) or die "Cannot create file\n";
		foreach( sort( keys( %gcContent2 ) ) )
		{
			print OUT $_.' '.$gcContent2{ $_ }."\n";
		}
		close( OUT );
	}
	
	#write an R script and submit the LSF job to run it
	open( R, ">gc.R" ) or die "cannot open R script\n";
	print R 'gc_1<-read.table(\''.$gc_hist1.'\');';
        print R 'bitmap(\''.$gc_bmp.'\');';
        print R 'plot(gc_1[,1],gc_1[,2],xlab="%GC",ylab="Frequency",type="l");';
	
	if( length( $fastq2  ) > 0 )
	{
		print R 'gc_2<-read.table(\''.$gc_hist2.'\');';
		print R 'lines(gc_2[,1],gc_2[,2],type="l");';

		print "$fastq2 Average GC: $avgContent2\n";
	}
        print R 'dev.off();';
	close( R );
	
        system( $R_LOCATION.' CMD BATCH gc.R;'.$IMAGE_MAGIC.' gc.bmp gc.gif;rm gc.bmp' );
	#unlink( "$fastq1._gc.R" );
	
	print "$fastq1 Average GC: $avgContent1\n";
}

sub duplicateReads
{
	croak "Usage: duplicateReads fastq read_prefix_size min_hits output_file" unless @_ == 4;
	my $fastq = shift;
	my $prefix_size = shift;
	my $minHits = shift;
	my $output = shift;
	
	my %reads;
	
	my $totalBp = 0;
	my $seqLength = 0;
	my $totalReads = 0;
	
	if( $fastq =~ /.*\.fastq\.gz$/ )
	{
		open( FASTQ, "gunzip -c $fastq |" ) or die "Cannot open gzipped fastq1\n";
	}
	else
	{
		open( FASTQ, "$fastq" ) or die "Cannot open fastq1\n";
	}
	
	while( <FASTQ> )
	{
		chomp;
		my $name = substr( $_, 1 );
		my $seq = <FASTQ>;
		chomp( $seq );
		
		$totalBp += length( $seq );
		$totalReads ++;
		if( $prefix_size > 0 )
		{
			my $p1 = substr( $seq, 0, $prefix_size );
			my $p2 = substr( $seq, (length( $seq ) / 2), $prefix_size );
			$seq = $p1.$p2;
		}
		
		<FASTQ>;
		<FASTQ>;
		
		if( defined $reads{ $seq } )
		{
			$reads{ $seq } ++;
		}
		else
		{
			$reads{ $seq } = 1;
		}
	}
	close( FASTQ );
	
	#sort by the max number of hits
	my $count = 0;
	my $readLength = $totalBp / $totalReads;
	open( OUT, ">$output" ) or die "Cannot create output file: $output\n";
	foreach( sort{ $reads{ $b } <=> $reads{ $a } } keys %reads )
	{
		if( $reads{ $_ } > $minHits )
		{
			my $percent = ( ( $reads{ $_ } * $readLength ) / $totalBp ) * 100;
			$percent = sprintf("%.2f", $percent);
			print OUT $_.' '.$reads{ $_ }.' '.$percent."%\n";
		}
	}
	close( OUT );
}

sub splitSolexaReadsMaq
{
	croak "Usage: splitSolexaReadsMaq in.fastq" unless @_ == 1;
	
	open( READS, "$_[ 0 ]" ) or die "cannot open reads file\n";
	open( L, ">$_[ 0 ].left" ) or die "Cannot create LEFT file\n";
	open( R, ">$_[ 0 ].right" ) or die "Cannot create LEFT file\n";
	while( <READS> )
	{
		chomp;
		my $name = $_;
		print L $name."\n";
		print R $name."\n";
		my $t = <READS>;
		chomp( $t );
		print L substr( $t, 0, length( $t ) / 2 )."\n";
		print R substr( $t, length( $t ) / 2 )."\n";
		$t = <READS>;
		print L $t;
		print R $t;
		$t = <READS>;
		chomp( $t );
		print L substr( $t, 0, length( $t ) / 2 )."\n";
		print R substr( $t, length( $t ) / 2 )."\n";
	}
	close( READS );
	close( L );
	close( R );
}

sub mergeMaqMaps
{
	croak "Usage: mergeMaqMaps directory map-prefix output_map" unless @_ == 3;
	my $dir = shift;
	my $prefix = shift;
	my $output = shift;
	
	if( -f $output ){unlink($output);}
	
	#pre: first need to count the total number of reads (mapstat)
	my $totalPreReads = 0;
	
	my $files = '';
	my $cmd = $MAQ_CMD.' mapmerge '.$output.' ';
	opendir( DIR, $dir ) or die "Cannot open directory of map files\n";
	foreach( readdir( DIR ) )
	{
		if( $_ =~ /^$prefix.*\.[0-9]+\.map$/ )
		{
			$totalPreReads += `$MAQ_CMD mapstat $_ | grep "Total number of reads" | awk '{print $6}'`;
			
			$files .= $_.' ';
		}
	}
	close( DIR );
	
	#$cmd .= $files.';rm '.$files;
	$cmd .= $files;
	print $cmd."\n";
	system( $cmd );
	
	#post: count the total number of reads after operation
	my $totalPostReads = `$MAQ_CMD mapstat $_ | grep "Total number of reads" | awk '{print $6}'`;
	
	if( $totalPreReads != $totalPostReads )
	{
		print "ERROR: Numbers of reads after merge not equal to total reads before merge: $totalPreReads vs. $totalPostReads\n";
	}
}
#=pod
sub readCoverageDistribution
{
	croak "Usage: readCoverageDistribution maq-map ref.bfa max_coverage" unless @_ == 3;
	my $map = shift;
	my $bfa = shift;
	my $max = shift;
	
	if( ! -f $map ){print "Cant find map file\n";exit;}
	
	#convert the maq map to a pileup - and then read off the coverage
	my %coverageHistogram;
	my %chrSizes;
	my $currentSize = 0;
	my $currentChr = '';
	open( PILEUP, $MAQ_CMD." pileup $bfa $map|" ) or die "Cannot read pileup\n";
	while( <PILEUP> )
	{
		chomp;
		
		my @s = split( /\s+/, $_ );
		
		if( length( $currentChr ) == 0 )
		{
			$currentChr = $s[ 0 ];
			$currentSize ++;
		}
		elsif( $currentChr ne $s[ 0 ] )
		{
			$chrSizes{ $currentChr } = $currentSize;
			
			$currentChr = $s[ 0 ];
			$currentSize = 1;
		}
		else
		{
			$currentSize ++;
		}
		
		if( defined( $coverageHistogram{ $s[ 0 ] } ) )
		{
			$coverageHistogram{ $s[ 0 ] } += $s[ 3 ];
		}
		else
		{
			$coverageHistogram{ $s[ 0 ] } = $s[ 3 ];
		}
	}
	close( PILEUP );
	
	#normalise the coverage per kb
	my %normalisedCoverage;
	foreach( keys( %chrSizes ) )
	{
		$normalisedCoverage{ $_ } = ( $coverageHistogram{ $_ } / $chrSizes{ $_ } ) * 1000;
	}
	
	open( DIST, ">$map.cov_dist" ) or die "Cannot create cov distribution file\n";
	foreach( sort( keys( %coverageHistogram ) ) )
	{
		if( $_ =~ /^\d+$/ )
		{
			print DIST $_." ".$coverageHistogram{ $_ }."\n";
		}
	}
	close( DIST );
	
	#print an R script and graph the distribution
	open( R, ">$map.cov_dist.R" ) or die "Cannot create R script\n";
	print R "cov_dist <- read.table('$map.cov_dist');";
	print R "cov_dist$max <- cov_dist[,2][cov_dist[,2]<$max];";
	print R "bitmap('$map.cov_dist.bmp');";
	print R "plot(density(cov_dist$max);";
	print R "dev.off()";
	close( R );
	
	system( $R_LOCATION.' CMD BATCH '.$map.'.cov_dist.R' );
	
	#print the normalised coverage per chr
	open( DIST, ">$map.cov_per_kb" ) or die "Cannot create cov_per_kb distribution file\n";
	foreach( sort( keys( %normalisedCoverage ) ) )
	{
		if( $_ =~ /^\d+$/ )
		{
			print DIST $_." ".$normalisedCoverage{ $_ }."\n";
		}
	}
	close( DIST );
	
	#print an R script and graph the distribution
	open( R, ">$map.cov_per_kb.R" ) or die "Cannot create R script\n";
	print R "cov_per_kb <- read.table('$map.cov_per_kb');";
	print R "bitmap('$map.cov_dist.bmp');";
	print R "plot(cov_per_kb[,1], cov_per_kb[,2],xlab='Chromosome',ylab='Coverage per kb',main='Coverage per kb');";
	print R "dev.off()";
	close( R );
	
	system( $R_LOCATION.' CMD BATCH '.$map.'.cov_per_kb.R' );
}
#=cut

#needs to be fully fixed for paired fastq's!
sub generateStatsReport
{
	croak "Usage: generateStatsReport mapViewFile rawMapViewFile outcsv dupPrefixSize fastq1 [fastq2]" unless @_ == 5 || @_ == 6;
	my $mapView = shift;
	my $rawMapView = shift;
	my $csv = shift;
	my $prefixSize = shift;
	my $fastq1 = shift;
	
	my $fastq2 = '';
	if( @_ == 1 )
	{
		$fastq2 = shift;
	}
	
	croak "Usage: generateStatsReport mapViewFile rawMapViewFile outcsv dupPrefixSize fastq1 [fastq2]" unless $prefixSize =~ /^\d+$/;
	
	#get the num reads, pairs, lengths, and num bases
	if( $fastq1 =~ /.*\.fastq\.gz$/ )
	{
		open( FASTQ1, "gunzip -c $fastq1 |" ) or die "Cannot open gzipped fastq1\n";
		
		if( length( $fastq2 ) > 0 )
		{
			open( FASTQ2, "gunzip -c $fastq2 |" ) or die "Cannot open gzipped fastq2\n";
		}
	}
	else
	{
		open( FASTQ1, "$fastq1" ) or die "Cannot open fastq1\n";
		
		if( length( $fastq2 ) > 0 )
		{
			open( FASTQ2, "$fastq2" ) or die "Cannot open fastq2\n";
		}
	}
	
	my $numBases = 0;
	my $numPairs = 0;
	my $readLength = 0;
	my $linkerReads = 0;
	my %readBaseComposition;
	my $polyA;
	my $avgContent1 = 0;
	my $totalReads = 0;
	my $avgContent2 = 0;
	while( <FASTQ1> )
	{
		chomp;
		my $seq1 = <FASTQ1>;
		chomp( $seq1 );
		$seq1 = uc( $seq1 );
		$numBases += length( $seq1 );
		<FASTQ1>;
		<FASTQ1>;
		
		<FASTQ2>;
		my $seq2 = <FASTQ2>;
		chomp( $seq2 );
		$numBases += length( $seq2 );
		<FASTQ2>;
		<FASTQ2>;
		
		my $gc = 0;
		my @seq1_ = split( //, $seq1 );
		foreach( @seq1_ )
		{
			if( $_ eq 'G' || $_ eq 'g' || $_ eq 'C' || $_ eq 'c' )
			{
				$gc ++;
			}
		}
		$gc = $gc / @seq1_;
		$avgContent1 = ( ( $avgContent1 * $numPairs ) + $gc ) / ( $numPairs + 1 );
		
		$gc = 0;
		my @seq2_ = split( //, $seq2 );
		foreach( @seq2_ )
		{
			if( $_ eq 'G' || $_ eq 'g' || $_ eq 'C' || $_ eq 'c' )
			{
				$gc ++;
			}
		}
		$gc = $gc / @seq2_;
		$avgContent2 = ( ( $avgContent2 * $numPairs ) + $gc ) / ( $numPairs + 1 );
		
		#find if the linker is in the seqs
		if( $prefixSize > 0 )
		{
			if( substr( $seq1, 0, $prefixSize ) =~ /^$SLX_LINKER/ )
			{
				$linkerReads ++;
			}
			
			if( substr( $seq2, 0, $prefixSize ) =~ /^$SLX_LINKER/ )
			{
				$linkerReads ++;
			}
		}
		else
		{
			if( $seq1 =~ /^$SLX_LINKER/ )
			{
				$linkerReads ++;
			}
			
			if( $seq2 =~ /^$SLX_LINKER/ )
			{
				$linkerReads ++;
			}
		}
		
		my @s = split( //, $seq1 );
		my $numA = 0;
		foreach( @s )
		{
			$numA ++ unless lc( $_ ) ne 'A';
		}
		$polyA ++ unless ( $numA / length( $seq1 ) ) < 0.95;
		
		@s = split( //, $seq2 );
		$numA = 0;
		foreach( @s )
		{
			$numA ++ unless lc( $_ ) ne 'A';
		}
		$polyA ++ unless ( $numA / length( $seq2 ) ) < 0.95;
		
		for( my $i = 0; $i < @s; $i ++ )
		{
			if( ! defined $readBaseComposition{ $i } )
			{
				$readBaseComposition{ $i } = [ 0, 0, 0, 0 ];
			}
			
			$readBaseComposition{ $i }[ 0 ] ++ unless $_ ne 'A';
			$readBaseComposition{ $i }[ 1 ] ++ unless $_ ne 'C';
			$readBaseComposition{ $i }[ 2 ] ++ unless $_ ne 'G';
			$readBaseComposition{ $i }[ 3 ] ++ unless $_ ne 'T';
		}
		
		$readLength = length( $seq1 ) unless $readLength > 0;
		$numPairs ++;
	}
	close( FASTQ1 );
	close( FASTQ2 );
	
	my $a = 'a = c( ';
	my $c = 'c = c( ';
	my $g = 'g = c( ';
	my $t = 't = c( ';
	foreach( sort( keys( %readBaseComposition ) ) )
	{
		$a .= $readBaseComposition{ $_ }[ 0 ].', ';
		$c .= $readBaseComposition{ $_ }[ 1 ].', ';
		$g .= $readBaseComposition{ $_ }[ 2 ].', ';
		$t .= $readBaseComposition{ $_ }[ 3 ].', ';
	}
	
	$a = substr( $a, 0, length( $a ) - 2 ).' );';
	$c = substr( $c, 0, length( $c ) - 2 ).' );';
	$g = substr( $g, 0, length( $g ) - 2 ).' );';
	$t = substr( $t, 0, length( $t ) - 2 ).' );';
	
	#create a plot of the A,C,G,T bases in each position
	open( R, ">baseFreqs.R" ) or die "Cannot create R script for barplot\n";
	print R "$a\n";
	print R "$t\n";
	print R "$g\n";
	print R "$c\n";
	print R "bases = 1:$readLength;";
	print R "bmp( 'baseFreqs.bmp' );";
	print R "plot( a, bases, type='l', xlab='Base Position', ylab='Frequency', col=c(\"red\") );";
	print R "lines( c, bases, type='l', col=c(\"green\"));";
	print R "lines( g, bases, type='l', col=c(\"blue\"));";
	print R "lines( t, bases, type='l', col=c(\"yellow\"));";
	#print R "legend(350,300,lty=1,col=c("green", "red"),legend=c("LHS", "RHS"))";
	print R "dev.off();";
	close( R );
	
	system( $R_LOCATION.' CMD BATCH baseFreqs.R' );
	
	#get the num reads mapped at Q0 and Q30
	my $numMapped = 0;
	my $numMappedQ30 = 0;
	my %mappedQHist;
	
	if( $mapView =~ /\.gz$/ )
	{
		open( VIEW, "gunzip -c $mapView |" ) or die "Cannot open gzipped mapView\n";
	}
	else
	{
		open( VIEW, "$mapView" ) or die "Cannot open mapView\n";
	}
	
	while( <VIEW> )
	{
		chomp;
		my @s = split( /\s+/, $_ );
		
		if( $s[ 5 ] != 192 )
		{
			$numMapped ++;
		}
		
		if( $s[ 6 ] >= 30 )
		{
			$numMappedQ30 ++;
		}
		
		if( defined( $mappedQHist{ $s[ 6 ] } ) )
		{
			$mappedQHist{ $s[ 6 ] } ++;
		}
		else
		{
			$mappedQHist{ $s[ 6 ] } ++
		}
	}
	close( VIEW );
	
	my $rawMapped = 0;
	if( $rawMapView =~ /\.gz$/ )
	{
		$rawMapped = `zcat $rawMapView | awk '$6!=192' | wc -l`;
	}
	else
	{
		$rawMapped = `cat $rawMapView | awk '$6!=192' | wc -l`;
	}
	my $numDuplicates = $rawMapped - $numMapped;
	
	my $exists = -f $csv ? 1 : 0;
	open( STATS, ">>$csv" ) or die "Cannot create stats file\n";
	print STATS "Location,Pairs,Length,Bases,Linker,%,PolyA,%,Duplicates,%,Mapped,%,Q30,%\n" unless $exists == 1;
	my $cwd = getcwd();
	my $numReads = $numPairs * 2;
	print STATS "$cwd,$numPairs,$readLength,$numBases,$linkerReads,".($linkerReads/$numReads)."$polyA,".($polyA/$numReads).",$numDuplicates,".($numDuplicates/$numReads).",$numMapped,".($numMapped/$numReads).",$numMappedQ30,".($numMappedQ30/$numReads)."\n";
	close( STATS );
}

sub simpleQCStats
{
	croak "Usage: generateStatsReport mapViewFile rawMapViewFile outcsv fastqcheck1 [fastqcheck2]" unless @_ == 4 || @_ == 5;
	my $mapView = shift;
	my $rawMapView = shift;
	my $csv = shift;
	my $fastqCheck1 = shift;
	
	my $numBases1 = 0;
	my $numBases2 = 0;
	
	my $fastqCheck2 = '';
	my @fqc2;
	if( @_ == 1 )
	{
		$fastqCheck2 = shift;
		croak "cant find fastqcheck file2" unless -f $fastqCheck2;
		
		open( FQC, $fastqCheck2 ) or die "Cant open fastqcheck2";
		my $t = <FQC>;
		chomp( $t );
		@fqc2 =  split( /\s+/, $t );
		close( FQC );
		
		$numBases2 = $fqc2[ 2 ];
	}
	
	croak "cant find fastqcheck file1" unless -f $fastqCheck1;
	
	open( FQC, $fastqCheck1 ) or die "Cant open fastqcheck1";
	my $t = <FQC>;
	chomp( $t );
	my @fqc1 =  split( /\s+/, $t );
	close( FQC );
	
	$numBases1 = $fqc1[ 2 ];
	my $numPairs = $fqc1[ 0 ];
	my $readLength = $fqc1[ 7 ];
	
	#get the num reads mapped at Q0 and Q30
	my $numMapped = 0;
	my $numMappedQ30 = 0;
	my %mappedQHist;
	my $numMappedInPairs = 0;
	
	if( $mapView =~ /\.gz$/ )
	{
		open( VIEW, "gunzip -c $mapView |" ) or die "Cannot open gzipped mapView\n";
	}
	else
	{
		open( VIEW, "$mapView" ) or die "Cannot open mapView\n";
	}
	
	while( <VIEW> )
	{
		chomp;
		my @s = split( /\s+/, $_ );
		
		if( $s[ 5 ] != 192 )
		{
			$numMapped ++;
		}
		
		if( $s[ 6 ] >= 30 )
		{
			$numMappedQ30 ++;
		}
		
		if( defined( $mappedQHist{ $s[ 6 ] } ) )
		{
			$mappedQHist{ $s[ 6 ] } ++;
		}
		else
		{
			$mappedQHist{ $s[ 6 ] } ++
		}
	}
	close( VIEW );
	
	my $rawMapped = 0;
	if( $rawMapView =~ /\.gz$/ )
	{
		$rawMapped = `zcat $rawMapView | awk '{if(\$6!=192){print}}' | wc -l`;
	}
	else
	{
		$rawMapped = `cat $rawMapView | awk '\$6!=192' | wc -l`;
	}
	my $numDuplicates = $rawMapped - $numMapped;
	
	my $exists = -f $csv ? 1 : 0;
	open( STATS, ">$csv" ) or die "Cannot create stats file\n";
	my $cwd = getcwd();
	my $numReads = $numPairs * 2;
	my $numBases = $numBases1 + $numBases2;
	print STATS "$cwd,$numPairs,$readLength,$numBases,$numDuplicates,".($numDuplicates/$numReads).",$numMapped,".($numMapped/$numReads).",$numMappedQ30,".($numMappedQ30/$numReads)."\n";
	close( STATS );
}

sub createInsertGraphLowMem
{
	croak "Usage: createInsertGraph mapview maxInsert output_file" unless @_ == 3;
	my $mapview = shift;
	my $maxInsert = shift;
	my $output_file = shift;
	
	croak "Cant find mapview file: $mapview" unless -f $mapview;
	croak "Usage: createInsertGraph mapview maxInsert" unless $maxInsert =~ /^\d+$/;
	
	if( $mapview =~ /\.gz$/ )
	{
		open( MV, "gunzip -c $mapview |" ) or die "Cannot open file: $mapview\n";
	}
	else
	{
		open( MV, "$mapview" ) or die "Cannot open file: $mapview\n";
	}
	
	my %inserts;
	$inserts{ 0 } = 0;
	while( <MV> )
	{
		chomp;
		my @s = split( /\s+/, $_ );
		
		if( $s[ 5 ] == 130 || $s[ 5 ] == 18 || $s[ 5 ] == 20 || ( $maxInsert > 2000 && $s[ 5 ] == 4 ) )
		{
			my $insert = abs( $s[ 4 ] );
			if( $insert < $maxInsert )
			{
				if( defined $inserts{ $insert } )
				{
					$inserts{ $insert } ++;
				}
				else
				{
					$inserts{ $insert } = 1;
				}
			}
		}
		elsif( $s[ 5 ] != 192 ) #i.e. reads mapped but not paired
		{
			#$inserts{ 0 } ++;
		}
	}
	close( MV );
	
	open( I, ">$mapview.inserts" )  or die "cannot create inserts file\n";
	foreach( sort {$a <=> $b}( keys( %inserts ) ) )
	{
		print I $_." ".$inserts{ $_ }."\n";
	}
	close( I );
	
	open( R, ">$mapview.insert.R" ) or die "Cannot create R script\n";
	print R "inserts<-read.table('$mapview.inserts');\n";
	print R "bitmap('$output_file');\n";
	print R "plot(inserts[,1],inserts[,2],main='Insert Sizes',xlab='Insert Size', ylab='Frequency', type='p', axes='F' );\n";
	print R "at.labels=axis(2);\n";
	print R "at.labels=axis(1)\n";
	print R "axis(1,at=seq(at.labels[1], at.labels[length(at.labels)], (at.labels[2]-at.labels[1])/4), labels=F )\n";
	print R "dev.off();\n";
	close( R );
	
	my $output_file_gif = $output_file;
	$output_file_gif =~ s/bmp/gif/;
	system( $R_LOCATION." CMD BATCH $mapview.insert.R;$IMAGE_MAGIC $output_file $output_file_gif" );
}

sub createInsertGraph
{
	croak "Usage: createInsertGraph mapview maxInsert output_file" unless @_ == 3;
	my $mapview = shift;
	my $maxInsert = shift;
	my $output_file = shift;
	
	croak "Cant find mapview file: $mapview" unless -f $mapview;
	croak "Usage: createInsertGraph mapview maxInsert" unless $maxInsert =~ /^\d+$/;
	
	if( $mapview =~ /\.gz$/ )
	{
		system( 'zcat '.$mapview.' | awk \'{if($6==130||$6==18||$6==20){if($5<0){print $1,-$5}else{print $1,$5}}}\' | sort | uniq > '.$mapview.'.inserts' );
	}
	else
	{
		system( 'cat '.$mapview.' | awk \'{if($6==130||$6==18||$6==20){if($5<0){print $1,-$5}else{print $1,$5}}}\' | sort | uniq > '.$mapview.'.inserts' );
	}
	
	open( R, ">$mapview.insert.R" ) or die "Cannot create R script\n";
	print R "inserts<-read.table('$mapview.inserts');\n";
	print R "inserts$maxInsert <- inserts[,2][inserts[,2]<$maxInsert];\n";
	print R "bitmap('$output_file');\n";
	print R "plot(density(inserts$maxInsert),main='Inserts $mapview');\n";
	print R "dev.off();\n";
	close( R );
	
	system( $R_LOCATION." CMD BATCH $mapview.insert.R" );
}

sub fastq2fasta
{
	croak "Usage: fastq2fasta file_name output_file_name" unless @_ == 2;
	my $file = shift;
	my $output_file = shift;
	
	if( ! (-f $file) ){croak "Cannot find file: $file\n";}
	if( ! (-s $file) ){croak "Empty file: $file\n";}
	
	open( READS_FILE, $file ) or die "Error: Cannot open fastq file\n";
	open( FASTA_FILE, ">$output_file" ) or die "Error: Cannot create fasta file\n";
	
	my $line;
	while( <READS_FILE> )
	{
		chomp;
		my $name = substr( $_, 1 );
		my $seq = <READS_FILE>;
		chomp( $seq );
		print FASTA_FILE ">$name\n$seq\n";
		<READS_FILE>;
		<READS_FILE>;
	}
	close( READS_FILE );
	close( FASTA_FILE );
}

1;
