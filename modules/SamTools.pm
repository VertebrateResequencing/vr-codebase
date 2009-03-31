package SamTools;
use strict;
use Carp;
use File::Spec;
use Cwd;
use Digest::MD5  qw(md5 md5_hex md5_base64);

use AssemblyTools;

=pod
	A method to convert SSAHA cigar output to SAM format for unpaired reads
	Assumes that the cigar filea are ordered the same as the fastq files
=cut
my $MIN_HIT_LENGTH_454 = 30;
sub ssaha2samUnpaired
{
	croak "Usage: ssaha2sam fastq1 ssaha_cigar1 read_group output_sam_file\n" unless @_ == 4;
	my $fastq1 = shift;
	my $cigar1 = shift;
	my $read_group = shift;
	my $sam = shift;
	
	croak "Cant find fastq1 file\n" unless -f $fastq1;
	croak "Cant find SSAHA cigar1 file\n" unless -f $cigar1;
	
	if( $fastq1 =~ /\.gz$/ )
	{
		open( READS1, "gunzip -c $fastq1 |" ) or die "Cannot open gzipped fastq1 file\n";
	}
	else
	{
		open( READS1, $fastq1 ) or die "Failed to open reads file";
	}
	
	open( SAM, "|gzip -c > $sam" ) or die "Cannot create SAM file: $!\n";
	
	#this could be a problem due to memory!
	my $ref = AssemblyTools::hash_ssaha_cigar_output( $cigar1 );
	print "Hashed cigar1\n";
	my %cigar1Hash = %{ $ref };
	
	while( <READS1> )
	{
		chomp;
		my $name1 = substr( $_, 1 );
		$name1 = (split( /\s+/, $name1 ) )[ 0 ];
		my $seq1 = <READS1>;
		chomp( $seq1 );
		<READS1>;
		my $quals1 = <READS1>;
		chomp( $quals1 );
		
		my @hits1 = ();
		if( defined( $cigar1Hash{ $name1 } ) )
		{
			@hits1 = @{ $cigar1Hash{ $name1 } };
		}
		
		if( @hits1 > 0 )
		{
			my $samCigar1 = ssahaCigar2SamCigar( $hits1[ 0 ], $seq1 );
			
			my $mapScore = detemineMappingScore( $cigar1Hash{ $name1 } );
			
			if( $mapScore > -1 )
			{
				my @s1 = split( /\s+/, $hits1[ 0 ] );
				
				my $flag = determineSamFlagUnpaired( $hits1[ 0 ] );
				
				$seq1 = AssemblyTools::revCompDNA( $seq1 ) if( $s1[ 4 ] eq '-' );
				$quals1 = reverse( $quals1 ) if( $s1[ 4 ] eq '-' );
				
				print SAM "$name1\t".$flag."\t".$s1[ 5 ]."\t".$s1[ 6 ]."\t".$mapScore."\t".$samCigar1."\t*\t0\t0\t".$seq1."\t".$quals1."\tRG:Z:".$read_group."\n";
			}
		}
	}
	close( READS1 );
	close( SAM );
}

=pod
	A method to convert SSAHA cigar output to SAM format
	Assumes that the cigar filea are ordered the same as the fastq files
=cut
sub ssaha2samPaired
{
	croak "Usage: ssaha2sam fastq1 ssaha_cigar1 fastq2 ssaha_cigar2 insert_size read_group output_sam_file\n" unless @_ == 7;
	my $fastq1 = shift;
	my $cigar1 = shift;
	my $fastq2 = shift;
	my $cigar2 = shift;
	my $insert = shift;
	my $read_group = shift;
	my $sam = shift;
	
	croak "Cant find fastq1 file\n" unless -f $fastq1;
	croak "Cant find fastq2 file\n" unless -f $fastq2;
	croak "Cant find SSAHA cigar1 file\n" unless -f $cigar1;
	croak "Cant find SSAHA cigar2 file\n" unless -f $cigar2;
	
	if( $fastq1 =~ /\.gz$/ )
	{
		open( READS1, "gunzip -c $fastq1 |" ) or die "Cannot open gzipped fastq1 file\n";
		open( READS2, "gunzip -c $fastq2 |" ) or die "Cannot open gzipped fastq2 file\n";
	}
	else
	{
		open( READS1, $fastq1 ) or die "Failed to open reads file";
		open( READS2, $fastq2 ) or die "Failed to open reads file";
	}
	
	open( SAM, "|gzip -c > $sam" ) or die "Cannot create SAM file: $!\n";
	
	#this could be a problem due to memory!
	my $ref = AssemblyTools::hash_ssaha_cigar_output( $cigar1 );
	print "Hashed cigar1\n";
	my %cigar1Hash = %{ $ref };
	$ref = AssemblyTools::hash_ssaha_cigar_output( $cigar2 );
	print "Hashed cigar2\n";
	my %cigar2Hash = %{ $ref };
	
	while( <READS1> )
	{
		chomp;
		my $name1 = substr( $_, 1 );
		$name1 = (split( /\s+/, $name1 ) )[ 0 ];
		my $seq1 = <READS1>;
		chomp( $seq1 );
		<READS1>;
		my $quals1 = <READS1>;
		chomp( $quals1 );
		
		my $name2 = <READS2>;
		chomp( $name2 );
		$name2 = (split( /\s+/, $name2 ) )[ 0 ];
		$name2 = substr( $name2, 1 );
		my $seq2 = <READS2>;
		chomp( $seq2 );
		<READS2>;
		my $quals2 = <READS2>;
		chomp( $quals2 );
		
		my @hits1 = ();
		if( defined( $cigar1Hash{ $name1 } ) )
		{
			@hits1 = @{ $cigar1Hash{ $name1 } };
		}
		
		my @hits2 = ();
		if( defined( $cigar2Hash{ $name2 } ) )
		{
			@hits2 = @{ $cigar2Hash{ $name2 } };
		}
		
		if( @hits2 > 0 && @hits1 > 0 )
		{
			my $mapQual = 0;
			my $flag = 0;
			my $hit1;
			my $hit2;
			
			#from the list of read1 hits and read2 hits - get 1 hit per read
			
			my @hits = @{ determineHits( \@hits1, \@hits2, $insert ) };
			
			my $samCigar1 = ssahaCigar2SamCigar( $hits[ 0 ], $seq1 );
			my $samCigar2 = ssahaCigar2SamCigar( $hits[ 1 ], $seq2 );
			my @s1 = split( /\s+/, $hits[ 0 ] );
			my @s2 = split( /\s+/, $hits[ 1 ] );
			
			my @flags = @{ determineSamFlagPaired( $hits[ 0 ], $hits[ 1 ], $insert ) };
			
			my $mapScore1 = detemineMappingScore( \@hits1 );
			$mapScore1 = 0 if( $mapScore1 < 0 );
			my $mapScore2 = detemineMappingScore( \@hits2 );
			$mapScore2 = 0 if( $mapScore2 < 0 );
			
			if( $flags[ 2 ] == 1 && $mapScore1 > -1 && $mapScore2 > -1 ) #if paired consistently - then double the mapping scores (RD suggest)
			{
				$mapScore1 = $mapScore1 * 2;
				$mapScore2 = $mapScore2 * 2;
				
				$mapScore1 = 255 if( $mapScore1 > 255 );
				$mapScore2 = 255 if( $mapScore2 > 255 );
			}
			
			$seq1 = AssemblyTools::revCompDNA( $seq1 ) if( $s1[ 4 ] eq '-' );
			$quals1 = reverse( $quals1 ) if( $s1[ 4 ] eq '-' );
			$seq2 = AssemblyTools::revCompDNA( $seq2 ) if( $s2[ 4 ] eq '-' );
			$quals2 = reverse( $quals2 ) if( $s2[ 4 ] eq '-' );
			
			if( $mapScore1 > -1 && $mapScore2 > -1 )
			{
				print SAM "$name1\t".$flags[ 0 ]."\t".$s1[ 5 ]."\t".$s1[ 6 ]."\t".$mapScore1."\t".$samCigar1."\t".($s2[ 4 ] eq $s1[ 4 ] ? "=" : $s2[ 4 ])."\t".$s2[ 5 ]."\t".( determineInsertSize( $s1[ 6], $s1[ 7 ], $s2[ 6 ], $s2[ 7 ] ) )."\t".$seq1."\t".$quals1."\tRG:Z:".$read_group."\n";
				print SAM "$name2\t".$flags[ 1 ]."\t".$s2[ 5 ]."\t".$s2[ 6 ]."\t".$mapScore2."\t".$samCigar2."\t".($s1[ 4 ] eq $s2[ 4 ] ? "=" : $s1[ 4 ])."\t".$s1[ 5 ]."\t".( determineInsertSize( $s1[ 6 ], $s1[ 7 ], $s2[ 6 ], $s2[ 7 ] ) )."\t".$seq2."\t".$quals2."\tRG:Z:".$read_group."\n";
			}
			elsif( $mapScore1 > -1 )
			{
				print SAM "$name1\t".$flags[ 0 ]."\t".$s1[ 5 ]."\t".$s1[ 6 ]."\t".$mapScore1."\t".$samCigar1."\t*\t0\t0\t".$seq1."\t".$quals1."\tRG:Z:".$read_group."\n";
			}
			elsif( $mapScore2 > -1 )
			{
				print SAM "$name2\t".$flags[ 1 ]."\t".$s2[ 5 ]."\t".$s2[ 6 ]."\t".$mapScore2."\t".$samCigar2."\t*\t0\t0\t".$seq2."\t".$quals2."\tRG:Z:".$read_group."\n";
			}
		}
		elsif( @hits1 > 0 )
		{
			my $samCigar1 = ssahaCigar2SamCigar( $hits1[ 0 ], $seq1 );
			my @s1 = split( /\s+/, $hits1[ 0 ] );
			
			my $flag = determineSamFlagUnpaired( $hits1[ 0 ] );
			$flag += hex( '0x001' );
			$flag += hex( '0x008' );
			
			my $mapScore = detemineMappingScore( \@hits1 );
			
			if( $mapScore > -1 )
			{
				$seq1 = AssemblyTools::revCompDNA( $seq1 ) if( $s1[ 4 ] eq '-' );
				$quals1 = reverse( $quals1 ) if( $s1[ 4 ] eq '-' );
				
				print SAM "$name1\t".$flag."\t".$s1[ 5 ]."\t".$s1[ 6 ]."\t".$mapScore."\t".$samCigar1."\t*\t0\t0\t".$seq1."\t".$quals1."\tRG:Z:".$read_group."\n";
			}
			
			if( $name2 ne $name1 )
			{
				while( <READS2> )
				{
					chomp;
					if( $_ =~ /^\@$name1/ )
					{
						<READS2>;
						<READS2>;
						<READS2>;
						last;
					}
				}
			}
		}
		elsif( @hits2 > 0 )
		{
			my $samCigar2 = ssahaCigar2SamCigar( $hits2[ 0 ], $seq2 );
			my @s2 = split( /\s+/, $hits2[ 0 ] );
			
			my $flag = determineSamFlagUnpaired( $hits2[ 0 ] );
			$flag += hex( '0x001' );
			$flag += hex( '0x008' );
			
			my $mapScore = detemineMappingScore( \@hits2 );
			
			if( $mapScore > -1 )
			{
				$seq2 = AssemblyTools::revCompDNA( $seq2 ) if( $s2[ 4 ] eq '-' );
				$quals2 = reverse( $quals2 ) if( $s2[ 4 ] eq '-' );
				
				print SAM "$name2\t".$flag."\t".$s2[ 5 ]."\t".$s2[ 6 ]."\t".$mapScore."\t".$samCigar2."\t*\t0\t0\t".$seq2."\t".$quals2."\tRG:Z:".$read_group."\n";
			}
			
			if( $name1 ne $name2 )
			{
				while( <READS1> )
				{
					chomp;
					if( $_ =~ /^\@$name2/ )
					{
						<READS1>;
						<READS1>;
						<READS1>;
						last;
					}
				}
			}
		}
		else
		{
			#print "No hits for: $name1 and $name2\n";
		}
	}
	close( READS1 );
	close( READS2 );
	close( SAM );
}

sub ssahaCigar2SamCigar
{
	croak "Usage: ssahaCigar2SamCigar ssaha_cigar read_sequence\n" unless @_ == 2;
	my $ssahaCigar = shift;
	my $seq = shift;
	
	my @s = split( /\s+/, $ssahaCigar );
	
	#work out the soft clip of the reads (if any)
	my $softClip = '';
	if( $s[ 4 ] eq '+' && $s[ 2 ] > 1 )
	{
		$softClip = ($s[ 2 ] - 1).'S';
	}
	elsif( $s[ 4 ] eq '-' && $s[ 2 ] < length( $seq ) )
	{
		$softClip = (length( $seq ) - $s[ 2 ]).'S';
	}
	
	my $endSoftClip = '';
	if( $s[ 4 ] eq '+' && $s[ 3 ] < length( $seq ) )
	{
		$endSoftClip = (length( $seq ) - $s[ 3 ]).'S';
	}
	elsif( $s[ 4 ] eq '-' && $s[ 3 ] > 1 )
	{
		$endSoftClip = ($s[ 3 ] - 1).'S';
	}
	
	#get the cigar string - need to switch the M/I/D for SAM
	my @cigar = splice( @s, 10 );
	for( my $i = 0; $i < @cigar; $i += 2 )
	{
		my $t = $cigar[ $i ];
		$cigar[ $i ] = $cigar[ $i + 1 ];
		$cigar[ $i + 1 ] = $t;
	}
	#print "SAM: ".$softClip.join( '', @cigar );exit;
	return $softClip.join( '', @cigar ).$endSoftClip;
}

=pod
	Takes 2 arrays of SSAHA cigar hits and tries to determine 1 hit per read
	Returns a ref to array: [ hit1, hit2 ]
=cut
sub determineHits
{
	croak "Usage: determineHits hits1 hits2 expected_insert_size\n" unless @_ == 3;
	my $hits1ref = shift;
	my $hits2ref = shift;
	my $expectedInsert = shift;
	
	my @hits1 = @{ $hits1ref };
	my @hits2 = @{ $hits2ref };

	if( @hits1 == 1 && @hits2 == 1 )
	{
		my @hit1Split = split( /\s+/, $hits1[ 0 ] );
		my @hit2Split = split( /\s+/, $hits2[ 0 ] );
		my $insert = abs( determineInsertSize( $hit1Split[ 2 ], $hit1Split[ 3 ], $hit2Split[ 2 ], $hit2Split[ 3 ] ) );
		my @out = ( $hits1[ 0 ], $hits2[ 0 ] );
		return \@out;
	}
	elsif( @hits1 == 1 && @hits2 > 1 )
	{
		my @hit1Split = split( /\s+/, $hits1[ 0 ] );
		my $hit1Chr = $hit1Split[ 5 ];
		
		my $closestOffset = -1;
		my $closestHit = $hits2[ 0 ];
		foreach( @hits2 )
		{
			chomp;
			my @hit2Split = split( /\s+/, $_ );
			if( $_ =~ /\s+$hit1Chr\s+/ && ( abs( determineInsertSize( $hit1Split[ 2 ], $hit1Split[ 3 ], $hit2Split[ 2 ], $hit2Split[ 3 ] ) ) - $expectedInsert < $closestOffset || $closestOffset == -1 ) )
			{
				$closestOffset = abs( determineInsertSize( $hit1Split[ 2 ], $hit1Split[ 3 ], $hit2Split[ 2 ], $hit2Split[ 3 ] ) ) - $expectedInsert;
				my $closestHit = $_;
			}
		}
		my @out = ( $hits1[ 0 ], $closestHit );
		return \@out;
	}
	elsif( @hits1 > 1 && @hits2 == 1 )
	{
		my @hit2Split = split( /\s+/, $hits2[ 0 ] );
		my $hit2Chr = $hit2Split[ 5 ];
		
		my $closestOffset = -1;
		my $closestHit = $hits1[ 0 ];
		foreach( @hits1 )
		{
			chomp;
			my @hit1Split = split( /\s+/, $_ );
			if( $_ =~ /\s+$hit2Chr\s+/ && ( abs( determineInsertSize( $hit1Split[ 2 ], $hit1Split[ 3 ], $hit2Split[ 2 ], $hit2Split[ 3 ] ) ) - $expectedInsert < $closestOffset || $closestOffset == -1 ) )
			{
				$closestOffset = abs( determineInsertSize( $hit1Split[ 2 ], $hit1Split[ 3 ], $hit2Split[ 2 ], $hit2Split[ 3 ] ) ) - $expectedInsert;
				my $closestHit = $_;
			}
		}
		my @out = ( $closestHit, $hits2[ 0 ] );
		return \@out;
	}
	else
	{
		#pick the best combination of hits on same chr closest to insert size
		my $closestOffset = -1;
		my $hit1 = '';
		my $hit2 = '';
		foreach( @hits1 )
		{
			my @hit1Split = split( /\s+/, $_ );
			my $h1 = $_;
			foreach( @hits2 )
			{
				my @hit2Split = split( /\s+/, $_ );
				if( $hit1Split[ 5 ] eq $hit2Split[ 5 ] && ( abs( determineInsertSize( $hit1Split[ 2 ], $hit1Split[ 3 ], $hit2Split[ 2 ], $hit2Split[ 3 ] ) ) - $expectedInsert < $closestOffset || $closestOffset == -1 ) )
				{
					$closestOffset = abs( determineInsertSize( $hit1Split[ 2 ], $hit1Split[ 3 ], $hit2Split[ 2 ], $hit2Split[ 3 ] ) ) - $expectedInsert;
					$hit1 = $h1;
					$hit2 = $_;
				}
			}
		}
		
		if( length( $hit1 ) == 0 )
		{
			$hit1 = $hits1[ 0 ];
			$hit2 = $hits2[ 0 ];
		}
		
		my @out = ( $hit1, $hit2 );
		return \@out;
	}
}

=pod
	A function to determine a mapping score for ssaha cigar lines
=cut
sub detemineMappingScore
{
	croak "Usage: detemineMappingScore \@hits\n" unless @_ == 1;
	my $ref = shift;
	my @hits = @{ $ref };
	
	my @hit1 = split( /\s+/, $hits[ 0 ] );
	if( abs( $hit1[ 2 ] - $hit1[ 3 ] ) < $MIN_HIT_LENGTH_454 )
	{
		return -1;
	}
	
	if( @hits == 1 )
	{
		return 255;
	}
	else
	{
		my @hit2 = split( /\s+/, $hits[ 1 ] );
		
		my $d = ( $hit1[ 9 ] - $hit2[ 9 ] ) - ( 2 * scalar( @hits ) );
		
		return $d < 256 ? $d : 255;
	}
}

sub determineInsertSize
{
	croak "Usage: determineInsertSize st1 end1 st2 end2\n" unless @_ == 4;
	my $st1 = shift;
	my $end1 = shift;
	my $st2 = shift;
	my $end2 = shift;
	
	my $insert = abs( $st1 - $st2 );
	my $posNegInsert = 0;
	if( abs( $st1 - $end2 ) > $insert )
	{
		$insert = abs( $st1 - $end2 );
		$posNegInsert = $st1 - $end2;
	}
	
	if( abs( $end1 - $st2 ) > $insert )
	{
		$insert = abs( $end1 - $st2 );
		$posNegInsert = $end1 - $st2;
	}
	
	if( abs( $end1 - $end2 ) > $insert )
	{
		$insert = abs( $end1 - $end2 );
		$posNegInsert = $end1 - $end2;
	}
	
	return $posNegInsert;
}

sub determineSamFlagPaired
{
	croak "Usage: determineSamFlags cigar1 cigar2 expected_insert\n" unless @_ == 3;
	my $cigar1 = shift;
	my $cigar2 = shift;
	my $expected = shift;
	
	my $flag1 = hex( '0x0001' );
	my $flag2 = hex( '0x0001' );
	my @s1 = split( /\s+/, $cigar1 );
	my @s2 = split( /\s+/, $cigar2 );
	
	my $consistent = 0;
	
	if( length( $cigar1 ) > 0 && length( $cigar2 ) > 0 )
	{
		if( $s1[ 5 ] eq $s2[ 5 ] && abs( abs( determineInsertSize( $s1[ 2 ], $s1[ 3 ], $s2[ 2 ], $s2[ 3 ] ) ) - $expected ) < $expected * 3 && $s1[ 4 ] eq $s2[ 4 ] ) #if within 3 times the insert size and on same strand (454) - then paired normally
		{
			$flag1 += hex( '0x0002' );
			$flag2 += hex( '0x0002' );
		}
		$flag1 += hex( '0x0010' ) if( $s1[ 4 ] eq '-' );
		$flag2 += hex( '0x0010' ) if( $s2[ 4 ] eq '-' );
		$flag1 += hex( '0x0040' );
		$flag2 += hex( '0x0080' );
		
		$consistent = 1;
	}
	elsif( length( $cigar1 ) > 0 )
	{
		$flag1 += hex( '0x008' );
		$flag1 += hex( '0x0010' ) if( $s1[ 4 ] eq '-' );
	}
	elsif( length( $cigar2 ) > 0 )
	{
		$flag2 += hex( '0x008' );
		$flag2 += hex( '0x0010' ) if( $s2[ 4 ] eq '-' );
	}
	
	my $ref = [ $flag1, $flag2, $consistent ];
	return $ref;
}

sub determineSamFlagUnpaired
{
	croak "Usage: determineSamFlags cigar\n" unless @_ == 1;
	my $cigar = shift;
	my @s = split( /\s+/, $cigar );
	
	my $flag = 0;
	$flag += hex( '0x0010' ) if( $s[ 4 ] eq '-' );
	
	return $flag;
}

1;
