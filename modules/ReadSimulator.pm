package ReadSimulator;
use strict;
use Carp;
use Digest::MD5 'md5_hex';
#use Math::Random::OO::Normal;

use AssemblyTools;

sub generatePerfectReadsUnpaired
{
	croak "Usage: generatePerfectReadsUnpaired reference_fasta read_length num_reads output_fastq" unless @_ == 4; 
	my $reference_fasta = shift;
	my $read_length = shift;
	my $num_reads = shift;
	my $output = shift;
	
	my $read_quals = '';
	for( my $i = 0; $i < $read_length; $i ++ )
	{
		$read_quals .= '>';
	}
	
	my $ref = AssemblyTools::createHash( $reference_fasta );
	my %reference = %{ $ref };
	my @reference_keys = keys( %reference );

	foreach( @reference_keys )
	{
		$reference{ $_ } = (split( /\n/, $reference{ $_ } ) )[ 1 ];
	}
	@reference_keys = keys( %reference );
	
	open( OUT, ">$output" ) or die "Cannot create output file\n";
	my $chr = 0;
	for( my $i = 0; $i < $num_reads; $i ++ )
	{
		$chr = int( rand( @reference_keys ) );
		
		my $start = int( rand( length( $reference{ $reference_keys[ $chr ] } ) ) );
		my $read = uc( substr( $reference{ $reference_keys[ $chr ] }, $start, $read_length ) );
		
		if( length( $read ) == $read_length )
		{
			#reverse complement some of the reads...
			if( rand( 1 ) > 0.5 )
			{
				#rev comp the read
				my $new_read = AssemblyTools_unstable::revCompDNA( $read );
				#print OUT '@'.$reference_keys[ $chr ].'_'.($start + 35).'_'.$start.'_R'."\n".$new_read."\n+".$reference_keys[ $chr ].'_'.($start + 35).'_'.$start."_R\n".$read_quals."\n";
				print OUT '@'.$reference_keys[ $chr ].'_'.($start + $read_length).'_'.$start.'_R'."\n".$new_read."\n+\n".$read_quals."\n";
			}
			else
			{
				#print OUT '@'.$reference_keys[ $chr ].'_'.$start.'_'.($start + 35)."_F\n".$read."\n+".$reference_keys[ $chr ].'_'.$start.'_'.($start + 35)."_F\n".$read_quals."\n";
				print OUT '@'.$reference_keys[ $chr ].'_'.$start.'_'.($start + $read_length).'_F'."\n".$read."\n+\n".$read_quals."\n";
			}
		}
		else
		{
			$i --;
		}
	}
	close( OUT );
}

sub generateSNPReads
{
	croak "Usage: generateSNPReads reference_fasta read_length num_reads snps_per_read output_fastq" unless @_ == 5; 
	my $reference_fasta = shift;
	my $read_length = shift;
	my $num_reads = shift;
	my $snps_per_read = shift;
	my $output = shift;
	
	my $read_quals = '';
	for( my $i = 0; $i < $read_length; $i ++ )
	{
		$read_quals .= '>';
	}
	
	my $ref = AssemblyTools::createHash( $reference_fasta );
	my %reference = %{ $ref };
	my @reference_keys = keys( %reference );

	foreach( @reference_keys )
	{
		$reference{ $_ } = (split( /\n/, $reference{ $_ } ) )[ 1 ];
	}
	@reference_keys = keys( %reference );
	
	open( OUT, ">$output" ) or die "Cannot create output file\n";
	my $chr = 0;
	for( my $i = 0; $i < $num_reads; $i ++ )
	{
		$chr = int( rand( @reference_keys ) );
		
		my $start = int( rand( length( $reference{ $reference_keys[ $chr ] } ) ) );
		my $read = uc( substr( $reference{ $reference_keys[ $chr ] }, $start, $read_length ) );
		
		if( length( $read ) == $read_length )
		{
			#add the snps to the read
			my $positionsString;
			for( my $j = 0; $j < $snps_per_read; $j ++ )
			{
				my $pos = int( rand( $read_length ) );
				
				$read = substr( $read, 0, $pos ).randomNucleotideExclude( substr( $read, $pos, 1 ) ).substr( $read, $pos + 1 );
				
				$positionsString .= '_'.$pos;
			}
			
			#reverse complement some of the reads...
			if( rand( 1 ) > 0.5 )
			{
				#rev comp the read
				my $new_read = AssemblyTools_unstable::revCompDNA( $read );
				print OUT '@'.$reference_keys[ $chr ].'_'.($start + $read_length).'_'.$start.$positionsString.'_R'."\n".$new_read."\n+\n".$read_quals."\n";
			}
			else
			{
				print OUT '@'.$reference_keys[ $chr ].'_'.$start.'_'.($start + $read_length).$positionsString.'_F'."\n".$read."\n+\n".$read_quals."\n";
			}
		}
		else
		{
			$i --;
		}
	}
	close( OUT );
}

sub randomNucleotideExclude
{
	croak "Usage: randomNucleotide excludedNucleotide" unless @_ == 1; 
	
	while()
	{
		my $nuc = int( rand( 4 ) );
		
		if( $nuc == 0 && $_[ 0 ] ne 'A' ){return 'A'}
		elsif( $nuc == 1 && $_[ 0 ] ne 'C' ){return 'C'}
		elsif( $nuc == 2 && $_[ 0 ] ne 'G' ){return 'G'}
		elsif( $_[ 0 ] ne 'T' ){return 'T';}
	}
	return 'A';
}

sub generateDeletionReads
{
	croak "Usage: generateDeletionReads reference_fasta read_length num_reads del_size output_fastq" unless @_ == 5; 
	my $reference_fasta = shift;
	my $read_length = shift;
	my $num_reads = shift;
	my $del_size = shift;
	my $output = shift;
	
	my $read_quals = '';
	for( my $i = 0; $i < $read_length - $del_size; $i ++ )
	{
		$read_quals .= '>';
	}
	
	my $ref = AssemblyTools::createHash( $reference_fasta );
	my %reference = %{ $ref };
	my @reference_keys = keys( %reference );

	foreach( @reference_keys )
	{
		$reference{ $_ } = (split( /\n/, $reference{ $_ } ) )[ 1 ];
	}
	@reference_keys = keys( %reference );
	
	open( OUT, ">$output" ) or die "Cannot create output file\n";
	my $chr = 0;
	for( my $i = 0; $i < $num_reads; $i ++ )
	{
		$chr = int( rand( @reference_keys ) );
		
		my $start = int( rand( length( $reference{ $reference_keys[ $chr ] } ) ) );
		my $read = uc( substr( $reference{ $reference_keys[ $chr ] }, $start, $read_length ) );
		
		if( length( $read ) == $read_length )
		{
			#generate the random insert sequence
			my $dpos = int( rand( length( $read ) ) );
			while( $dpos - $del_size < 0 )
			{
				$dpos = int( rand( length( $read ) ) );
			}
			$read = substr( $read, 0, $dpos - $del_size ).substr( $read, $dpos );
			
			
			#reverse complement some of the reads...
			if( rand( 1 ) > 0.5 )
			{
				#rev comp the read
				my $new_read = AssemblyTools_unstable::revCompDNA( $read );
				print OUT '@'.$reference_keys[ $chr ].'_'.($start + $read_length).'_'.$start.'_'.$dpos.'_'.$del_size.'_R'."\n".$new_read."\n+\n".$read_quals."\n";
			}
			else
			{
				print OUT '@'.$reference_keys[ $chr ].'_'.$start.'_'.($start + $read_length).'_'.$dpos.'_'.$del_size.'_F'."\n".$read."\n+\n".$read_quals."\n";
			}
		}
		else
		{
			$i --;
		}
	}
	close( OUT );
}

sub generateInsertReads
{
	croak "Usage: generateInsertReads reference_fasta read_length num_reads insert_size output_fastq" unless @_ == 5; 
	my $reference_fasta = shift;
	my $read_length = shift;
	my $num_reads = shift;
	my $insert_size = shift;
	my $output = shift;
	
	my $read_quals = '';
	for( my $i = 0; $i < $read_length + $insert_size; $i ++ )
	{
		$read_quals .= '>';
	}
	
	my $ref = AssemblyTools::createHash( $reference_fasta );
	my %reference = %{ $ref };
	my @reference_keys = keys( %reference );

	foreach( @reference_keys )
	{
		$reference{ $_ } = (split( /\n/, $reference{ $_ } ) )[ 1 ];
	}
	@reference_keys = keys( %reference );
	
	open( OUT, ">$output" ) or die "Cannot create output file\n";
	my $chr = 0;
	for( my $i = 0; $i < $num_reads; $i ++ )
	{
		$chr = int( rand( @reference_keys ) );
		
		my $start = int( rand( length( $reference{ $reference_keys[ $chr ] } ) ) );
		my $read = uc( substr( $reference{ $reference_keys[ $chr ] }, $start, $read_length ) );
		
		if( length( $read ) == $read_length )
		{
			#generate the random insert sequence
			my $ipos = int( rand( length( $read ) ) );
			$read = substr( $read, 0, $ipos ).randomNucleotides( $insert_size ).substr( $read, $ipos );
			
			#reverse complement some of the reads...
			if( rand( 1 ) > 0.5 )
			{
				#rev comp the read
				my $new_read = AssemblyTools_unstable::revCompDNA( $read );
				print OUT '@'.$reference_keys[ $chr ].'_'.($start + $read_length).'_'.$start.'_'.$ipos.'_'.$insert_size.'_R'."\n".$new_read."\n+\n".$read_quals."\n";
			}
			else
			{
				print OUT '@'.$reference_keys[ $chr ].'_'.$start.'_'.($start + $read_length).'_'.$ipos.'_'.$insert_size.'_F'."\n".$read."\n+\n".$read_quals."\n";
			}
		}
		else
		{
			$i --;
		}
	}
	close( OUT );
}

sub randomNucleotides
{
	croak "Usage: randomNucleotide size" unless @_ == 1;
	
	my $insert = '';
	for( my $i = 0; $i < $_[ 0 ]; $i ++ )
	{
		$insert .= randomNucleotide();
	}
	return $insert;
}

sub randomNucleotide
{
	while()
	{
		my $nuc = int( rand( 4 ) );
		
		if( $nuc == 0 ){return 'A'}
		elsif( $nuc == 1 ){return 'C'}
		elsif( $nuc == 2 ){return 'G'}
		else{return 'T';}
	}
	return 'A';
}

=pod
	A function to determine % uniqueness of a genome at a given read length
	Generates all possible reads at given read length and counts how many are unique
=cut

sub generateUniqueReadsStats
{
	croak "Usage: generateUniqueReadsStats fasta read_length skip_size" unless @_ == 3;
	my $contigs = $_[ 0 ];
	my $readLength = $_[ 1 ];
	my $skip = $_[ 2 ];

	croak "Usage: generateUniqueReadsStats fasta read_length skip" unless -s $contigs;
	croak "Usage: generateUniqueReadsStats fasta read_length skip" unless $readLength > 0;
	croak "Usage: generateUniqueReadsStats fasta read_length skip" unless $skip > 0;
	
	my %readCounts;
	my $totalReads = 0;
	my $currentChr = '';
	open( IN, $contigs ) or die "Cannot open contigs\n";
	while( <IN> )
	{
		chomp;
		if( $_ =~ /^>/ && length( $currentChr ) > 0 )
		{
			print "Starting $_\n";
			for( my $i = 0; $i < length( $currentChr ); $i += $skip )
			{
				if( length( $currentChr ) - $i < $readLength )
				{
					next;
				}
				$totalReads += 2;
				
				my $read = uc( substr( $currentChr, $i, $readLength ) );
				
				my $revComp = uc( reverse( $read ) );
				$revComp =~ tr/ACGT/TGCA/;
				
				if( ! defined $readCounts{ $read } )
				{
					$readCounts{ $read } = 1;
				}
				else
				{
					$readCounts{ $read } ++;
				}
				
				if( ! defined $readCounts{ $revComp } )
				{
					$readCounts{ $revComp } = 1;
				}
				else
				{
					$readCounts{ $revComp } ++;
				}
			}
			
			$currentChr = '';
		}
		else
		{
			$currentChr .= $_;
		}
	}
	close( IN );
	
			for( my $i = 0; $i < length( $currentChr ); $i += $skip )
			{
				if( length( $currentChr ) - $i < $readLength )
				{
					next;
				}
				$totalReads ++;
				
				my $read = substr( $currentChr, $i, $readLength );
				if( ! defined $readCounts{ md5_hex( $read ) } )
				{
					$readCounts{ md5_hex( $read ) } = 1;
				}
				else
				{
					$readCounts{ md5_hex( $read ) } ++;
				}
			}
	
	print "Total Reads: $totalReads\n";
	
	#count unique reads
	my $numUnique = 0;
	my @keys = keys( %readCounts );
	foreach( @keys )
	{
		if( $readCounts{ $_ } == 1 )
		{
			$numUnique ++;
		}
	}
	
	my $proportion = $numUnique / $totalReads;
	print "UNIQ: $readLength\t$proportion\n";
}

sub generatePerfectReadsPaired
{
	croak "Usage: generatePerfectReadsPaired reference_fasta read_length total_bases insert_size output_fastq" unless @_ == 5; 
	my $reference_fasta = shift;
	my $readLength = shift;
	my $total_bases = shift;
	my $insert_size = shift;
	my $output = shift;
	
	#uniform base quality - Q30
	my $quals = '';
	for( my $i = 0; $i < $readLength; $i ++ )
	{
		$quals .= '>';
	}
	
	my $ref = AssemblyTools::fastaHash( $reference_fasta );
	my %reference = %{ $ref };
	my @reference_keys = keys( %reference );

	foreach( @reference_keys )
	{
		$reference{ $_ } = (split( /\n/, $reference{ $_ } ) )[ 1 ];
	}
	@reference_keys = keys( %reference );
	
	open( OUT, ">$output" ) or die "Cannot create output file: $!\n";
	my $totalBp = 0;

	#pick a random entry and generate a read pair from it
	my @keys = keys( %reference );
	my $insertGenerator = Math::Random::OO::Normal->new( $insert_size, $insert_size * 0.2 );
	$insertGenerator -> seed( rand() );
	while( $totalBp < $total_bases )
	{
		my $key = int( rand( scalar( @keys ) ) );
		my $chr = $reference{ $keys[ $key ] };

		my $startPos = int( rand( length( $chr ) ) );
		my $insertSize = int( $insertGenerator -> next );
		#my $insertSize = int( rand( $insert_size * 0.4 ) ) + ( $insert_size - ($insert_size * 0.2) );
		
		if( $startPos + $insertSize > length( $chr ) )
		{
			next;
		}
		else
		{
			#generate a read pair
			#switch the orientation of the reads randomly....
			if( rand( 1 ) > 0.5 )
			{
				my $read1 = substr( $chr, $startPos, $readLength );
				my $read2 = AssemblyTools::revCompDNA( substr( $chr, $startPos + $insertSize - $readLength, $readLength ) );
				next if ( $read1 =~ /^N+$/ || $read2 =~ /^N+$/ );
				
				print OUT "@".$key."_".$startPos.".p1k\n".$read1."\n+\n$quals\n"."@".$key."_".($startPos+$insertSize).".q1k\n".$read2."\n+\n$quals\n";
			}
			else
			{
				my $read1 = AssemblyTools::revCompDNA( substr( $chr, $startPos, $readLength ) );
				my $read2 =  substr( $chr, $startPos + $insertSize - $readLength, $readLength );
				next if ( $read1 =~ /^N+$/ || $read2 =~ /^N+$/ );
				
				print OUT "@".$key."_".$startPos.".p1k\n".$read1."\n+\n$quals\n"."@".$key."_".($startPos+$insertSize).".q1k\n".$read2."\n+\n$quals\n";
			}
			$totalBp += $readLength * 2;
		}
	}
	close( OUT );
}

sub generatePerfectReadsPairedParallel
{
	croak "Usage: generatePerfectReadsPairedParallel reference_fasta read_length total_bases insert_size output_fastq chunk_size" unless @_ == 6; 
	my $reference_fasta = shift;
	my $readLength = shift;
	my $total_bases = shift;
	my $insert_size = shift;
	my $output = shift;
	my $chunkSize = shift;
	
	croak "Chunk size must be less than total bases!" unless $chunkSize < $total_bases;
	
	my $total = 0;
	my $i =0;
	my $rand = int( rand( 10000000 ) );
	while()
	{
		if( $total_bases - $total > $chunkSize )
		{
			my $cmd = 'bsub -q basement -P hpag-pipeline -o o -e e -J generate.'.$rand.'.'.$i.' perl -w -e "use ReadSimulator;ReadSimulator::generatePerfectReadsPaired( \''.$reference_fasta.'\', \''.$readLength.'\', \''.$chunkSize.'\', \''.$insert_size.'\',\'chunk.'.$rand.'.'.$i.'.fastq\');"';
			system( $cmd );
			print $cmd."\n";
			$total += $chunkSize;
		}
		else
		{
			my $lastChunk = $total_bases - $total;
			my $cmd = 'bsub -q basement -P hpag-pipeline -o o -e e -J generate.'.$rand.'.'.$i.' perl -w -e "use ReadSimulator;ReadSimulator::generatePerfectReadsPaired( \''.$reference_fasta.'\', \''.$readLength.'\', \''.$lastChunk.'\', \''.$insert_size.'\',\'chunk.'.$rand.'.'.$i.'.fastq\');"';
			system( $cmd );
			print $cmd."\n";
			last;
		}
		$i ++;
	}
	
	if( $i > 0 )
	{
		my $cmd = 'bsub -o o -e e -w "done(generate.'.$rand.'.*)" -q basement -P hpag-pipeline -o o -e e "cat chunk.'.$rand.'.*.fastq > '.$output.'; rm chunk.'.$rand.'.*.fastq"';
		system( $cmd );
		print $cmd."\n";
	}
}
1;
