package AssemblyTools;
use strict;
use Carp;
use Cwd;
use File::Basename;

#a perl module of many commonly used functions during assembly analysis

sub extractSubsetFasta
{
	croak "Usage: extractSubsetFasta seq_list fasta_file_name output_file_name" unless @_ == 3;
	my $list = shift;
	my $file = shift;
	my $output_file = shift;
	
	if( ! (-f $file) && ! -l $file ){croak "Cannot find file: $file\n";}
	if( ! (-s $file) ){croak "Empty file: $file\n";}
	
	my %extract;
	if( -f $list )
	{
		open( REMOVE, $list ) or die "Cannot open reads/contigs to be extracted file\n";
		while( <REMOVE> )
		{
			chomp;
			$extract{ $_ } = 1;
		}
		close( REMOVE );
	}
	else
	{
		my @split = split( /,/, $list );
		foreach( @split )
		{
			$extract{ $_ } = 1;
		}
	}
	
	if( $file =~ /\.gz$/ )
	{
		open( SEQS, "gunzip -c $file |" ) or die "Cannot open gzipped fastq file\n";
	}
	else
	{
		open( SEQS, $file ) or die "Failed to open reads file";
	}
	
	open( OUT, ">$output_file" ) or die "Cannot create output file\n";
	my $currentSeqName = '';
	my $currentSeq = '';
	while( <SEQS> )
	{
		chomp;
		if( $_ =~ /^>/ )
		{
			if( length( $currentSeqName ) > 0 )
			{
				if( defined( $extract{ $currentSeqName } ) )
				{
					print OUT $currentSeq."\n";
				}
			}
			
			$currentSeq = $_."\n";
			$currentSeqName = substr( ( split( /\s+/, $_ ) )[ 0 ], 1 );
		}
		else
		{
			$currentSeq .= $_;
		}
	}
	
	if( length( $currentSeqName ) > 0 )
	{
		if( defined( $extract{ $currentSeqName } ) )
		{
			print OUT $currentSeq."\n";
		}
	}
	close( SEQS );
	close( OUT );
}

#takes either a comma separated list OR file of read names
sub extractSubsetFastq
{
	croak "Usage: extractSubsetReadsFastq reads_list fastq_file_name output_file_name" unless @_ == 3;
	my $list = shift;
	my $file = shift;
	my $output_file = shift;
	
	if( ! (-f $file) && ! -l $file ){croak "Cannot find file: $file\n";}
	if( ! (-s $file) ){croak "Empty file: $file\n";}
	
	my %extract;
	if( -f $list )
	{
		open( REMOVE, $list ) or die "Cannot open reads/contigs to be extracted file\n";
		while( <REMOVE> )
		{
			chomp;
			$extract{ $_ } = 1;
		}
		close( REMOVE );
	}
	else
	{
		my @split = split( /,/, $list );
		foreach( @split )
		{
			$extract{ $_ } = 1;
		}
	}
	
	my %reads_file;
	if( $file =~ /\.gz$/ )
	{
		open( READS, "gunzip -c $file |" ) or die "Cannot open gzipped fastq file\n";
	}
	else
	{
		open( READS, $file ) or die "Failed to open reads file";
	}
	
	open( OUT, ">$output_file" ) or die "Cannot create output file\n";
	while( <READS> )
	{
		chomp;
		my $key = substr( $_, 1 ); #remove the @ sign
		my $sequence = <READS>;
		chomp( $sequence );
		my $qual_name = <READS>;
		chomp( $qual_name );
		my $quals = <READS>;
		chomp( $quals );
		
		if( length( $sequence ) != length( $quals ) )
		{
			print "WARNING: Number of quals not equal number of bases: $key\n";
			exit;
		}
		
		if( defined $extract{ $key } )
		{
			print OUT "$_\n"."$sequence\n"."$qual_name\n"."$quals\n";
		}
	}
	close( READS );
	close( OUT );
}

#takes either a comma separated list OR file of read names
sub excludeSubsetReadsFastq
{
	croak "Usage: excludeSubsetReadsFastq reads_list fastq_file_name output_file_name" unless @_ == 3;
	my $list = shift;
	my $file = shift;
	my $output_file = shift;
	
	if( ! (-f $file) && ! -l $file ){croak "Cannot find file: $file\n";}
	if( ! (-s $file) ){croak "Empty file: $file\n";}
	
	my %extract;
	if( -f $list )
	{
		open( EX, $list ) or die "Cannot open reads/contigs to be extracted file\n";
		while( <EX> )
		{
			chomp;
			$extract{ $_ } = 1;
		}
		close( EX );
	}
	else
	{
		my @split = split( /,/, $list );
		foreach( @split )
		{
			$extract{ $_ } = 1;
		}
	}
	
	my %reads_file;
	if( $file =~ /\.gz$/ )
	{
		open( READS, "gunzip -c $file |" ) or die "Cannot open gzipped fastq file\n";
	}
	else
	{
		open( READS, $file ) or die "Failed to open reads file";
	}
	
	open( OUT, ">$output_file" ) or die "Cannot create output file\n";
	while( <READS> )
	{
		chomp;
		my $key = substr( $_, 1 ); #remove the @ sign
		my $sequence = <READS>;
		chomp( $sequence );
		my $qual_name = <READS>;
		chomp( $qual_name );
		my $quals = <READS>;
		chomp( $quals );
		
		if( length( $sequence ) != length( $quals ) )
		{
			print "WARNING: Number of quals not equal number of bases: $key\n";
			exit;
		}
		
		if( ! defined $extract{ $key } )
		{
			print OUT "$_\n"."$sequence\n"."$qual_name\n"."$quals\n";
		}
	}
	close( READS );
	close( OUT );
}

sub printSingleHitReadsCigar
{
	croak "Usage: printSingleHitReadsCigar cigar.in cigar.out\n" unless @_ == 2;
	my $in = shift;
	my $out = shift;
	
	croak "Cant find cigar file\n" unless -f $in;
	
	my $currentRead = '';
	my $currentHit = '';
	my $numHits = 0;
	my $singles = 0;
	
	if( $in =~ /\.gz$/ )
	{
		open( IN, "gunzip -c $in |" ) or die "Cannot open gzipped cigar file\n";
	}
	else
	{
		open( IN, $in ) or die "cannot open cigar file\n";
	}
	open( CIGAR, ">$out" ) or die "cant create cigar output\n";
	while( <IN> )
	{
		if( $_ =~ /^cigar::/ )
		{
			my @s = split( /\s+/, $_ );
			if( length( $currentRead ) == 0 )
			{
				$currentRead = $s[ 1 ];
				$currentHit = $_;
				$numHits ++;
			}
			elsif( $currentRead ne $s[ 1 ] )
			{
				if( $numHits == 1 )
				{
					print CIGAR $currentHit;
				}
				
				$currentRead = $s[ 1 ];
				$currentHit = $_;
				$numHits = 1;
			}
			else
			{
				$numHits ++;
			}
		}
	}
	close( IN );
	close( CIGAR );
}

sub splitFastaChunks
{
	croak "Usage: splitFastaChunks fasta size_of_chunks_in_bytes output_files_prefix output_directory" unless @_ == 4;
	my $fasta = shift;
	my $size = shift;
	my $prefix = shift;
	my $dir = shift;
	
	croak "Cant find fasta file\n" unless -f $fasta;
	croak "Chunk size must be numerical\n" unless $size =~ /^\d+$/;
	
	if( ! -d $dir )
	{
		mkdir( $dir );
	}
	
	my $bytes = 0;
	my $currentChunk = '';
	my $fileCounter = 1;
	open( IN, $fasta ) or die "Cannot open fasta file\n";
	while( <IN> )
	{
		chomp;
		if( $_ =~ /^>/ )
		{
			if( $bytes > $size )
			{
				open( OUT, ">$dir/$prefix.$fileCounter.fasta" ) or die "Cannot create output file\n";
				print OUT $currentChunk."\n";
				close( OUT );
				$fileCounter ++;
				
				$currentChunk = $_."\n";
				$bytes = 0;
			}
			else
			{
				$bytes += length( $_ );
				$currentChunk .= $_."\n";
			}
		}
		else
		{
			$bytes += length( $_ );
			$currentChunk .= $_."\n";
		}
	}
	close( IN );
}

sub fastaHash
{
	croak "Usage: fastaHash reference_fasta\n" unless @_ == 1;
	my $fasta = shift;
	
	my %reference;
	my $currentSequence = '';
	my $currentName = '';
	open( IN, "$fasta" ) or die "Cannot open $fasta\n";
	while( <IN> )
	{
		chomp;
		if( $_ =~ /^>/ )
		{
			if( length( $currentSequence ) > 0 )
			{
				$reference{ $currentName } = $currentSequence;
			}
			
			$currentName = substr( (split( /\s+/, $_ ) )[ 0 ], 1 );
			$currentSequence = $_."\n";
		}
		else
		{
			$currentSequence .= $_;
		}
	}
	close( IN );
	
	if( length( $currentSequence ) > 0 )
	{
		$reference{ $currentName } = $currentSequence;
	}
	
	return \%reference;
}

=pod

=head1 NAME

createFastqHash - 

=head1 SYNOPSIS

createFastqHash file_name

=head1 ARGUMENTS

The function takes the following arguments:

B<file_name> : the fastq file of reads

The function returns a hash of the fastq file - keys are the read names.

=head1 DESCRIPTION

A subroutine to create a hash table from a fastq file keys are the read names
entries are the sequences and quality values

=head1 AUTHOR

Thomas Keane I<tk2@sanger.ac.uk>

=cut

sub createFastqHash
{
	croak "Usage: createFastqHash file_name" unless @_ == 1; 
	my $file = shift;
	
	if( ! (-f $file) ){croak "Cannot find file: $file\n";}
	if( ! (-s $file) ){croak "Empty file: $file\n";}
	
	my %reads_file;
	if( $file =~ /\.gz$/ )
	{
		open( READS, "gunzip -c $file |" ) or die "Cannot open gzipped fastq file\n";
	}
	else
	{
		open( READS, $file ) or die "Failed to open reads file";
	}
	
	while( <READS> )
	{
		chomp;
		my @t = split( /\s+/, $_ );
		my $key = substr( $t[ 0 ], 1, length( $t[ 0 ] ) - 1 ); #remove the @ sign
		my $sequence = <READS>;
		chomp( $sequence );
		my $qual_name = <READS>;
		chomp( $qual_name );
		my $quals = <READS>;
		chomp( $quals );
		
		croak "Number of quals not equal to number of bases: $key\n" unless length( $quals ) == length( $sequence );
		
		$reads_file{ $key } = "$_\n"."$sequence\n"."$qual_name\n"."$quals";
	}
	close( READS );
	
	return \%reads_file;
}

#a function to split the sanger solexa fastq's into a format compatible for Jared's assembler....and possibly other programs
#gzip input file = gzipped output files
sub sanger2SplitFastq
{
	croak "Usage: sanger2SplitFastq inputfastq fastq1name fastq2name" unless @_ == 3; 
	my $fastq = shift;
	my $fastq1 = shift;
	my $fastq2 = shift;
	
	croak "Cant find fastq file: $fastq\n" unless -f $fastq;
	
	if( $fastq =~ /.*\.fastq\.gz$/ )
	{
		open( L, "| gzip -c > $fastq1" ) or die "Cannot create $fastq1\n";
		open( R, "| gzip -c > $fastq2" ) or die "Cannot create $fastq2\n";
		open( IN, "gunzip -c $fastq |" ) or die "Cannot read input fastq: $fastq\n";
	}
	else
	{
		open( L, ">$fastq1" ) or die "Cannot create $1"."_1.fastq\n";
		open( R, ">$fastq2" ) or die "Cannot create $2"."_2.fastq\n";
		open( IN, "$fastq" ) or die "Cannot read input fastq: $fastq\n";
	}
	
	while( <IN> )
	{
		chomp;
		
		croak "Bad head line in fastq: $_" unless $_ =~ /^@.+/;
		
		my $r1n = $_."/1";
		my $r2n = $_."/2";
		
		print L "$r1n\n";
		print R "$r2n\n";
		
		my $t = <IN>;
		chomp( $t );
		print L substr( $t, 0, length( $t ) / 2 )."\n";
		print R substr( $t, length( $t ) / 2 )."\n";
		$t = <IN>;
		chomp( $t );
		
		croak "Bad qual header line in fastq: $t" unless $t =~ /^\+/;
		
		print L $t."\n";
		print R $t."\n";
		$t = <IN>;
		chomp( $t );
		print L substr( $t, 0, length( $t ) / 2 )."\n";
		print R substr( $t, length( $t ) / 2 )."\n";
	}
	close( IN );
	close( L );
	close( R );
	
	return 1;
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
	open( FASTA_QUAL_FILE, ">$output_file.qual" ) or die "Error: Cannot create fasta QUAL file\n";
	
	while( <READS_FILE> )
	{
		chomp;
		my $key = substr( $_, 1 ); #remove the @ sign
		my $sequence = <READS_FILE>;
		chomp( $sequence );
		my $qual_name = <READS_FILE>;
		chomp( $qual_name );
		my $quals = <READS_FILE>;
		chomp( $quals );
		
		print FASTA_FILE ">".$key."\n".$sequence."\n";
		
		my @q = split( //, $quals );
		my $fasta_quals = '';
		foreach( @q )
		{
			$fasta_quals .= ( unpack( 'C', $_ ) - 33 ).' ';
		}
		print FASTA_QUAL_FILE ">".$key."\n".$fasta_quals."\n";
	}
	close( READS_FILE );
	close( FASTA_FILE );
	close( FASTA_QUAL_FILE );
}

sub splitUnpairedFastq
{
	croak "Usage: splitUnpairedFastq fastq1 prefix chunkBases outputdirectory" unless @_ == 4; 
	my $fastq1 = shift;
	my $prefix = shift;
	my $chunkBases = shift;
	my $outputDirectory = shift;
	
	croak "Cant find output directory: $outputDirectory\n" unless -d $outputDirectory;
	
	if( $fastq1 =~ /\.gz$/ )
	{
		open( READS1, "gunzip -c $fastq1 |" ) or die "Cannot open gzipped fastq1 file\n";
	}
	else
	{
		open( READS1, "$fastq1" ) or die "Cannot open fastq1 file\n";
	}
	
	my $fileCount = 0;
	my $baseCount = 0;
	open( L, ">$outputDirectory/$prefix$fileCount.fastq" ) or die "Cannot create LEFT file\n";
	while( <READS1> )
	{
		my $rn1 = $_;
		print L $rn1;
		
		my $seq1 = <READS1>;
		print L $seq1;
		
		my $qn1 = <READS1>;
		print L $qn1;
		
		my $q1 = <READS1>;
		print L $q1;
		
		$baseCount += length( $seq1 ) - 1;
		if( $baseCount > $chunkBases )
		{
			close( L );
			
			$fileCount ++;
			$baseCount = 0;
			open( L, ">$outputDirectory/$prefix$fileCount.fastq" ) or die "Cannot create LEFT file\n";
		}
	}
	close( L );
	close( READS1 );
	
	return 1;
}

=pod
A function to split paired fastq's into chunks with a defined prefix followed by sequential numbering (i.e. <prefix><num>.fastq)
=cut

sub splitPairedFastq
{
	croak "Usage: splitPairedFastq fastq1 fastq2 prefix chunkBases outputdirectory" unless @_ == 5; 
	my $fastq1 = shift;
	my $fastq2 = shift;
	my $prefix = shift;
	my $chunkBases = shift;
	my $outputDirectory = shift;
	
	croak "Cant find output directory: $outputDirectory\n" unless -d $outputDirectory;
	
	if( $fastq1 =~ /\.gz$/ )
	{
		open( READS1, "gunzip -c $fastq1 |" ) or die "Cannot open gzipped fastq1 file\n";
		open( READS2, "gunzip -c $fastq2 |" ) or die "Cannot open gzipped fastq2 file\n";
	}
	else
	{
		open( READS1, "$fastq1" ) or die "Cannot open fastq1 file\n";
		open( READS2, "$fastq2" ) or die "Cannot open fastq2 file\n";
	}
	
	my $fileCount = 0;
	my $baseCount = 0;
	open( L, ">$outputDirectory/$prefix$fileCount"."_1.fastq" ) or die "Cannot create LEFT file\n";
	open( R, ">$outputDirectory/$prefix$fileCount"."_2.fastq" ) or die "Cannot create RIGHT file\n";
	while( <READS1> )
	{
		my $rn1 = $_;
		my $rn2 = <READS2>;
		my $seq1 = <READS1>;
		my $seq2 = <READS2>;
		my $qn1 = <READS1>;
		my $qn2 = <READS2>;
		my $q1 = <READS1>;
		my $q2 = <READS2>;
		
		if( $seq1 =~ /^A+\n$/ || $seq1 =~ /^N+\n$/ || $seq2 =~ /^A+\n$/ || $seq2 =~ /^N+\n$/ )
		{
			next;
		}
=pod		
		#MAQ requires /1 before readname
		if( $rn1 !~ /\@.*\/1/ )
		{
			chomp( $rn1 );
			$rn1 .= "/1\n";
		}
=cut
		print L $rn1;
		print L $seq1;
		print L $qn1;
		print L $q1;
=pod		
		#MAQ requires /2 before readname
		if( $rn2 !~ /\@.*\/2/ )
		{
			chomp( $rn2 );
			$rn2 .= "/2\n";
		}
=cut		
		print R $rn2;
		print R $seq2;
		print R $qn2;
		print R $q2;
		
		$baseCount += length( $seq1 ) - 1;
		$baseCount += length( $seq2 ) - 1;
		if( $baseCount > $chunkBases )
		{
			close( L );
			close( R );
			
			$fileCount ++;
			$baseCount = 0;
			open( L, ">$outputDirectory/$prefix$fileCount"."_1.fastq" ) or die "Cannot create LEFT file\n";
			open( R, ">$outputDirectory/$prefix$fileCount"."_2.fastq" ) or die "Cannot create RIGHT file\n";
		}
		
		if( eof( READS2 ) && ! eof( READS1 ) )
		{
			croak "Fastq2 file truncated: $fastq2\n";
		}
	}
	close( L );
	close( R );
	close( READS1 );
	close( READS2 );
	
	if( ! eof( READS2 ) )
	{
		croak "Fastq1 file truncated: $fastq1\n";
	}
	
	return 1;
}

sub printFofnDifferences
{
	croak "Usage: printFofnDifferences reads_fofn reads_fofn output_fofn" unless @_ == 3; 
	my $file1 = shift;
	my $file2 = shift;
	my $output = shift;
	
	if( (! (-f $file1)) || (! (-f $file2)) ){croak "Cannot find input file\n";}
	if( (! (-s $file1)) || (! (-s $file2)) ){croak "Empty input file\n";}
	
	open( FILE1, $file1 ) or die "Cannot open file: $file1\n";
	open( FILE2, $file2 ) or die "Cannot open file: $file2\n";
	open( OUTPUT, ">$output" ) or die "Cannot create output file\n";
	
	my @fofn1;
	while( <FILE1> )
	{
		chomp;
		push( @fofn1, $_ );
	}
	close( FILE1 );
	
	my @fofn2;
	while( <FILE2> )
	{
		chomp;
		push( @fofn2, $_ );
	}
	close( FILE2 );
	
	print "Read in fofn files...starting comparison\n";
	
	my @isect = ();
	my @diff = ();
	my @union = ();
	my %count;
	my $e;
	
	foreach $e (@fofn1, @fofn2) { $count{$e}++ }
	
	foreach $e (keys %count)
	{
		push( @union, $e );
		push @{ $count{$e} == 2 ? \@isect : \@diff }, $e;
	}
	
	my $insect_size = @isect;
	print "There are $insect_size common reads\n";
	
	foreach( @diff )
	{
		print OUTPUT "$_\n";
	}

	close( OUTPUT );
}

sub printFofnOverlaps
{
	croak "Usage: printFofnDifferences reads_fofn reads_fofn output_fofn" unless @_ == 3; 
	my $file1 = shift;
	my $file2 = shift;
	my $output = shift;
	
	if( (! (-f $file1)) || (! (-f $file2)) ){croak "Cannot find input file\n";}
	if( (! (-s $file1)) || (! (-s $file2)) ){croak "Empty input file\n";}
	
	open( FILE1, $file1 ) or die "Cannot open file: $file1\n";
	open( FILE2, $file2 ) or die "Cannot open file: $file2\n";
	open( OUTPUT, ">$output" ) or die "Cannot create output file\n";
	
	my %fofn1;
	while( <FILE1> )
	{
		chomp;
		$fofn1{ $_ } = 1;
	}
	close( FILE1 );
	
	my %fofn2;
	while( <FILE2> )
	{
		chomp;
		$fofn2{ $_ } = 1;
	}
	close( FILE2 );
	
	print "Read in fofn files...starting comparison\n";
	
	my %isect;
	my @diff1;
	my @diff2;
	foreach( keys( %fofn1 ) )
	{
		if( ! defined $fofn2{ $_ } )
		{
			push( @diff1, $_ );
		}
		else
		{
			$isect{ $_ } = 1;
		}
	}
	
	foreach( keys( %fofn2 ) )
	{
		if( ! defined $fofn1{ $_ } )
		{
			push( @diff2, $_ );
		}
		else
		{
			$isect{ $_ } = 1;
		}
	}
	
	foreach( @diff1 )
	{
		print OUTPUT "UNIQ1: $_\n";
	}
	
	foreach( @diff2 )
	{
		print OUTPUT "UNIQ2: $_\n";
	}
	
	foreach( keys( %isect ) )
	{
		print OUTPUT "ISECT: $_\n";
	}
}

sub printFofnIntersection
{
	croak "Usage: printFofnIntersection reads_fofn reads_fofn output_fofn" unless @_ == 3; 
	my $file1 = shift;
	my $file2 = shift;
	my $output = shift;
	
	if( (! (-f $file1)) || (! (-f $file2)) ){croak "Cannot find input file\n";}
	if( (! (-s $file1)) || (! (-s $file2)) ){croak "Empty input file\n";}
	
	open( FILE1, $file1 ) or die "Cannot open file: $file1\n";
	open( FILE2, $file2 ) or die "Cannot open file: $file2\n";
	open( OUTPUT, ">$output" ) or die "Cannot create output file\n";
	
	my @fofn1;
	while( <FILE1> )
	{
		chomp;
		push( @fofn1, $_ );
	}
	close( FILE1 );
	
	my @fofn2;
	while( <FILE2> )
	{
		chomp;
		push( @fofn2, $_ );
	}
	close( FILE2 );
	
	print "Read in fofn files...starting comparison\n";
	
	my @isect = ();
	my @diff = ();
	my @union = ();
	my %count;
	my $e;
	
	foreach $e (@fofn1, @fofn2) { $count{$e}++ }
	
	foreach $e (keys %count)
	{
		push( @union, $e );
		push @{ $count{$e} == 2 ? \@isect : \@diff }, $e;
	}
	
	my $insect_size = @isect;
	print "There are $insect_size common reads\n";
	
	foreach( @isect )
	{
		print OUTPUT "$_\n";
	}

	close( OUTPUT );
}

sub revCompDNA
{
	croak "Usage: revCompDNA string\n" unless @_ == 1;
	my $seq = shift;
	$seq= uc( $seq );
	$seq=reverse( $seq );
	
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return $seq;
}

sub writeClipPointMeta
{
	croak "Usage: writeClipPointMeta reads_fastq meta_file read_num" unless @_ == 3;
	my $fastq = shift;
	my $meta_file = shift;
	my $read_num = shift;
	
	croak "Cant find reads file: $fastq\n" unless -f $fastq;
	croak "Cant find meta file: $meta_file\n" unless -f $meta_file;
	croak "Invalid read number: $read_num\n" unless $read_num > -1 && $read_num < 3;
	
	my $clip = determineClipPointMaq( $fastq );
	
	system( qq[ echo "clip$read_num:$clip" >> $meta_file ] );
}

sub determineClipPointMaq
{
	croak "Usage: determineClipPointMaq reads_fastq" unless @_ == 1;
	my $fastq = shift;
	
	croak "Cant find reads file: $fastq\n" unless -f $fastq;
	
	if( $fastq =~ /\.gz$/ )
	{
		open( READS, "gunzip -c $fastq |" ) or die "Cannot open gzipped fastq1 file\n";
	}
	else
	{
		open( READS, "$fastq" ) or die "Cannot open fastq1 file\n";
	}
	
	my $totalReads = 0;
	my %clipLengths;
	while( <READS> )
	{
		chomp;
		my $name = $_;
		my $seq = <READS>;
		chomp( $seq );
		
		if( $seq =~ /.*N.*/ )
		{
			my @s = split( //, $seq );
			my $nCount = 0;
			my $clipLength = @s;
			for( my $i = 0; $i < @s; $i ++ )
			{
				if( $s[ $i ] eq 'N' )
				{
					$nCount ++;
					if( $nCount == 2 )
					{
						$clipLength = $i;
						if( ! defined( $clipLengths{ $clipLength } ) )
						{
							$clipLengths{ $clipLength } = 1;
						}
						else
						{
							$clipLengths{ $clipLength } ++;
						}
						last;
					}
				}
			}
		}
		<READS>;
		<READS>;
		$totalReads ++;
	}
	close( READS );
	
	#determine the point where 90% of the reads are ok
	my $numReads = 0;
	foreach( sort { $a <=> $b }( keys( %clipLengths ) ) )
	{
		$numReads += $clipLengths{ $_ };
		
		if( $numReads > $totalReads * 0.9 )
		{
			print ($_ - 1);
			return ($_ - 1);
		}
	}
	print "-1";
	return -1;
}

sub extractContigSection
{
	croak "Usage: extractContigSection input_file output_file contig_name start_index end_index" unless @_ == 5;
	my $input = shift;
	my $output = shift;
	my $contigName = shift;
	my $start = shift;
	my $finish = shift;
	
	if( ! (-f $input) ){croak "Cannot find file: $input\n";}
	if( ! (-s $input) ){croak "Empty file: $input\n";}
	
	if( $start > $finish )
	{
		croak "Start index must be less than end index Start:$start End:$finish\n";
	}
	elsif( $start < 0 )
	{
		croak "Start index must be > 0 Start:$start End:$finish\n";
	}
	
	open( CONTIG, $input ) or die "Cannot open contig file\n";
	my @contig = ();
	my $name = <CONTIG>;
	my $sequence = "";
	my $record = 0;
	while( <CONTIG> )
	{
		chomp;
		if( $_ =~ /^>.*/ )
		{
			if( $record == 1 )
			{
				last;
			}
			
			my $name = (split( /\s+/, $_ ))[ 0 ];
			if( $name =~ /^>$contigName/ )
			{
				$record = 1;
			}
		}
		
		if( $record == 1 )
		{
			$sequence .= $_;
		}
	}
	
	if( $start > length( $sequence ) ) #something wrong if start is off end!
	{
		my $length = length( $sequence );
		print "WARNING Indexes are outside of sequence range. Start:$start End:$finish Length:$length \n";
		
		#return a section that goes up to the end of the contig
		my $size = $finish - $start;
		my $section = substr( $contig[ 1 ], $length - $size, length( $sequence ) - 1 );
		open( OUTPUT, ">$output" ) or die "Cannot open new contig file\n";
		print OUTPUT ">section\n$section";
		close( OUTPUT );
		return;
	}
	
	if(  $finish > length( $sequence ) )
	{
		$finish = length( $sequence ) - 1; #go up to the end of the contig
	}
	
	#extract the section
	my $section = substr( $sequence, $start, $finish - $start );
	my $length_ = length( $section );
	print "Extracted $length_ bases\n";
	open( OUTPUT, ">$output" ) or die "Cannot open new contig file\n";
	print OUTPUT ">section\n$section";
	
	close( CONTIG );
	close( OUTPUT );
}

sub verifyFastqFile
{
	croak "Usage: verifyFastqFile fastq" unless @_ == 1;
	my $input = shift;
	
	if( $input =~ /\.gz$/ )
	{
		open( IN, "gunzip -c $input |" ) or die "Cant open input gzip file: $!\n";
	}
	else
	{
		open( IN, $input ) or die "Cant open input file: $!\n";
	}
	
	while( <IN> )
	{
		chomp;
		my $name = $_;
		croak "Invalid name: $name\n" unless $name =~ /^@.+/;
		
		my $seq = <IN>;
		chomp( $seq );
		
		croak "Invalid sequence: $seq\n" unless uc( $seq ) =~ /^A|C|G|T|N$/;
		
		my $qname = <IN>;
		chomp( $qname );
		
		croak "Invalid qual tag: $qname\n" unless $qname =~ /^\+.*/;
		
		my $quals = <IN>;
		chomp( $quals );
		
		croak "No. quals not equal to length sequence: $name ".length( $seq )." ".length( $quals )."\n" unless length( $seq ) == length( $quals );
	}
	close( IN );
	
	print "$input Good\n";
}

=pod
	A function to take the output of two md5sum runs on a set of files and check the md5 vs. a meta index file
=cut
sub verifyMd5Outputs
{
	croak "Usage: verifyMd5 md5_output md5_output" unless @_ == 2;
	my $md5sum1 = shift;
	my $md5sum2 = shift;
	
	croak "Cant find input md5sum file\n" unless -f $md5sum1 && -f $md5sum2;
	
	my %md5;
	open( MD5, $md5sum1 ) or die "Cant open md5sum output\n";
	while( <MD5> )
	{
		chomp;
		$_ =~ /(.*)\s+(.*)/;
		
		$md5{ $1 } = basename( $2 );
	}
	close( MD5 );
	
	open( MD5, $md5sum2 ) or die "Cant open md5sum output\n";
	while( <MD5> )
	{
		chomp;
		
		$_ =~ /(.*)\s+(.*)/;
		my $md5_ = $1;
		my $file = basename( $2 );
		
		if( ! defined( $md5{ $md5_ } ) )
		{
			print "Not defined in md5_1 file: $2\n";
			next;
		}
		
		if( $md5{ $md5_ } !~ /$file$/ && $file !~ /$md5{ $md5_ }$/ )
		{
			print "Incorrect md5: $file $md5{ $md5_ }\n";
		}
		else
		{
			#print "correct: ".$_."\n";
		}
	}
	print "success\n";
}

sub filterOutShortReads
{
	croak "Usage: filterOutShortReads fastq minLength outputFile" unless @_ == 3;
	my $inFastq = shift;
	my $minLen = shift;
	my $outFastq = shift;
	
	if( $inFastq =~ /\.gz$/ )
	{
		open( IN, "gunzip -c $inFastq |" ) or die "Cant open input gzip file: $!\n";
		open( OUT, "| gzip -c > $outFastq" ) or die "Cant create output fastq: $!\n";
	}
	else
	{
		open( IN, $inFastq ) or die "Cant open input file: $!\n";
		open( OUT, ">$outFastq" ) or die "Cant create output fastq: $!\n";
	}
	
	while( <IN> )
	{
		chomp;
		my $name = $_;
		
		my $seq = <IN>;
		chomp( $seq );
		
		my $qname = <IN>;
		chomp( $qname );
		
		my $quals = <IN>;
		chomp( $quals );
		
		if( length( $seq ) >= $minLen )
		{
			print OUT "$name\n$seq\n$qname\n$quals\n";
		}
	}
	close( OUT );
	close( IN );
}

#convert one or more sff files to a fastq file with phred-like qualities
sub sff2fastq
{
	croak "Usage: sff2fastq fastq_output sff1 [sff2 sff3.....]" unless @_ > 1;
	
	my $fastq = shift;
	
	foreach( @_ )
	{
		if( -f $_ )
		{
			print "Processing sff file: $_\n";
			#rescore the sff file using sfffile to phred-like values
			#system( "sfffile -r " );
			
			print "Running sffinfo...\n";
			#convert the sff file to fasta
			system( "sffinfo -seq $_ > /tmp/$$.fasta" );
			system( "sffinfo -qual $_ > /tmp/$$.qual" );
			
			#convert the fasta file to fastq
			fasta2fastq(  "/tmp/$$.fasta", "/tmp/$$.fastq" );
			
			system( "cat /tmp/$$.fastq >> $fastq" );
		}
	}
	
	system( "rm /tmp/$$.fast* /tmp/$$.qual" );
}

sub fasta2fastq
{
	croak "Usage: fasta2fastq fasta_file output_fastq\n" unless @_ == 2;
	
	my $r = createFastaHash( $_[ 0 ] );
	my %fasta_seqs = %{ $r };
	my @fasta_seqs_keys = keys( %fasta_seqs );
	
	my %fasta_quals;
	my @fasta_quals_keys;
	$r = undef;
	
	my $prefix = (split( /\./, $_[0]))[0];
	
	if( -f  $_[ 0 ].'.qual' )
	{
		$r = createFastaQualHash( $_[ 0 ].'.qual' );
		%fasta_quals = %{ $r };
		@fasta_quals_keys = keys( %fasta_quals );
		
		if( @fasta_seqs_keys != @fasta_quals_keys )
		{
			croak "unequal sequences vs. quals!\n";
		}
	}
	elsif( -f $prefix.'.qual')
	{
		$r = createFastaQualHash( $prefix.'.qual' );
		%fasta_quals = %{ $r };
		@fasta_quals_keys = keys( %fasta_quals );
		
		if( @fasta_seqs_keys != @fasta_quals_keys )
		{
			croak "unequal sequences vs. quals!\n";
		}
	}
	else
	{
		print "No qual file supplied! Faking qualities\n";
	}
	
	open( OUT, ">$_[1]" ) or die "cannot create output fastq\n";
	foreach( @fasta_seqs_keys )
	{
		my @s1 = split( /\n/, $fasta_seqs{ $_ } );
		
		my @quals;
		my @s;
		if( @fasta_quals_keys > 0 )
		{
			@s = split( /\n/, $fasta_quals{ $_ } );
			@quals = split( /\s+/, $s[ 1 ] );
			if( @quals != length( $s1[ 1 ] ) )
			{
				my $n = @quals;
				croak "Read: $s[ 0 ] Bases: ".length( $s1[ 1 ] )." Quals: $n\n";
			}
		}
		
		my $name = ( split( /\s+/, substr( $fasta_seqs{ $_ }, 1 ) ) )[ 0 ];
		my $seq = ( split( /\n/, substr( $fasta_seqs{ $_ }, 1 ) ) )[ 1 ];
		print OUT "@".$name."\n".$seq."\n";
		print OUT "+".$name."\n";
		print OUT "!";
		
		if( @fasta_quals_keys > 0 )
		{
			for( my $i = 1; $i < @quals; $i ++ )
			{
				print OUT chr($quals[ $i ] + 33);
			}
		}
		else
		{
			for( my $i = 1; $i < length( $s1[1] ); $i ++ )
			{
				print OUT "I";
				#print OUT "&";
			}
		}
		print OUT "\n";
	}
	close( OUT );
}

sub createFastaQualHash
{
	croak "Usage: createFastaQualHash file_name" unless @_ == 1; 
	my $file = shift;
	
	if( ! ( -f $file ) )
	{
		croak "Cannot find file: $file\n";
	}
	
	my %reads_file;
	if( $file =~ /\.gz$/ )
	{
		open( READS, "gunzip -c $file |" ) or die "Cannot open gzipped fastq file\n";
	}
	else
	{
		open( READS, $file ) or die "Failed to open reads file";
	}
	
	my $read_name = <READS>; #first readname in first line of file
	chomp( $read_name );
	my $read = "$read_name\n";
	$read_name =~ s/^\s+//;
	$read_name =~ s/\s+$//;
	my $line = <READS>; #read next line
	chomp( $line );
	while() #read down file until hit start of next read
	{
		if( $line =~ /^>.*/ ) #if hit start of next read
		{
			$read_name = substr( $read_name, 1, length( $read_name ) - 1 ); #remove the > sign
			$read_name =~ s/^\s+//;
			$read_name =~ s/\s+$//;
			
			if( defined $reads_file{ $read_name } )
			{
				print "Warning: Read entry already exists in hash: $read_name\n";
			}
			else
			{
				$reads_file{ $read_name } = $read; #enter into the hash
			}
			
			$read_name = $line; #next read name is in the line variable
			$read = "$read_name\n";
		}
		else
		{
			if( length( $read ) > 0 )
			{
				$line =~ s/^\s+//;
				$line =~ s/\s+$//;
				if( $read !~ /\n$/ )
				{
					$read = $read." ".$line; #add to info for the current read
				}
				else
				{
					$read = $read.$line;
				}
			}
			else
			{
				$read = $line;
			}
		}
		
		$line = <READS>;
		if( ! defined $line )
		{
			last;
		}
		chomp( $line );
	}
	close( READS );
	
	#enter final value into the hash
	$read_name = substr( $read_name, 1, length( $read_name ) - 1 ); #remove the > sign
	$read_name =~ s/^\s+//;
	$read_name =~ s/\s+$//;
	$reads_file{ $read_name } = $read; #enter into the hash
	
	return \%reads_file;
}

sub createFastaHash
{
	croak "Usage: createFastaHash file_name" unless @_ == 1; 
	my $file = shift;
	
	if( ! ( -f $file ) )
	{
		croak "Cannot find file: $file\n";
	}
	
	my %reads_file;
	if( $file =~ /\.gz$/ )
	{
		open( READS, "gunzip -c $file |" ) or die "Cannot open gzipped fastq file\n";
	}
	else
	{
		open( READS, $file ) or die "Failed to open reads file";
	}
	
	my $read_name = <READS>; #first readname in first line of file
	chomp( $read_name );
	my $read = "$read_name\n";
	$read_name =~ s/^\s+//;
	$read_name =~ s/\s+$//;
	my $line = <READS>; #read next line
	chomp( $line );
	while() #read down file until hit start of next read
	{
		if( $line =~ /^>.*/ ) #if hit start of next read
		{
			$read_name = substr( $read_name, 1, length( $read_name ) - 1 ); #remove the > sign
			$read_name =~ s/^\s+//;
			$read_name =~ s/\s+$//;
			
			if( defined $reads_file{ $read_name } )
			{
				print "Warning: Read entry already exists in hash: $read_name\n";
			}
			else
			{
				$reads_file{ $read_name } = $read; #enter into the hash
			}
			
			$read_name = $line; #next read name is in the line variable
			$read = "$read_name\n";
		}
		else
		{
			if( length( $read ) > 0 )
			{
				$read = $read.$line; #add to info for the current read
			}
			else
			{
				$read = $line;
			}
		}
		
		$line = <READS>;
		if( ! defined $line )
		{
			last;
		}
		chomp( $line );
	}
	close( READS );
	
	#enter final value into the hash
	$read_name = substr( $read_name, 1, length( $read_name ) - 1 ); #remove the > sign
	$reads_file{ $read_name } = $read; #enter into the hash
	
	return \%reads_file;
}

1;
