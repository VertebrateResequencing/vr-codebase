package IndexFileUtilities;
use strict;
use Carp;
use File::Basename;
use File::Spec;
use Digest::MD5 'md5_hex';

=pod
	A function that takes an index file and replaces all the fastq names with the suffix: recal.fastq.gz
=cut
sub index2RecalIndexFile
{
	croak "Usage: index2RecalIndexFile pre_index_file post_index_file" unless @_ == 2;
	
	my $old = shift;
	my $new = shift;
	
	croak "Cant find pre index file!\n" unless -f $old;
	
	open( IN, $old ) or die "Cannot open index file\n";
	open( OUT, ">$new" ) or die "Cannot create new index file\n";
	while( <IN> )
	{
		chomp;
		if( $_ =~ /^FASTQ_FILE/ )
		{
			print OUT $_."\n";
			next;
		}
		
		my @s = split( /\t/, $_ );
		
		$s[ 0 ] =~ s/(.*)\.fastq\.gz/$1.recal.fastq.gz/;
		if( $_ =~ /PAIRED/ && defined( $s[ 19 ] ) )
		{
			$s[ 19 ] =~ s/(.*)\.fastq\.gz/$1.recal.fastq.gz/;
		}
		
		foreach( @s ){print OUT $_."\t";}
		print OUT "\n";
	}
	close( IN );
	close( OUT );
}

=pod
a function to uniquify the SRR numbers in an index file (needed for early index files with multiple SRR per run!)
=cut
sub uniquifySRRindexFile
{
	croak "Usage: uniquifySRRindexFile in_index_file out_index_file" unless @_ == 2;
	my $pre = shift;
	my $post = shift;
	
	croak "Cant find input index file\n" unless -f $pre;
	
	open( PRE, "$pre" ) or die "Cannot open input index file\n";
	open( POST, ">$post" ) or die "Cannot create output file: $!\n";
	my %srr; #already defined SRR's
	my %fastqSRR; #fastq file to SRR mapping
	while( <PRE> )
	{
		chomp;
		next if( $_ =~ /^\s+$/ || length( $_ ) == 0 || $_ =~ /^FASTQ_FILE/ );
		
		my @s = split( /\t+/, $_ );
		
		if( $_ =~ /PAIRED/ && defined( $s[ 19 ] ) )
		{
			if( defined( $fastqSRR{ $s[ 19 ] } ) ) #pair has a SRR assigned
			{
				print POST $s[ 0 ]."\t".$s[ 1 ]."\t".$fastqSRR{ $s[ 19 ] };
				for( my $i = 3; $i < @s; $i ++ )
				{
					print POST  "\t".$s[ $i ];
				}
				print POST "\n";
				next;
			}
			
			if( defined( $srr{ $s[ 2 ] } ) ) #SRR has already been assigned
			{
				my $i = 1;
				while( defined( $srr{ $s[ 2 ].'_'.$i } ) )
				{
					$i ++;
				}
				
				$srr{ $s[ 2 ].'_'.$i } = 1;
				$fastqSRR{ $s[ 0 ] } = $s[ 2 ].'_'.$i;
				
				print POST $s[ 0 ]."\t".$s[ 1 ]."\t".$s[ 2 ].'_'.$i;
				for( my $i = 3; $i < @s; $i ++ )
				{
					print POST  "\t".$s[ $i ];
				}
				print POST "\n";
			}
			else
			{
				print POST $_."\n";
				$srr{ $s[ 2 ] } = 1;
				$fastqSRR{ $s[ 0 ] } = $s[ 2 ];
			}
		}
		else
		{
			if( defined $srr{ $s[ 2 ] } )
			{
				my $i = 1;
				while( defined( $srr{ $s[ 2 ].'_'.$i } ) )
				{
					$i ++;
				}
				$srr{ $s[ 2 ].'_'.$i } = 1;
				
				print POST $s[ 0 ]."\t".$s[ 1 ]."\t".$s[ 2 ].'_'.$i;
				for( my $i = 3; $i < @s; $i ++ )
				{
					print POST "\t".$s[ $i ];
				}
				print POST "\n";
			}
			else
			{
				print POST $_."\n";
				$srr{ $s[ 2 ] } = 1;
			}
		}
	}
	close( PRE );
	close( POST );
}

=pod
	A function to take the output of md5sum run on a set of fastq's and check the md5 vs. a meta index file
=cut
sub verifyMd5
{
	croak "Usage: verifyMd5 index_file md5sum_output" unless @_ == 2;
	my $index = shift;
	my $md5sum = shift;
	
	croak "Cant find input index file\n" unless -f $index;
	
	my %md5;
	open( MD5, $md5sum ) or die "Cant open md5sum output\n";
	while( <MD5> )
	{
		chomp;
		$_ =~ /(.*)\s+(.*\.fastq.gz)$/;

		next if( ! defined $1 );
		
		$md5{ basename( $2 ) } = $1;
		$md5{ basename( $2 ) } =~ s/^\s+//;
		$md5{ basename( $2 ) } =~ s/\s+$//;
	}
	close( MD5 );
	
	open( IN, "$index" ) or die "Cannot open input index file\n";
	while( <IN> )
	{
		chomp;
		my @s = split( /\t/, $_ );
		
		my $fastq = basename( $s[ 0 ] );
		my $md5 = $s[ 1 ];
		$md5 =~ s/^\s+//;
		$md5 =~ s/\s+$//;
		
		if( ! defined( $md5{ $fastq } ) )
		{
			print "ERROR: Cant find fastq file in md5sum output: $fastq\n";
			next;
		}
		
		if( $md5{ $fastq } ne $md5 )
		{
			print "ERROR: MD5 incorrect for: $fastq\n";
			print $md5{ $fastq }." vs. $md5\n";
			exit;
		}
	}
	close( IN );
	
	print "success\n";
}

=pod
	A function that takes a directory of recalibrated files, uncalibated directory, and an index file. It prints which lines of the index file are not recalibrated.
=cut
sub reportUncalibrated
{
	croak "Usage: reportUncalibrated directoryCalibrated indexFile uncalDirectory" unless @_ == 3;
	my $dir = shift;
	my $indexFile = shift;
	my $uncalibratedDir = shift;
	
	my %index;
	open( IN, $indexFile ) or die "Cant open index file\n";
	while( <IN> )
	{
		chomp;
		my @s = split( /\t/, $_ );
		
		$s[ 0 ] =~ /(.*)\.(fastq)\.(gz)/;
		
		#print "Searching: $dir/$1.recal.fastq.gz\n";
		
		if( -f "$uncalibratedDir/$s[ 0 ]" && ! -f "$dir/$1.recal.fastq.gz" )
		{
			print $_."\n";
		}
	}
	close( IN );
}

=pod
	A function that takes an index file and a directory of fastq - prints out the entries for the fastq files
=cut
sub extractFastqEntries
{
	croak "Usage: extractEntries directory indexFile output_indexF" unless @_ == 3;
	my $dir = shift;
	my $indexFile = shift;
	my $outputIndexFile = shift;
	
	croak "Cant find index file: $indexFile\n" unless -f $indexFile;
	croak "Cant find fastq directory: $dir\n" unless -d $dir;
	
	open( IN, $indexFile ) or die "Cant open index file\n";
	open( OUT, ">$outputIndexFile" ) or die "Cant create output index file\n";
	while( <IN> )
	{
		chomp;
		my @s = split( /\t/, $_ );
		
		$s[ 0 ] =~ /(.*)\.(fastq)\.(gz)/;
		
		#print "Searching: $dir/$1.recal.fastq.gz\n";
		
		if( -f "$dir/$s[ 0 ]" )
		{
			print OUT $_."\n";
		}
	}
	close( IN );
	close( OUT );
}

sub indexFilesIntersectionDifference
{
	croak "Usage: indexFilesIntersectionDifference indexFile1 indexFile2 output_intersection output_differences" unless @_ == 4;
	my $indexFile1 = shift;
	my $indexFile2 = shift;
	my $intersection = shift;
	my $differences = shift;
	
	croak "Cant find index file: $indexFile1\n" unless -f $indexFile1;
	croak "Cant find index file: $indexFile2\n" unless -f $indexFile2;
	
	open( IN1, $indexFile1 ) or die "Cannot open index file1\n";
	open( IN2, $indexFile2 ) or die "Cannot open index file1\n";
	my %index1;
	my %index2;
	while( <IN1> )
	{
		chomp;
		my @s = split( /\t/, $_ );
		
		if( $s[ 0 ] =~ /.*\.fastq.gz/ )
		{
			$index1{ $s[ 0 ] } = $_;
		}
	}
	close( IN1 );
	
	while( <IN2> )
	{
		chomp;
		my @s = split( /\t/, $_ );
		
		if( $s[ 0 ] =~ /.*\.fastq.gz/ )
		{
			$index2{ $s[ 0 ] } = $_;
		}
	}
	close( IN2 );
	
	my @index1Keys = keys( %index1 );
	my @index2Keys = keys( %index2 );
	
	my @isect = ();
	my @diff = ();
	my @union = ();
	my %count;
	my $e;
	
	foreach $e (@index1Keys, @index2Keys) { $count{$e}++ }
	
	foreach $e (keys %count)
	{
		push( @union, $e );
		push @{ $count{$e} == 2 ? \@isect : \@diff }, $e;
	}
	
	open( ISECT, ">$intersection" ) or die "Cant create intersection file\n";
	foreach( @isect )
	{
		print ISECT "$index1{ $_ }\n";
	}
	close( ISECT );
	
	open( DIFF, ">$differences" ) or die "Cant create intersection file\n";
	foreach( @diff )
	{
		print DIFF "$index1{ $_ }\n" unless ! defined $index1{ $_ };
		print DIFF "$index2{ $_ }\n" unless ! defined $index2{ $_ };
	}
	close( DIFF );
}

#for downloading fastq files from an index file - with a command such as wget or scp etc. that dont require user interaction
sub downloadFastqFiles
{
	croak "Usage: downloadFastqFiles indexFile destinationDirectory downloadCommandPrefix" unless @_ == 3;
	my $indexFile = shift;
	my $destinationDir = shift;
	my $downloadCmdPrefix = shift;
	
	croak "Cant find index file: $indexFile\n" unless -f $indexFile;
	
	open( INDEX, $indexFile ) or die "Cant open index file: $!\n";
	
	chdir( $destinationDir );
	
	while( <INDEX> )
	{
		chomp;
		my @s = split( /\t/, $_ );
		
		if( $s[ 0 ] =~ /.*\.fastq.gz/ )
		{
			if( ! -f basename( $s[ 0 ] ) )
			{
				print "$downloadCmdPrefix/$s[ 0 ]\n";
				system( "$downloadCmdPrefix/$s[ 0 ]" );
			}
			else
			{
				print "Already downloaded: ".basename( $s[ 0 ] )."\n";
			}
		}
	}
	close( INDEX );
}

1;
