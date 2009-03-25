package G1KUtilities;
use strict;
use Carp;
use File::Spec;
use Cwd;

#a function to determin which ssaha disk hash table to use from the current location within the hierarchy
sub path2Gender
{
	croak "Usage: path2Gender gender_file" unless @_ == 1;
	my $genderFile = shift;
	
	my %genders;
	open( G, "$genderFile" ) or die "Cant open gender file: $genderFile\n";
	while( <G> )
	{
		my @s = split( /\s+/, $_ );
		$genders{ $s[ 0 ] } = $s[ 1 ];
	}
	close( G );
	
	my $cwd = getcwd;
	#print "cwd: $cwd\n";
	
	if( $cwd =~ /\/NA\d+/ )
	{
		my @s = split( /\//, $cwd );
		foreach( @s )
		{
			if( $_ =~ /^NA\d+$/ )
			{
				if( defined( $genders{ $_ } ) )
				{
					print $genders{ $_ };
					return $genders{ $_ };
				}
			}
		}
	}
	print 'unknown';
	return "unknown";
}

#change the NCBI G1K solid data into fastq for MAQ
sub solid2Fastq
{
	croak "Usage: solid2Fastq in.fastq out.fast" unless @_ == 2;
	my $input = shift;
	my $output = shift;
	
	if( $input =~ /\.gz$/ )
	{
		open( READS, "gunzip -c $input |" ) or die "Cannot open gzipped fastq file\n";
		open( OUT, "|gzip -c > $output" ) or die "Cannot create otuput file\n";
	}
	else
	{
		open( READS, $input ) or die "Failed to open reads file";
		open( OUT, ">$output" ) or die "Cannot create otuput file\n";
	}
	
	while( <READS> )
	{
		my $name = $_;
		my $seq = <READS>;
		$seq = substr( $seq, 2 );
		my $qname = <READS>;
		my $quals = <READS>;
		$quals = substr( $quals, 2 );
		
		#translate the bases from nums to nucleotides
		$seq =~ tr/0123./ACGTN/;
		
		print OUT "$name$seq$qname$quals";
	}
	close( READS );
	close( OUT );
}

1;
