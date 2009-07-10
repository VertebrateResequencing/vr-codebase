package QC;
use strict;
use Carp;
use File::Basename;
use File::Spec;
use Cwd;

use Mapping_SLX_Maq;
use Mapping_454_ssaha;
#use QC_SLX_Maq;
use QC_454_ssaha;

my $G1K = $ENV{ 'G1K' };
my $SAM_TOOLS = $G1K.'/bin/samtools';
my $R_LOCATION = '/software/R-2.7.1/bin/R';
my $IMAGE_MAGIC = '/usr/bin/convert';

=pod

=head1 NAME

QC

=head1 SYNOPSIS

A module that is that is used to for carrying out QC checks on sequencing data

=head1 REQUIRES

Perl5.8.8

=head1 DESCRIPTION

A module for carrying out QC on the Vertebrate Resequencing Informatics 
data hierarchy. The hierarchy can consist of multiple different sequencing 
technologies.

=head1 METHODS

=head2 qc

	Arg [1]    : root of the mapping hierarchy
	Arg [2]    : LSF queue for jobs
	Arg [3]    : file of lanes to be QC'd
	Example    : qc( '$G1K/MOUSE/MAPPING', 'normal', 'lanes.fofn');
	Description: QC the lanes (if not done already) specified in the lanes.fofn file
	Returntype : none

=cut

sub qc
{
	croak "Usage: rootDir data_root_dir mapping_root_dir lsf_queue lanes_fofn" unless @_ == 4;
	
	my $dRoot = shift;
	my $mRoot = shift;
	my $lsf_queue = shift;
	my $lanes_fofn = shift;
	
	croak "Cant find the lanes_fofn file: $lanes_fofn\n" unless -f $lanes_fofn;
	
	croak "Cant find the data root directory: $dRoot\n" unless -d $dRoot;
	croak "Cant find the mapping root directory: $mRoot\n" unless -d $mRoot;
	
	my %lanes;
	
	open( LANES, $lanes_fofn ) or die "Cant open lanes fofn file: $!\n";
	while( <LANES> )
	{
		chomp;
		
		my $dataDirLaneAbsPath = $_;
		my $mapDirLaneAbsPAth = $_;
		if( $_ !~ /^$dRoot.*/ )
		{
			$dataDirLaneAbsPath = $dRoot.'/'.$_;
			$mapDirLaneAbsPAth = $mRoot.'/'.$_;
		}
		
		if( ! -d $dataDirLaneAbsPath )
		{
			print "Cant find lane data directory: $dataDirLaneAbsPath\n";
			next;
		}
		
		if( ! -d $mapDirLaneAbsPAth )
		{
			print "Cant find lane mapping directory: $mapDirLaneAbsPAth\n";
			next;
		}
		
		my ($project, $sample, $platform, $lib, $lane ) = split "/", $_;
		
		#make the qc subdirectory
		if( ! -d $dataDirLaneAbsPath.'/qc' )
		{
			mkdir( $dataDirLaneAbsPath.'/qc' ) or die "Cannot create QC subdirectory in lane: $dataDirLaneAbsPath\n";
		}
		
		chdir( $dataDirLaneAbsPath.'/qc' );
		
		#create the GC plot
		if( ! -f "gc.gif" )
		{
			#get the fastq names
			my @fastq = HierarchyUtilities::getFastqNamesLane( $dataDirLaneAbsPath );
			
			#make the plot
			if( length( $fastq[ 2 ] ) > 0 )
			{
				my $cmd = 'bsub -q normal -o qc.o -e qc.e perl -w -e ';
				gcContentHistogram( 'gc.histogram.1', $fastq[ 1 ], 'gc.gif', 'gc.histogram.2', $fastq[ 2 ] );
			}
			else
			{
				
			}
		}
		
		if( $platform eq "SLX" )
		{
			if( Mapping_SLX_Maq::isLaneMapped( $mapDirLaneAbsPAth ) && ! Mapping_SLX_Maq::mappingInProgress( $mapDirLaneAbsPAth ) )
			{
				print "QC'ing Lane: $mapDirLaneAbsPAth\n";
				
				if( ! -f 'inserts.gif' )
				{
					#create the insert size plot
				}
				#QC_SLX_Maq::qcLane( $laneAbsPath, $lsf_queue );
				chromosomeDistribution( $lsf_queue, $dataDirLaneAbsPath, $mapDirLaneAbsPAth );
			}
		}
		elsif( $platform eq "454" )
		{
			if( Mapping_454_ssaha::isLaneMapped( $mapDirLaneAbsPAth ) && ! Mapping_454_ssaha::mappingInProgress( $mapDirLaneAbsPAth ) )
			{
				print "QC'ing Lane: $mapDirLaneAbsPAth\n";
				#QC_454_ssaha::qcLane( $laneAbsPath, $lsf_queue );
				chromosomeDistribution( $lsf_queue, $dataDirLaneAbsPath, $mapDirLaneAbsPAth );
			}
		}
	}
	close( LANES );
}

=head2 chromosomeDistribution

	Arg [1]    : lane directory
	Arg [2]    : LSF queue for jobs
	Example    : qc( '$G1K/MOUSE/MAPPING', 'normal', 'lanes.fofn');
	Description: Create a distribution of the reads over the chrs
	Returntype : none

=cut

sub chromosomeDistribution
{
	croak "Usage: chromosomeDistribution lsf_queue data_lane_dir mapping_lane_dir" unless @_ == 3;
	
	my $lsf_queue = shift;
	my $data_lane_dir = shift;
	my $map_lane_dir = shift;
	
	croak "Cant find lane directory: $data_lane_dir\n" unless -d $data_lane_dir;
	croak "Cant find lane directory: $map_lane_dir\n" unless -d $map_lane_dir;
	
	chdir( $data_lane_dir );
	
	croak "Cant rmdup.bam file in: $map_lane_dir\n" unless -f $map_lane_dir.'/rmdup.bam';
	
	mkdir( 'qc' ) unless -d 'qc';
	
	chdir( 'qc' );
	
	my %numReadsPerChr;
	open( BAM, "$SAM_TOOLS view $map_lane_dir/rmdup.bam |" ) or die "Cannot open bam file: $!\n";
	while( <BAM> )
	{
		chomp;
		my @s = split( /\s+/, $_ );
		if( defined( $numReadsPerChr{ $s[ 2 ] } ) )
		{
			$numReadsPerChr{ $s[ 2 ] } ++;
		}
		else
		{
			$numReadsPerChr{ $s[ 2 ] } = 1;
		}
	}
	close( BAM );
	
	my $gender = G1KUtilities::path2Gender();
	
	my %fai;
	if( $gender eq 'male' || $gender eq 'unknown' )
	{
		%fai = %{ readFAI( $Mapping_SLX_Maq::HUMAN_MALE_REF_FAI ) };
	}
	else
	{
		%fai = %{ readFAI( $Mapping_SLX_Maq::HUMAN_FEMALE_REF_FAI ) };
	}
	
	#normalise the coverage per kb and make a plot and table
	my %normalised;
	open( DIST, ">reads_per_chr_kb.tab" ) or die "cant create file in qc subdirectory: $data_lane_dir\n";
	foreach( sort( keys( %numReadsPerChr ) ) )
	{
		print DIST "$_\t";
	}
	print DIST "\n";
	
	foreach( sort( keys( %numReadsPerChr ) ) )
	{
		$normalised{ $_ } = $numReadsPerChr{ $_ } / ( int( $fai{ $_ } / 1000 ) );
		print DIST $normalised{ $_ }."\t";
	}
	print DIST "\n";
	close( DIST );
	
	#call R and make a plot....
}


sub gcContentHistogram
{
	croak "Usage: gcContentHistogram gc_hist_output1 fastq1 gc_gif_output [gc_hist_output2 fastq2]" unless @_ == 3 || @_ == 5;
	
	my $gc_hist1 = shift;
	my $fastq1 = shift;
	my $gc_gif = shift;
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
	print R 'bitmap(\'g.bmp\');';
	print R 'plot(gc_1[,1],gc_1[,2],xlab="%GC",ylab="Frequency",type="l");';
	
	if( length( $fastq2  ) > 0 )
	{
		print R 'gc_2<-read.table(\''.$gc_hist2.'\');';
		print R 'lines(gc_2[,1],gc_2[,2],type="l");';

		print "$fastq2 Average GC: $avgContent2\n";
	}
	print R 'dev.off();';
	close( R );
	
	system( $R_LOCATION.' CMD BATCH gc.R;'.$IMAGE_MAGIC.' g.bmp '.$gc_gif.';rm g.bmp' );
	print "$fastq1 Average GC: $avgContent1\n";
}

sub createInsertGraph
{
	croak "Usage: createInsertGraph sam_file maxInsert output_file_gif" unless @_ == 3;
	my $sam_file = shift;
	my $maxInsert = shift;
	my $output_file = shift;
	
	croak "Cant find sam file: $sam_file" unless -f $sam_file;
	croak "Usage: createInsertGraph sam_file maxInsert output_file\n" unless $maxInsert =~ /^\d+$/;
	
	if( $sam_file =~ /\.gz$/ )
	{
		open( MV, "gunzip -c $sam_file |" ) or die "Cannot open file: $sam_file\n";
	}
	else
	{
		open( MV, "$sam_file" ) or die "Cannot open file: $sam_file\n";
	}
	
	my %inserts;
	$inserts{ 0 } = 0;
	while( <MV> )
	{
		chomp;
		my @s = split( /\t/, $_ );
		
		my $insert = abs( $s[ 8 ] );
		if( $insert != 0 && $insert < $maxInsert )
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
	close( MV );
	
	open( I, ">$sam_file.inserts" )  or die "cannot create inserts file\n";
	foreach( sort {$a <=> $b}( keys( %inserts ) ) )
	{
		print I $_." ".$inserts{ $_ }."\n";
	}
	close( I );
	
	open( R, ">$sam_file.insert.R" ) or die "Cannot create R script\n";
	print R "inserts<-read.table('\$mapview.inserts');\n"; #*** used to say '$mapview' but that var doesn't exist
	print R "bitmap('t.bmp');\n";
	print R "plot(inserts[,1],inserts[,2],main='Insert Sizes',xlab='Insert Size', ylab='Frequency', type='p', axes='F' );\n";
	print R "at.labels=axis(2);\n";
	print R "at.labels=axis(1)\n";
	print R "axis(1,at=seq(at.labels[1], at.labels[length(at.labels)], (at.labels[2]-at.labels[1])/4), labels=F )\n";
	print R "dev.off();\n";
	close( R );
	
	system( $R_LOCATION." CMD BATCH $sam_file.insert.R;$IMAGE_MAGIC t.bmp $output_file;rm t.bmp" );
}

sub readFAI
{
	croak "Usage: readFAI FAI" unless @_ == 1;
	
	my $fai = shift;
	
	my %fai;
	open( FAI, "$fai" ) or die "Cannot open fai file";
	while( <FAI> )
	{
		chomp;
		my @s = split( /\s+/, $_ );
		$fai{ $s[ 0 ] } = $s[ 1 ];
	}
	close( FAI );
	
	return \%fai;
}
