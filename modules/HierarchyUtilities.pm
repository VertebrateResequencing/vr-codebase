package HierarchyUtilities;
use strict;
use Carp;
use File::Basename;
use File::Spec; 
use File::Copy;

=pod

=head1 NAME

=head1 SYNOPSIS

function to import the DCC tab index files into the hierarchy

=head1 ARGUMENTS

=head1 DESCRPTION

=head1 AUTHOR

Thomas Keane, I<tk2@sanger.ac.uk>

=cut

my $DFIND = '/software/solexa/bin/dfind';
my $MPSA_DOWNLOAD = '/software/solexa/bin/mpsa_download';
my $FASTQ_CHECK = '/software/solexa/bin/fastqcheck';
#my $LSF_QUEUE = 'normal';
my $LSF_QUEUE = 'small';
my $CLIP_POINT_SCRIPT='sh ~tk2/code/vert_reseq/user/tk2/mapping/slx/findClipPoints.sh';

sub importExternalData
{
	croak "Usage: importExternalData populations/samples_csv DCC_index hierarchy_parent_directory mapping_hierarchy_directory fastq_directory" unless @_ == 5;
	
	my $pop_csv = shift;
	my $index_file = shift;
	my $data_hierarchy_dir = shift;
	my $mapping_hierarchy_dir = shift;
	my $fastqDir = shift;
	
	croak "Can't find populations_csv file" unless -f $pop_csv;
	croak "Can't find index file" unless -f $index_file;
	croak "Can't find fastq directory" unless -d $fastqDir;
	croak "Can't find data hierarchy parent directory" unless -d $data_hierarchy_dir;
	croak "Can't find mapping hierarchy parent directory" unless -d $mapping_hierarchy_dir;
	
	my %study_vs_id = (
			28911   => 'LowCov',  # pilot 1
		    28919   => 'Trio',    # pilot 2
		    28917   => 'Exon',    # pilot 3
		    'SRP000031'   => 'LowCov',    # pilot 1
		    'SRP000032'   => 'Trio',    # pilot 2
		    'SRP000033'   => 'Exon',    # pilot 3
		    'Pilot1'   => 'LowCov',    # pilot 1
		    'Pilot2'   => 'LowCov',    # pilot 1
		    'Pilot3'   => 'LowCov',    # pilot 1
		    );
	
	#read in the project spreadsheet with individuals
	open( PCSV, $pop_csv ) or die "Failed to open population csv sheet: $!\n";
	my %populations;
	while( <PCSV> )
	{
		chomp;
		next unless length( $_ ) > 0 && $_ !~ /^\s+$/ && $_ !~ /^,+$/;
		
		my @s = split( /,/, $_ );
		$populations{ $s[ 1 ] } = $s[ 5 ];
	}
	close( PCSV );
	
	my %exp_unpaired;
	my %exp_read1;
	my %exp_read2;
	my %data_exp_path; #directory within data hierarchy
	my %mapping_exp_path; #mapping directory
	
	my %srrCount;
	
	open( LCSV, $index_file ) or die "Failed to open lane index file: $!\n";
	while( <LCSV> )
	{
		chomp;
		my $original = $_;
		
		next if( $_ =~ /^\s+$/ || length( $_ ) == 0 || $_ =~ /^FASTQ_FILE/ );
		
		my @info = split( /\t/, $_ );
		croak "Can't find study for lane ID: $info[ 3 ]\n" if( ! defined $study_vs_id{ $info[ 3 ] } );
		croak "Can't find population for lane ID: $info[ 9 ]\n" if( ! defined $populations{ $info[ 9 ] } );
		
		#check SRR is unique
		if( defined( $srrCount{ $info[ 2 ] } ) )
		{
			if( ( $srrCount{ $info[ 2 ] } == 2 && $info[ 12 ] =~ /.*[solexa|illumina].*/ ) || ( $srrCount{ $info[ 2 ] } == 3 && $info[ 12 ] =~ /.*[roche|454].*/ ) || $info[ 19 ] eq 'SINGLE' )
			{
				print "ERROR: Non-unique SRR found: $info[ 2 ]\n";
				exit;
			}
			else
			{
				$srrCount{ $info[ 2 ] } ++;
			}
		}
		else
		{
			$srrCount{ $info[ 2 ] } = 1;
		}
		
		#unpaired
		if( $info[ 18 ] eq 'SINGLE' || ( $info[ 18 ] eq 'PAIRED' && @info < 20 ) )
		{
			print "UP: $info[ 0 ] $info[ 18 ]\n";
			$exp_unpaired{ $info[ 2 ] } = basename( $info[ 0 ] );
		}
		elsif( $info[ 18 ] eq 'PAIRED' )
		{
			print "P: $info[ 0 ] $info[ 18 ]\n";
			if( ! defined( $exp_read1{ $info[ 2 ] } ) )
			{
				$exp_read1{ $info[ 2 ] } = basename( $info[ 0 ] );
			}
			else
			{
				$exp_read2{ $info[ 2 ] } = basename( $info[ 0 ] );
			}
		}
		
		my $dpath = $data_hierarchy_dir.'/'.$study_vs_id{ $info[ 3 ] }.'-'.$populations{ $info[ 9 ] }.'/'.$info[ 9 ];
		my $mpath = $mapping_hierarchy_dir.'/'.$study_vs_id{ $info[ 3 ] }.'-'.$populations{ $info[ 9 ] }.'/'.$info[ 9 ];
		
		#technology
		if( $info[ 12 ] =~ /454|roche/i )
		{
			$dpath = $dpath.'/454';
			$mpath = $mpath.'/454';
		}
		elsif( $info[ 12 ] =~ /solid/i )
		{
			$dpath = $dpath.'/SOLID';
			$mpath = $mpath.'/SOLID';
		}
		elsif( $info[ 12 ] =~ /solexa|illumina/i )
		{
			$dpath = $dpath.'/SLX';
			$mpath = $mpath.'/SLX';
		}
		else
		{
			croak "Cannot determine sequencing technology: $info[ 12 ] in\n$original\n";
		}
		$dpath .= '/'.$info[ 14 ].'/'.$info[ 2 ];
		$mpath .= '/'.$info[ 14 ].'/'.$info[ 2 ];
		
		#check that the data has not already been imported
		if( -d $dpath )
		{
			print "Already Imported: $info[ 0 ]\n";
			next;
		}
		
		$data_exp_path{ $info[ 2 ] } = $dpath;
		$mapping_exp_path{ $info[ 2 ] } = $mpath;
	}
	close( LCSV );
	
	foreach( keys( %data_exp_path ) )
	{
		#check the fastq's exist
		if( defined( $exp_read1{ $_ } ) && ! defined( $exp_read2{ $_ } ) )
		{
			#turn it into an unpaired entry
			$exp_unpaired{ $_ } = $exp_read1{ $_ };
			delete( $exp_read1{ $_ } );
			print "No read2 entry so making read unpaired: $exp_unpaired{ $_ }\n";
			
			#print "ERROR: Failed to find read2 for: $_\n";
			#next;
		}
		
		if( defined( $exp_read1{ $_ } ) && ! -f $fastqDir.'/'.$exp_read1{ $_ } )
		{
			print "WARNING: Failed on $_: Failed to read file $exp_read1{ $_ }\n";
			next;
		}
		elsif( defined( $exp_read2{ $_ } ) && ! -f $fastqDir.'/'.$exp_read2{ $_ } )
		{
			print "WARNING: Failed on $_: Failed to read file $exp_read2{ $_ }\n";
			next;
		}
		elsif( defined( $exp_unpaired{ $_ } ) && ! -f $fastqDir.'/'.$exp_unpaired{ $_ } )
		{
			print "WARNING: Failed on $_: Failed to read file $exp_unpaired{ $_ }\n";
			next;
		}
		
		my $pdirName = $data_exp_path{ $_ };
		my $mdirName = $mapping_exp_path{ $_ };
		
		#make the directory in the hierarchy
		if( ! -d $pdirName )
		{
			print "Making path: $pdirName\n";
			system( "mkdir -p ".$pdirName );
			if( ! -d $pdirName )
			{
				croak "ERROR: Failed to create directory: ".$pdirName."\n";
			}
			
			print "Making path: $mdirName\n";
			system( "mkdir -p ".$mdirName );
			if( ! -d $mdirName )
			{
				croak "ERROR: Failed to create directory: ".$mdirName."\n";
			}
		}
		
		#write a meta file with the read file names
		open( META, ">$mdirName/meta.info" ) or die "Can't open meta.info file in $mdirName\n";
		open( META1, ">$pdirName/meta.info" ) or die "Can't open meta.info file in $pdirName\n";
		
		#copy the files into the directory
		if( defined( $exp_unpaired{ $_ } ) && ! -l $pdirName.'/'.$exp_unpaired{ $_ } && ! -f $pdirName.'/'.$exp_unpaired{ $_ } )
		{
			copy( $fastqDir.'/'.$exp_unpaired{ $_ }, $pdirName.'/'.$exp_unpaired{ $_ } ) or die "Failed to copy file: ".$exp_unpaired{ $_ };
			symlink( $pdirName.'/'.$exp_unpaired{ $_ }, $mdirName.'/'.$exp_unpaired{ $_ } );
			system( 'bsub -q small -o fastqcheck.o -e fastqcheck.e "zcat '.$pdirName.'/'.$exp_unpaired{ $_ }.' | awk \'{print \$1}\' | /software/solexa/bin/fastqcheck > '.$pdirName.'/'.$exp_unpaired{ $_ }.'.fastqcheck; ln -fs '.$pdirName.'/'.$exp_unpaired{ $_ }.'.fastqcheck '.$mdirName.'/'.$exp_unpaired{ $_ }.'.fastqcheck"' );
			print META "read0:".$exp_unpaired{ $_ }."\n";
			print META1 "read0:".$exp_unpaired{ $_ }."\n";
		}
		
		if( defined( $exp_read1{ $_ } ) && ! -l $pdirName.'/'.$exp_read1{ $_ } && ! -f $pdirName.'/'.$exp_read1{ $_ } )
		{
			copy( $fastqDir.'/'.$exp_read1{ $_ }, $pdirName.'/'.$exp_read1{ $_ } ) or die "Failed to copy file: ".$exp_read1{ $_ };
			symlink( $pdirName.'/'.$exp_read1{ $_ }, $mdirName.'/'.$exp_read1{ $_ } ) or die "Failed to sym link read into mapping directory: ".$mdirName.'/'.$exp_read1{ $_ };
			system( 'bsub -q small -o fastqcheck.o -e fastqcheck.e "zcat '.$pdirName.'/'.$exp_read1{ $_ }.' | awk \'{print \$1}\' | /software/solexa/bin/fastqcheck > '.$pdirName.'/'.$exp_read1{ $_ }.'.fastqcheck; ln -fs '.$pdirName.'/'.$exp_read1{ $_ }.'.fastqcheck '.$mdirName.'/'.$exp_read1{ $_ }.'.fastqcheck"' );
			print META "read1:".$exp_read1{ $_ }."\n";
			print META1 "read1:".$exp_read1{ $_ }."\n";
		}
		
		if( defined( $exp_read2{ $_ } ) && ! -l $pdirName.'/'.$exp_read2{ $_ } && ! -f $pdirName.'/'.$exp_read2{ $_ } )
		{
			copy( $fastqDir.'/'.$exp_read2{ $_ }, $pdirName.'/'.$exp_read2{ $_ } ) or die "Failed to copy file: ".$exp_read2{ $_ };
			symlink( $pdirName.'/'.$exp_read2{ $_ }, $mdirName.'/'.$exp_read2{ $_ } ) or die "Failed to sym link read into mapping directory: ".$mdirName.'/'.$exp_read2{ $_ };
			system( 'bsub -q small -o fastqcheck.o -e fastqcheck.e "zcat '.$pdirName.'/'.$exp_read2{ $_ }.' | awk \'{print \$1}\' | /software/solexa/bin/fastqcheck > '.$pdirName.'/'.$exp_read2{ $_ }.'.fastqcheck; ln -fs '.$pdirName.'/'.$exp_read2{ $_ }.'.fastqcheck '.$mdirName.'/'.$exp_read2{ $_ }.'.fastqcheck"' );
			print META "read2:".$exp_read2{ $_ }."\n";
			print META1 "read2:".$exp_read2{ $_ }."\n";
		}
		
		close( META );
		close( META1 );
	}
}


=head2 importInternalData

  Arg [1]    : file of NPG project names, one per line
  Arg [2]    : Data directory to build hierarchy in
  Arg [3]    : Analysis directory to build hierarchy in
  Example    : importInternalData( 'mouse.proj', '$G1K/MOUSE/DATA', '$G1K/MOUSE/MAPPING');
  Description: Imports all new sequence from NPG into a hierarchy.  Existing fastq are skipped.
  Returntype : none

=cut

sub importInternalData {
    croak "Usage: importInternalData project_list_file data_hierarchy_parent_directory analysis_hierarchy_parent_directory" unless @_ == 3;
    my $projectFile = shift;
    my $dhierarchyDir = shift;
    my $ahierarchyDir = shift;

    croak "Can't find data hierarchy directory\n" unless -d $dhierarchyDir;
    croak "Can't find analysis hierarchy directory\n" unless -d $ahierarchyDir;
    croak "Can't find projects file\n" unless -f $projectFile;
    
    my %projects;
    open( my $PROJS, "$projectFile" ) or die "Cannot open projects file\n";
    while(my $project =  <$PROJS> ) {
	chomp $project;
	$projects{$project} = {};

	open(my $LIBS, q[-|], qq[$DFIND -project "$project" -libraries] ) or die "Cannot run dfind on project: $project\n";
	while(my $lib = <$LIBS> ) {
	    chomp $lib;
	    $projects{$project}{$lib} = [];

	    open my $LANES, q[-|], qq[$DFIND -project "$project" -library "$lib" -filetype fastq] or die "Can't get lanes for $project $lib: $!\n";
	    while(my $lane = <$LANES> ) {
		chomp $lane;
		push @{$projects{$project}{$lib}}, $lane;
	    }
	    close( $LANES );
	}
	close( $LIBS );
    }
    close( $PROJS )
    &buildInternalHierarchy(\%projects,$dhierarchyDir,$ahierarchyDir);
}



sub buildInternalHierarchy {
    my $projecthash = shift; # ref to hash of projectnames->samplenames->lists of fastq
    my $dhierarchyDir = shift;
    my $ahierarchyDir = shift;

    croak "Can't find data hierarchy directory\n" unless -d $dhierarchyDir;
    croak "Can't find analysis hierarchy directory\n" unless -d $ahierarchyDir;

    foreach my $project (keys %$projecthash){
	print "Updating project: $project\n";
	
	my $project1 = $project;
	$project1 =~ s/\W+/_/g;
	
	my $projPath = $dhierarchyDir.'/'.$project1;
	if( ! -d $projPath )
	{
		mkdir $projPath or die "Cannot create directory $projPath\n";
	}
	
	my $aprojPath = $ahierarchyDir.'/'.$project1;
	if( ! -d $aprojPath )
	{
		mkdir $aprojPath or die "Cannot create directory $aprojPath\n";
	}
	
	my $numLibraries = 0;
	
	foreach my $sample (keys %{$projecthash->{$project}}){
	    chomp;
	    print "Updating library: $sample\n";
	    
	    #hack for G1K where sample starts with the individual
	    my $individual = 1;
	    if( $sample =~ /^NA\d+/ )
	    {
		    $individual = (split( /-/, $sample ) )[ 0 ];
	    }
	    
	    my $sample1 = $sample;
	    $sample1 =~ s/\W+/_/g;
	    
	    my $path = $projPath.'/'.$individual;
	    if( ! -d $path )
	    {
		    mkdir $path or die "Cannot create directory $path\n";
	    }
	    $path .= '/SLX';
	    if( ! -d $path )
	    {
		    mkdir $path or die "Cannot create directory $path\n";
	    }
	    $path .= '/'.$sample1;
	    if( ! -d $path )
	    {
		    mkdir $path or die "Cannot create directory $path\n";
	    }
	    
	    my $apath = $aprojPath.'/'.$individual;
	    if( ! -d $apath )
	    {
		    mkdir $apath or die "Cannot create directory $apath\n";
	    }
	    $apath .= '/SLX';
	    if( ! -d $apath )
	    {
		    mkdir $apath or die "Cannot create directory $apath\n";
	    }
	    $apath .= '/'.$sample1;
	    if( ! -d $apath )
	    {
		    mkdir $apath or die "Cannot create directory $apath\n";
	    }
	    
	    my @dirs;
	    foreach my $fastq (@{$projecthash->{$project}{$sample}}){
		chomp;
		my $fastq = basename( $fastq );
		
		my @s = split( /\./, $fastq );
		my @s1 = split( '_', $s[ 0 ] );
		
		my $lPath = '';
		my $alPath = '';
		if( $s1[ 1 ] eq 's' ) #hack for old read file names with s character
		{
			$lPath = $path.'/'.$s1[ 0 ].'_'.$s1[ 2 ];
			$alPath = $apath.'/'.$s1[ 0 ].'_'.$s1[ 2 ];
		}
		else
		{
			$lPath = $path.'/'.$s1[ 0 ].'_'.$s1[ 1 ];
			$alPath = $apath.'/'.$s1[ 0 ].'_'.$s1[ 1 ];
		}
		
		push( @dirs, $lPath );
		
		if( ! -d $lPath )
		{
			mkdir $lPath or die "Cannot create directory $lPath\n";
		}
		
		if( ! -d $alPath )
		{
			mkdir $alPath or die "Cannot create directory $alPath\n";
		}
		
		#check if the fastq file already exists
		my $fqPath = $lPath.'/'.$fastq;
		
		#check if the gzipped reads are already on disk
		if( $s1[ 1 ] eq 's' ) #hack for old read file names with s character
		{
			my $fastq1 = $s1[ 0 ].'_'.$s1[ 2 ].'_1.fastq.gz'; #read1 fastq
			
			if(! -s $fqPath && ! -s $lPath.'/'.$fastq1 )
			{
				chdir( $lPath );
				
				print "Requesting file $fastq from MPSA.....\n";
				my $random = int(rand( 100000000 ));
				
				#use mpsadownload to get fastq
				# jws 2009-02-11 - do this inline, and then do any parallelisation outside this module
				#my $cmd = qq[bsub -J mpsa.$random -o import.o -e import.e -q $LSF_QUEUE "$MPSA_DOWNLOAD -c -f $fastq > $fastq"];
				my $cmd = "$MPSA_DOWNLOAD -c -f $fastq > $fastq 2> import.e";
				system( $cmd );
				unless (-s $fastq){
				    print "Error retrieving $fastq with mpsa_download: $!\n";
				    next;
				}
				
				$fastq =~ /(\d+)_s_(\d+)\.fastq/;
				
				my $fastq1 = $1.'_'.$2.'_1.fastq';
				my $fastq2 = $1.'_'.$2.'_2.fastq';
				
				#write the meta.info file
				open( META, ">$lPath/meta.info" ) or die "Cannot create meta.info file\n";
				print META "read1:$fastq1.gz\n";
				print META "read2:$fastq2.gz\n";
				close( META );
				
				print "Writing meta info for: $fastq1.gz\n";
				print "Writing meta info for: $fastq2.gz\n";
				
				#link the meta file into mapping hierarchy
				system( "ln -fs $lPath/meta.info $alPath/meta.info" );

				#create a bsub job to split the fastq
				$cmd = qq[bsub -J split.$random -o import.o -e import.e -q $LSF_QUEUE perl -w -e "use AssemblyTools;AssemblyTools::sanger2SplitFastq( '$fastq', '$fastq1', '$fastq2');unlink '$fastq';"];
				system( $cmd );
				
				#work out the clip points (if any)
				$cmd = qq[bsub -J clip.$random -w "done(split.$random)" -o import.o -e import.e -q $LSF_QUEUE "$CLIP_POINT_SCRIPT $fastq1 meta.info 1;$CLIP_POINT_SCRIPT $fastq2 meta.info 2"];
				system( $cmd );	

				#run fastqcheck
				$cmd = qq[bsub -J fastqcheck.$random.1 -o import.o -e import.e -q $LSF_QUEUE -w "done(split.$random)" "cat $fastq1 | $FASTQ_CHECK > $fastq1.gz.fastqcheck;ln -fs $lPath/$fastq1.gz.fastqcheck $alPath/$fastq1.gz.fastqcheck"];
				system( $cmd );
				
				$cmd = qq[bsub -J fastqcheck.$random.2 -o import.o -e import.e -q $LSF_QUEUE -w "done(split.$random)" "cat $fastq2 | $FASTQ_CHECK> $fastq2.gz.fastqcheck;ln -fs $lPath/$fastq2.gz.fastqcheck $alPath/$fastq2.gz.fastqcheck"];
				system( $cmd );
				
				#gzip the split fastq files
				$cmd = qq[bsub -q $LSF_QUEUE -w "done(clip.$random)&&done(fastqcheck.$random.*)" -o import.o -e import.e "gzip $fastq1 $fastq2; ln -fs $lPath/$fastq1.gz $alPath/$fastq1.gz;ln -fs $lPath/$fastq2.gz $alPath/$fastq2.gz"];
				system( $cmd );

			}
			else
			{
				print "Lane already in hierarchy: $fastq\n";
			}
		}
		elsif( ! -s $fqPath.'.gz' && ! -s $fqPath )
		{
			print "Requesting file $fastq from MPSA.....\n";
			
			chdir( $lPath );
			my $random = int(rand( 100000000 ));
			
			#use mpsadownload to get fastq.
			# jws 2009-02-11 - do this inline, and then do any parallelisation outside this module
			#my $cmd = qq[bsub -J mpsa.$random -o import.o -e import.e -q $LSF_QUEUE "$MPSA_DOWNLOAD -c -f $fastq > $fastq"];
			my $cmd = "$MPSA_DOWNLOAD -c -f $fastq > $fastq 2> import.e";
			system( $cmd );
			unless ( -s $fastq){	# TODO: add md5 check if mpsa ever add it
			    print "Error retrieving $fastq with mpsa_download: $!\n";
			    next;
			}
			
			#run fastqcheck
			$cmd = qq[bsub -J fastqcheck.$random -o import.o -e import.e -q $LSF_QUEUE "cat $fastq | $FASTQ_CHECK > $fastq.gz.fastqcheck; ln -fs $lPath/$fastq.gz.fastqcheck $alPath/$fastq.gz.fastqcheck"];
			system( $cmd );
			
			print "Writing meta.info for: $fastq.gz\n";
			
			#then write the meta info file
			open( META, ">$lPath/meta.info" ) or die "Cannot create meta file in $lPath\n";
			
			#figure out if its paired or unpaired read
			if( $fastq =~ /^\d+_\d+\.fastq$/ ) #unpaired
			{
				print META "read0:$fastq.gz\n";
			}
			elsif( $fastq =~/^\d+_\d+_1.fastq$/ )
			{
				print META "read1:$fastq.gz\n";
			}
			elsif( $fastq =~/^\d+_\d+_2.fastq$/ )
			{
				print META "read2:$fastq.gz\n";
			}
			else
			{
				print "ERROR: Can't determine read information: $fastq\n";
				exit;
			}
			close( META );
			
			#link the meta file into mapping hierarchy
			system( "ln -fs $lPath/meta.info $alPath/meta.info" );
			
			#work out the clip point (if any)
			my $readNum = 1;
			if( $fastq =~ /.*_2\.fastq/ )
			{
				$readNum = 2;
			}
			$cmd = qq[bsub -J clip.$random -o import.o -e import.e -q $LSF_QUEUE "$CLIP_POINT_SCRIPT $fastq $lPath/meta.info $readNum"];
			system( $cmd );
			
			#gzip the fastq file
			
			$cmd = qq[bsub -w "done(fastqcheck.$random)&&done(clip.$random)" -q $LSF_QUEUE -o import.o -e import.e "gzip $fastq; ln -fs $lPath/$fastq.gz $alPath/$fastq.gz"];
			system( $cmd );

		}
		else
		{
			print "Fastq already in hierarchy: $fastq\n";
		}

	    }
	    $numLibraries ++;
	}
	
	print "WARNING: No libraries found for $project\n" unless $numLibraries > 0;
    }

}



=head2 checkInternalFastq

  Arg [1]    : DATA lane directory to check
  Example    : my $is_ok = checkInternalFastq('/lustre/sf4/1kgenomes/G1K/SANGER/1000Genomes_B1_TOS/NA20534/SLX/NA20534_TOS_1/1513_2');
  Description: compares NPG/MPSA fastqcheck to fastqcheck(s) for retrieved fastq.  Handles split fastq as well as unsplit.  Note, uses dfind to retrieve the NPG fastqcheck file(s).
  Returntype : 1 if fastqcheck counts are the same.

=cut

sub checkInternalFastq {
    croak "Usage: checkInternalFastql hierarchy_fastq_dir" unless @_ == 1;
    my $fastqdir = shift;
    $fastqdir =~ s|/$||;    # remove any trailing separator
    my ($run,$lane) = $fastqdir=~m|.*/(\d+)_(\d+)$|;
    my $check_bad;
    unless ($run && $lane){
	carp "Can't find run and lane from $fastqdir";
	return undef;
    }

    open my $FQCHECK, q[-|], qq[$DFIND -run $run -lane $lane -filetype fastqcheck];
    while (<$FQCHECK>){
	chomp;
	my ($orig_seqs,$orig_len) = getSeqAndLengthFromFastqcheck($_);
	my ($lanename) = m|.*/(\w+)\.fastqcheck$|;
	my ($new_seqs,$new_len);
	if ($lanename eq "${run}_s_${lane}"){
	    # old fastq that needed splitting.
	    # should have two corresponding fastqchecks from the split files
	    foreach my $ori (1,2){
		my $fq = "${run}_${lane}_${ori}";
		my ($seqs, $len) = getSeqAndLengthFromFastqcheck("$fastqdir/$fq.fastq.gz.fastqcheck");
		$new_seqs += $seqs;
		$new_len  += $len;
	    }
	    $new_seqs = $new_seqs/2;	# should be same number of seqs in both
					# split files
	}
	elsif ($lanename =~ /^${run}_${lane}(_[12])?$/){
	    # files are pre-split, so just need corresponding fastqcheck
	    ($new_seqs, $new_len) = getSeqAndLengthFromFastqcheck("$fastqdir/$lanename.fastq.gz.fastqcheck");
	}
	else {
	    warn "$_ unrecognised fastqcheck name\n";
	    next;
	}
	unless ($orig_seqs == $new_seqs && $orig_len == $new_len){
	    $check_bad = 1;
	}
    }
    my $check_ok = $check_bad ? 0 : 1;
    return $check_ok;
}


sub getSeqAndLengthFromFastqcheck {
    my $fastqcheck = shift;
    open (my $FQC, $fastqcheck) or croak "Can't open $fastqcheck: $!\n";
    my $header = <$FQC>;
    close $FQC;
    my ($seqcount,undef,$length) = split /\s+/,$header;
    return ($seqcount,$length);
}

=pod
	Function to delete lanes from a hierarchy.
	Input is a meta index file
=cut
sub deleteIndexLanes
{
	croak "Usage: deleteIndexLanes samples_csv index_file hierarchy_root_dir" unless @_ == 3;
	
	my $pop_csv = shift;
	my $index_file = shift;
	my $hierarchy_dir = shift;
	
	croak "Can't find populations_csv file" unless -f $pop_csv;
	croak "Can't find index file" unless -f $index_file;
	croak "Can't find hierarchy parent directory" unless -d $hierarchy_dir;
	
	my %study_vs_id = (28911   => 'LowCov',  # pilot 1
		    28919   => 'Trio',    # pilot 2
		    28917   => 'Exon',    # pilot 3
		    'SRP000031'   => 'LowCov',    # pilot 1
		    'SRP000032'   => 'Trio',    # pilot 2
		    'SRP000033'   => 'Exon',    # pilot 3
		    'Pilot1'   => 'LowCov',    # pilot 1
		    'Pilot2'   => 'LowCov',    # pilot 1
		    'Pilot3'   => 'LowCov',    # pilot 1
		    );
	
	#read in the project spreadsheet with individuals
	open( PCSV, $pop_csv ) or die "Failed to open population csv sheet: $!\n";
	my %populations;
	while( <PCSV> )
	{
		chomp;
		next unless length( $_ ) > 0 && $_ !~ /^\s+$/ && $_ !~ /^,+$/;
		
		my @s = split( /,/, $_ );
		$populations{ $s[ 1 ] } = $s[ 5 ];
	}
	close( PCSV );
	
	my %exp_path; #within hierarchy
	
	croak "Can't find input index file\n" unless -f $index_file;
	
	open( IN, "$index_file" ) or die "Cannot open input index file\n";
	while( <IN> )
	{
		chomp;
		next if( $_ =~ /^\s+$/ || length( $_ ) == 0 || $_ =~ /^FASTQ_FILE/ );
		
		my @info = split( /\t+/, $_ );
		croak "Can't find study for lane ID: $info[ 3 ]\n" if( ! defined $study_vs_id{ $info[ 3 ] } );
		croak "Can't find population for lane ID: $info[ 9 ]\n" if( ! defined $populations{ $info[ 9 ] } );
		
		my @s = split( /\t+/, $_ );
		
		my $path = $hierarchy_dir.'/'.$study_vs_id{ $info[ 3 ] }.'-'.$populations{ $info[ 9 ] }.'/'.$info[ 9 ];
			
		#technology
		if( $info[ 12 ] =~ /454/ )
		{
			$path = $path.'/454';
		}
		elsif( $info[ 12 ] =~ /solid/i )
		{
			$path = $path.'/SOLID';
		}
		elsif( $info[ 12 ] =~ /solexa/i || $info[ 12 ] =~ /illumina/i )
		{
			$path = $path.'/SLX';
		}
		else
		{
			croak "Cannot determine sequencing technology: $info[ 12 ]\n";
		}
		$path .= '/'.$info[ 15 ];
		
		#see if we can find the lane fastq file
		my @fastq = <$path/*/*.fastq.gz>;
		my $alreadyImported = 0;
		my $fastqFile = $info[ 0 ];
		my $found = 0;
		foreach( @fastq )
		{
			chomp;
			if( basename( $_ ) eq $fastqFile )
			{
				#delete the lane directory
				my $dirname = dirname( $_ );
				system( "rm -rf $dirname" );
				print "Deleteing: $dirname\n";
				$found = 1;
				last;
			}
		}
		print "Cannot locate: $fastqFile\n" if( $found == 0 );
	}
	close( IN );
}

=pod
	A function to merge a hierarchy (e.g. a set of recalibrated files) into another hierarchy

sub mergeHierarchies
{
	croak "Usage: mergeHierarchies destinationHierarchy destinationHierarhcy_root_dir originHierarchy1_root_dir [originHierachy2_root_dir...]" unless @_ < 3;
	
	my $destinationHierarchy = shift;
	my $destinationRoot = shift;
	my $originHierarchy = shift;
	
	croak "Cannot find Destination Rood directory!" unless -d $destinationRoot;
	
	
}
=cut

sub pathmk {
   my @parts = File::Spec->splitdir( shift() );
   my $nofatal = shift;
   my $pth = $parts[0];
   my $zer = 0;
   if(!$pth) {
      $pth = File::Spec->catdir($parts[0],$parts[1]);
      $zer = 1;
   }
   my $DirPerms = '';
   for($zer..$#parts) {
      $DirPerms = oct($DirPerms) if substr($DirPerms,0,1) eq '0';
      mkdir($pth,$DirPerms) or return if !-d $pth && !$nofatal;
      mkdir($pth,$DirPerms) if !-d $pth && $nofatal;
      $pth = File::Spec->catdir($pth, $parts[$_ + 1]) unless $_ == $#parts;
   }
   1;
}

1;
