package HierarchyUtilities;
use strict;
use Carp;
use Utils;
use File::Basename;
use File::Spec; 
use File::Copy;
use Cwd qw(getcwd abs_path);
use Utils;
use File::Copy;

my $DFIND = '/software/solexa/bin/dfind';
my $MPSA_DOWNLOAD = '/software/solexa/bin/mpsa_download';
my $LSF_QUEUE = 'normal';

=pod

=head1 NAME

=head1 SYNOPSIS


=head1 ARGUMENTS

=head1 DESCRPTION

=head1 AUTHORS

Thomas Keane, I<tk2@sanger.ac.uk>
Petr Danecek, I<pd3@sanger.ac.uk>

=cut

=head1 METHODS

=head2 lane_info

    Arg [1]     : the lane path
    Returntype  : the hash with the following fields: 
                    project
                    sample
                    technology
                    library
                    lane
                    bwa_ref             .. e.g. /nfs/sf7/MOUSE/ref/NCBIM37_um
                    fa_ref              .. e.g. /nfs/sf7/MOUSE/ref/NCBIM37_um.fa
                    fai_ref             .. e.g. /nfs/sf7/MOUSE/ref/NCBIM37_um.fa.fai
                    dict_ref            .. e.g. /nfs/sf7/MOUSE/ref/NCBIM37_um.dict
                    ref_name            .. e.g. NCBIM37 (official name of the assembly)
                    snps                .. e.g. /nfs/sf7/MOUSE/ref/hapmap_mm9_matrix.snps.bin
                    genotype            .. e.g. 129S1_SvImJ
                    gtype_confidence    .. the expected likelihood ratio for genotype checks
                    insert_size         .. expected insert size, 0 for unpaired

=cut

sub lane_info
{
    my ($lane) = @_;

    # Get rid of multiple slashes // and the slash at the end dir/. Pathological cases
    #   like ./././../../ are not treated.
    #
    $lane =~ s{/+}{/};
    $lane =~ s{/$}{};
    my @items = split m{/}, $lane;
    if ( scalar @items < 5 ) { Utils::error("Wrong lane path: \"$lane\".\n") }

    my $info = 
    {
        'project'     => $items[-5],
        'sample'      => $items[-4],
        'technology'  => $items[-3],
        'library'     => $items[-2],
        'lane'        => $items[-1],

        'insert_size' => 500,
    };
	
    # This should be done differently in the future - the DB should tell us.
    if ( $$info{'project'} =~ /mouse/i or $$info{'project'} =~ /mice/i ) 
    { 
        $$info{'bwa_ref'}  = $ENV{MOUSE}.'/ref/NCBIM37_um';
        $$info{'fa_ref'}   = $ENV{MOUSE}.'/ref/NCBIM37_um.fa';
        $$info{'fai_ref'}  = $ENV{MOUSE}.'/ref/NCBIM37_um.fa.fai';
        $$info{'dict_ref'} = $ENV{MOUSE}.'/ref/NCBIM37_um.dict';
        $$info{'ref_name'} = 'NCBIM37';
        $$info{'snps'}     = '/nfs/sf7/MOUSE/ref/mousehapmap.snps.bin';

        my $genotype = 
        {
            '129P2_Mouse_Genome'         => '129P2_OlaHsD',
            '129S1_SvImJ_Mouse_Genome'   => '129S1_SvImJ',
            'AKR_J_Mouse_Genome'         => 'AKR_J',
            'A_J_Mouse_Genome'           => 'A_J',
            'BALBc_J_Mouse_Genome'       => 'BALB_cJ',
            'C3H_HeJ_Mouse_Genome'       => 'C3H_HeJ',
            'C57BL_6N_Mouse_Genome'      => 'C57BL_6NJ',
            'CAST_Ei_Mouse_Genome'       => 'CAST_EiJ',
            'CBA_J_Mouse_Genome'         => 'CBA_J',
            'DBA_2J_Mouse_Genome'        => 'DBA_2J',
            'LP_J_Mouse_Genome'          => 'LP_J',
            'NOD_Mouse_Genome'           => 'NOD_LtJ',
            'NZO_Mouse_Genome'           => 'NZO_HlLtJ',
            'PWK_Ph_Mouse_Genome'        => 'PWK_PhJ',
            'Spretus_Ei_Mouse_Genome'    => 'SPRET_EiJ',
            'WSB_Ei_Mouse_Genome'        => 'WSB_EiJ',
        };
        my $gtype_confidence = 
        {
            '129P2_Mouse_Genome'         => 5.0,
            '129S1_SvImJ_Mouse_Genome'   => 5.0,
            'AKR_J_Mouse_Genome'         => 5.0,
            'A_J_Mouse_Genome'           => 5.0,
            'BALBc_J_Mouse_Genome'       => 5.0,
            'C3H_HeJ_Mouse_Genome'       => 5.0,
            'C57BL_6N_Mouse_Genome'      => 5.0,
            'CAST_Ei_Mouse_Genome'       => 5.0,
            'CBA_J_Mouse_Genome'         => 5.0,
            'DBA_2J_Mouse_Genome'        => 5.0,
            'LP_J_Mouse_Genome'          => 5.0,
            'NOD_Mouse_Genome'           => 5.0,
            'NZO_Mouse_Genome'           => 5.0,
            'PWK_Ph_Mouse_Genome'        => 5.0,
            'Spretus_Ei_Mouse_Genome'    => 5.0,
            'WSB_Ei_Mouse_Genome'        => 5.0,
        };

        if ( exists($$genotype{$$info{'project'}}) )
        {
            $$info{'genotype'} = $$genotype{$$info{'project'}};
        }
        if ( exists($$gtype_confidence{$$info{'project'}}) )
        {
            $$info{'gtype_confidence'} = $$gtype_confidence{$$info{'project'}};
        }
    }
    elsif ( $$info{'sample'} =~ /^NA\d+$/  )   # something like NA18942
    {
        my $gender = `grep $$info{'sample'} $ENV{'G1K'}/ref/genders.txt`;
        chomp($gender);
        if ( $gender && $gender=~/\s+female/ )
        {
            $$info{'bwa_ref'}  = $ENV{G1K}.'/ref/human_b36_female';
            $$info{'fa_ref'}   = $ENV{G1K}.'/ref/human_b36_female.fa';
            $$info{'fai_ref'}  = $ENV{G1K}.'/ref/human_b36_female.fa.fai';
            $$info{'dict_ref'} = $ENV{G1K}.'/ref/human_b36_female.dict';
        }
        elsif ( $gender && $gender=~/\s+male/ )
        {
            $$info{'bwa_ref'}  = $ENV{G1K}.'/ref/human_b36_male';
            $$info{'fa_ref'}   = $ENV{G1K}.'/ref/human_b36_male.fa';
            $$info{'fai_ref'}  = $ENV{G1K}.'/ref/human_b36_male.fa.fai';
            $$info{'dict_ref'} = $ENV{G1K}.'/ref/human_b36_male.dict';
        }
        else
        {
            Utils::error("FIXME: $ENV{'G1K'}/ref/genders.txt not accessible or no $$info{'sample'} in there?\n");
        }
        $$info{'ref_name'} = 'NCBI36';
        $$info{'snps'}     = '/nfs/sf8/G1K/ref/snps/hapmap3.snps.bin';
        $$info{'genotype'} = $$info{sample};
        $$info{'gtype_confidence'} = 1.2;
    }
    else
    {
        Utils::error("FIXME: no reference data for the project $$info{'project'} .. $lane\n");
    }

    return $info;
}

=head2 getProjectDirectory

  Arg [1]    : seq index line
  Example    : getProjectDirectory( 'a seq index line string');
  Description: Figures out what is the name of the project directory (includes a nice hack for G1K)
  Returntype : string

=cut

sub getProjectDirectory
{
	croak "Usage: getProjectDirectory( 'seq_index_line' )" unless @_ == 1;
	
	my @s = split( /\t/, $_[ 0 ] );
	
	if( $s[ 3 ] eq 'SRP000031' || $s[ 3 ] eq 'SRP000032' || $s[ 3 ] eq 'SRP000033' || $s[ 4 ] =~ /1000Genomes/ )
	{
		my $samples_csv = '/nfs/sf8/G1K/meta-data/G1K_samples.txt';
		
		croak "Cant find G1K samples csv file: $samples_csv\n" unless -f $samples_csv;
		
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
		open( PCSV, $samples_csv ) or die "Failed to open population csv sheet: $!\n";
		my %populations;
		while( <PCSV> )
		{
			chomp;
			next unless length( $_ ) > 0 && $_ !~ /^\s+$/ && $_ !~ /^,+$/;
			
			my @s = split( /,/, $_ );
			$populations{ $s[ 1 ] } = $s[ 5 ];
		}
		close( PCSV );
		
		croak "Cant find study for lane ID: $s[ 3 ]\n" if( ! defined $study_vs_id{ $s[ 3 ] } );
		croak "Cant find population for lane ID: $s[ 9 ]\n" if( ! defined $populations{ $s[ 9 ] } );
		
		return $study_vs_id{ $s[ 3 ] }.'-'.$populations{ $s[ 9 ] };
	}
	else
	{
		$s[ 4 ] =~ s/\W+/_/g;
		
		croak "Study name not set: $_\n" unless $_ =~ /.+/;
		
		return $s[ 4 ];
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

sub importInternalData 
{
	croak "Usage: importInternalData project_list_file data_hierarchy_parent_directory analysis_hierarchy_parent_directory" unless @_ == 3;
    my $projectFile = shift;
    my $dhierarchyDir = shift;
    my $ahierarchyDir = shift;
	
    croak "Can't find data hierarchy directory\n" unless -d $dhierarchyDir;
    croak "Can't find analysis hierarchy directory\n" unless -d $ahierarchyDir;
    croak "Can't find projects file\n" unless -f $projectFile;
    
    my %projects;
    open( my $PROJS, "$projectFile" ) or die "Cannot open projects file\n";
    while(my $project =  <$PROJS> ) 
	{
		chomp $project;
		$projects{$project} = {};

		open(my $LIBS, q[-|], qq[$DFIND -project "$project" -libraries] ) or die "Cannot run dfind on project: $project\n";
		while(my $lib = <$LIBS> ) 
		{
			chomp $lib;
			$projects{$project}{$lib} = [];
			
			open my $LANES, q[-|], qq[$DFIND -project "$project" -library "$lib" -filetype fastq] or die "Can't get lanes for $project $lib: $!\n";
			while(my $lane = <$LANES> ) 
			{
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

=head2 importInternalData

  Arg [1]    : file of sanger lanes ids in format <run>_<lane>
  Arg [2]    : Data directory to build hierarchy in
  Arg [3]    : Mapping directory to build hierarchy in
  Example    : importInternalData( 'mouse.proj', '$G1K/MOUSE/DATA', '$G1K/MOUSE/MAPPING');
  Description: Imports all new sequence from the lanes listed into a hierarchy.  Existing fastq are skipped.
  Returntype : none

=cut

sub importInternalLanes
{
	croak "Usage: importInternalData lanes_fofn data_hierarchy_parent_directory analysis_hierarchy_parent_directory" unless @_ == 3;
    my $lanes = shift;
    my $dhierarchyDir = shift;
    my $ahierarchyDir = shift;
	
    croak "Can't find data hierarchy directory\n" unless -d $dhierarchyDir;
    croak "Can't find analysis hierarchy directory\n" unless -d $ahierarchyDir;
    croak "Can't find lanes file\n" unless -f $lanes;
	
	my $projectsHash;
	open( my $fh, $lanes ) or die "Cannot open lanes file: $!";
	while( <$fh>)
	{
		chomp;
		my $lane = $_;
		my @s = split( /_/, $_ );
		
		print "SEARCHING: $lane\n";
		open( my $ffh, q[-|], qq/$DFIND -file $s[ 0 ]_$s[ 1 ]_1.fastq/ ) or die "Cannot run dfind on lane $_\n";
		my $lib = '';
		my $proj = '';
		while( <$ffh> )
		{
			chomp;
			if( $_ =~ /^Project Tracking/ )
			{
				my @s1 = split( /\t/, $_ );
				$proj = $s1[ 1 ];
			}
			elsif( $_ =~ /Library Tracking/ )
			{
				my @s1 = split( /\t/, $_ );
				$lib = $s1[ 1 ];
			}
		}
		
		if( length( $proj ) == 0 || length( $lib ) == 0 )
		{
			#try the older hack of using _s_
			open( my $ffh, q[-|], qq/$DFIND -file $s[ 0 ]_s_$s[ 1 ].fastq/ ) or die "Cannot run dfind on lane $_\n";
			my $lib = '';
			my $proj = '';
			while( <$ffh> )
			{
				chomp;
				if( $_ =~ /^Project Tracking/ )
				{
					my @s1 = split( /\t/, $_ );
					$proj = $s1[ 1 ];
				}
				elsif( $_ =~ /Library Tracking/ )
				{
					my @s1 = split( /\t/, $_ );
					$lib = $s1[ 1 ];
				}
			}
			
			if( length( $proj ) == 0 || length( $lib ) == 0 )
			{
				print "WARNING: No proj/library info for $lane\n";
				next;
			}
			
			if( defined( $$projectsHash{ $proj }{ $lib } ) )
			{
				push( @{ $$projectsHash{ $proj }{ $lib } }, qq/$s[ 0 ]_s_$s[ 1 ].fastq/ );
			}
			else
			{
				$$projectsHash{ $proj }{ $lib } = [ qq/$s[ 0 ]_s_$s[ 1 ].fastq/ ];
			}
		}
		else
		{
			if( defined( $$projectsHash{ $proj }{ $lib } ) )
			{
				push( @{ $$projectsHash{ $proj }{ $lib } }, qq/$s[ 0 ]_$s[ 1 ]_1.fastq/ );
				push( @{ $$projectsHash{ $proj }{ $lib } }, qq/$s[ 0 ]_$s[ 1 ]_2.fastq/ );
			}
			else
			{
				$$projectsHash{ $proj }{ $lib } = [ qq/$s[ 0 ]_$s[ 1 ]_1.fastq/, qq/$s[ 0 ]_$s[ 1 ]_2.fastq/ ];
			}
		}
		close( $ffh );
	}
	
	foreach( keys( %{$projectsHash}) )
	{
		my $p = $_;
		foreach( keys( %{ $$projectsHash{ $_ } } ) )
		{
			my $l = $_;
			foreach( @{$$projectsHash{ $p }{ $l }} )
			{
				print qq/$p -> $l -> $_\n/;
			}
		}
	}
	&buildInternalHierarchy($projectsHash,$dhierarchyDir,$ahierarchyDir);
}


=head2 getMpsaFastq

  Arg [1]    : fastq name to retrieve
  Example    : my $dl_ok = getMpsaFastq ('1513_2_1.fastq');
  Description: retrieves a fastq from the mpsa and writes it to the same name in the current directory.  Checks the md5 of the fastq.
 
  If an error occurs or the md5 don't match, warns error and removes the fastq.
  Returntype : 1 if there were no errors, otherwise undef.

=cut

sub getMpsaFastq {
    croak "Usage: getMpsaFastq fastq" unless @_ == 1;
    my $fastq = shift;

    # jws 2009-02-11 - do this inline, and then do any parallelisation outside
    # this module
    my $cmd = "$MPSA_DOWNLOAD -c -f $fastq > $fastq 2> import.e";
    system( $cmd );
    my $ran_ok;
    if ( -s $fastq){
    # MD5 check
    my $mpsa_md5 = `$MPSA_DOWNLOAD -m -f $fastq|awk '{print \$1}'`;
    my $dl_md5 = `md5sum $fastq|awk '{print \$1}'`;
    if ($mpsa_md5 && $dl_md5){
        if($mpsa_md5 eq $dl_md5){
        $ran_ok = 1;
        }
        else {
        warn "md5 of $fastq did not match that in MPSA\n";
        $ran_ok = 0;
        }
    }
    else {
        warn "Can't retrieve md5 of $fastq\n";
        $ran_ok = 0;
    }
    }
    else {
    warn "Error retrieving $fastq with mpsa_download: $!\n";
    $ran_ok = 0;
    }

    unless ($ran_ok){
    unlink($fastq);
    }
    return $ran_ok;

}
 


=head2 checkInternalFastq

  Arg [1]    : DATA lane directory to check
  Example    : my $is_ok = checkInternalFastq('/lustre/sf4/1kgenomes/G1K/SANGER/1000Genomes_B1_TOS/NA20534/SLX/NA20534_TOS_1/1513_2');
  Description: compares NPG/MPSA fastqcheck to fastqcheck(s) for retrieved fastq.  Handles split fastq as well as unsplit.  Note, uses dfind to retrieve the NPG fastqcheck file(s).
  Returntype : 1 if fastqcheck counts are the same.

=cut

sub checkInternalFastq {
    croak "Usage: checkInternalFastq hierarchy_fastq_dir" unless @_ == 1;
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
      $new_seqs = $new_seqs/2;  # should be same number of seqs in both
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

sub getSeqAndLengthFromFastqcheck 
{
    my $fastqcheck = shift;
    open (my $FQC, $fastqcheck) or croak "Can't open $fastqcheck: $!\n";
    my $header = <$FQC>;
    close $FQC;
    my ($seqcount,undef,$length) = split /\s+/,$header;
    return ($seqcount,$length);
}

=head2 deleteIndexLanes

  Arg [1]    : G1K samples csv file
  Arg [2]    : index file of lanes to delete
  Arg [3]    : root directory of hierarchy
  Example    : deleteIndexLanes( '$G1K/ref/G1K_samples.txt', 'index.file', 'hier_parent_dir' );"
  Description: Function to delete lanes from a hierarchy.
  Returntype : none

=cut
sub deleteIndexLanes
{
  croak "Usage: deleteIndexLanes samples_csv index_file hierarchy_root_dir" unless @_ == 3;
  
  my $pop_csv = shift;
  my $index_file = shift;
  my $hierarchy_dir = shift;
  
  croak "Cant find populations_csv file" unless -f $pop_csv;
  croak "Cant find index file" unless -f $index_file;
  croak "Cant find hierarchy parent directory" unless -d $hierarchy_dir;
  
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
  
  croak "Cant find input index file\n" unless -f $index_file;
  
  open( IN, "$index_file" ) or die "Cannot open input index file\n";
  while( <IN> )
  {
    chomp;
    next if( $_ =~ /^\s+$/ || length( $_ ) == 0 || $_ =~ /^FASTQ_FILE/ );
    
    my @info = split( /\t+/, $_ );
    croak "Cant find study for lane ID: $info[ 3 ]\n" if( ! defined $study_vs_id{ $info[ 3 ] } );
    croak "Cant find population for lane ID: $info[ 9 ]\n" if( ! defined $populations{ $info[ 9 ] } );
    
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

sub syncIndexFile
{
    croak "Usage: syncIndexFile index_file lanes_fofn dataRootDir mappingRootDir" unless @_ == 4;
    
    my $index_file = shift;
    my $lanesfofn = shift;
    my $dRoot = shift;
    my $mRoot =  shift;
    
    croak "Cant find index file" unless -f $index_file;
    croak "Cant find lanes file" unless -f $lanesfofn;
    croak "Cant find data root directory" unless -d $dRoot;
    croak "Cant find mapping root directory" unless -d $mRoot;
    
    my $base = basename( $index_file );
    
    my %lanesInHierarchy;
    open( LANES, $lanesfofn ) or die "Cannot open lanes fofn: $lanesfofn\n";
    while( <LANES> )
    {
        chomp;
        
        my $id = basename( $_ );
        
        $lanesInHierarchy{ $id } = $_;
    }
    close( LANES );
    
    open( LIBMIS, ">$base.lib.missing.index" ) or die $!;
    open( FMIS, ">$base.fastq.missing.index" ) or die $!;
    open( LMIS, ">$base.lane.missing.index" ) or die $!;
    #open( LM, ">$base.lane.missing.index" ) or die $!;
    open( WITH, ">remove_withdrawn.sh" ) or die $!;
    open( INCOMPLETE, ">remove_incomplete.sh" ) or die $!;
    open( IN, $index_file ) or die "Cannot open index file\n";
    open( EXTRA, ">remove_extra_lanes.sh" ) or die $!;
    my %index_ids;
    my %withdrawn;
    while( <IN> )
    {
        chomp;
        
        next unless $_ !~ /FASTQ/;
        
        my @paths = @{ indexLineToPaths( $_ ) };
        my @s = split( /\t/, $_ );
        
        #check if the lane has been withdrawn
        if( $s[ 20 ] == 1 )
        {
            if( -d $dRoot."/".$paths[ 1 ] )
            {
                print "WITHDRAWN: $dRoot/$paths[ 1 ]\n";
                
                #delete the directory in the mapping and data hierarchies
                print WITH "rm -rf $dRoot/$paths[ 1 ] $mRoot/$paths[ 1 ]\n";
                $withdrawn{ $s[ 2 ] } = 1;
            }
            
            next;
        }
        
        $index_ids{ $s[ 2 ] } = 1;
        
        if( ! -d $dRoot."/".$paths[ 0 ] )
        {
            if( defined( $lanesInHierarchy{ $s[ 2 ] } ) )
            {
                print "LIBRARY_SWAP: $lanesInHierarchy{ $s[ 2 ] } $_\n";
                next;
            }
            else
            {
                print "LIBRARY_NOT_FOUND: $dRoot/$paths[ 0 ]\n";
                print LIBMIS $_."\n";
                next;
            }
        }
        
        if( ! -d $dRoot."/".$paths[ 1 ] )
        {
            #check the file has not been withdrawn
            if( $s[ 20 ] == 0 )
            {
                if( defined( $lanesInHierarchy{ $s[ 2 ] } ) )
                {
                    print "LANE_SWAP: $lanesInHierarchy{ $s[ 2 ] } $_\n";
                    next;
                }
                else
                {
                    print "LANE_NOT_FOUND: $_\t$dRoot/$paths[ 1 ]\n";
                    print LMIS $_."\n";
                }
            }
        }
        elsif( ! -f $dRoot."/".$paths[ 2 ] )
        {
            #check the file has not been withdrawn
            if( $s[ 20 ] == 0 )
            {
                if( defined( $lanesInHierarchy{ $s[ 2 ] } ) )
                {
                    print "FASTQ_SWAP: $lanesInHierarchy{ $s[ 2 ] } $_\n";
                    next;
                }
                else
                {
                    print "FASTQ_NOT_FOUND: $_\t$dRoot/$paths[ 2 ]\n";
                    print INCOMPLETE "rm -rf $dRoot/$paths[ 1 ] $mRoot/$paths[ 1 ]\n";
                    print FMIS $_."\n";
                }
            }
        }
        else
        {
            print "AGREE: $dRoot/$paths[ 1 ]\n";
        }
    }
    
    foreach( keys( %lanesInHierarchy ) )
    {
        if( ! defined( $index_ids{ $_ } ) && ! defined( $withdrawn{ $_ } ) )
        {
            print "EXTRA_LANE: $_\n";
            print EXTRA "$_\n";
        }
    }
    
    close( IN );
    close( LIBMIS );
    close( LMIS );
    close( FMIS );
    close( WITH );
    close( INCOMPLETE );
    close( EXTRA );
}

#takes an index line and returns 3 paths (the library path, lane path, fastq path)
sub indexLineToPaths
{
  croak "Usage: indexLineToPath index_line" unless @_ == 1;
  
  my $index_line = shift;
  
  my @s = split( /\t/, $index_line );
  
  my @paths;
  
  my $individualPath = determineIndividualProjectDir( $s[ 9 ], $s[ 3 ] )."/".$s[ 9 ];
  
  #technology
  if( $s[ 12 ] =~ /454|roche/i )
  {
    $individualPath = $individualPath.'/454';
  }
  elsif( $s[ 12 ] =~ /solid/i )
  {
    $individualPath = $individualPath.'/SOLID';
  }
  elsif( $s[ 12 ] =~ /solexa|illumina/i )
  {
    $individualPath = $individualPath.'/SLX';
  }
  else
  {
    croak "Cannot determine sequencing technology: $s[ 12 ]\n";
  }
  
  $s[ 14 ] =~ s/\s+/_/g;
  
  my $libraryPath = $individualPath."/".$s[ 14 ];
  push( @paths, $libraryPath );
  
  my $lanePath = $libraryPath.'/'.$s[ 2 ];
  push( @paths, $lanePath );
  
  my $fastqPath = $lanePath.'/'.basename( $s[ 0 ] );
  push( @paths, $fastqPath );
  
  return \@paths;
}

my $SAMPLES_CSV = $ENV{ 'G1K' }.'/meta-data/G1K_samples.txt';
my %G1K_SAMPLES;
sub determineIndividualProjectDir
{
  croak "Usage: determineIndividualProject individual proj_code" unless @_ == 2;
  
  my $individual = shift;
  my $proj_code =  shift;
  
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
  if( keys( %G1K_SAMPLES ) == 0 )
  {
    open( PCSV, $SAMPLES_CSV ) or die "Failed to open population csv sheet: $!\n";
    while( <PCSV> )
    {
      chomp;
      next unless length( $_ ) > 0 && $_ !~ /^\s+$/ && $_ !~ /^,+$/;
      
      my @s = split( /,/, $_ );
      $G1K_SAMPLES{ $s[ 1 ] } = $s[ 5 ];
    }
    close( PCSV );
  }
  
  if( ! defined $study_vs_id{ $proj_code } )
  {
    croak "Undefined project code: $proj_code\n";
  }
  
  if( ! defined $G1K_SAMPLES{ $individual } )
  {
    croak "Undefined individual: $individual"."x\n";
  }
  
  return $study_vs_id{ $proj_code }.'-'.$G1K_SAMPLES{ $individual };
}

=head2 markFailedGenotypeLanes

    Arg [1]    : file of lane accessions
    Arg [2]    : hierarchy root directory
    Arg [3]    : sequence.index file
    Example    : markFailedGenotypeLanes( 'accessions.fofn', '$G1K/MOUSE/DATA', 'sequence.index');
    Description: Takes a list of lane accessions and marks them (i.e. in meta.info) as being genotype failed
    Returntype : none
=cut
sub markFailedGenotypeLanes
{
  croak "Usage: markBadLanes fileOfAccessions data_root index_file" unless @_ == 3;
  
  my $accessions = shift;
  my $dRoot =  shift;
  my $indexF = shift;
  
  croak "Cant find accessions file: $accessions\n" unless -f $accessions;
  croak "Cant find index file: $indexF\n" unless -f $indexF;
  croak "Cant find data root directory: $dRoot\n" unless -d $dRoot;
  
  my %accessions;
  open( A, $accessions ) or die "failed to open accessions file: $accessions $!\n";
  while( <A> )
  {
    chomp;
    my $acc = (split( /\t/, $_))[ 0 ];
    $accessions{ $acc } = 1;
  }
  close( A );

  open( IN, $indexF ) or die "Cant open index file: $indexF\n";
  while( <IN> )
  {
    chomp;
    my $acc = (split( /\t/, $_))[ 2 ];

    if( defined( $accessions{ $acc } ) )
    {
      my @paths = @{ indexLineToPaths( $_ ) };
      
      my $meta = $dRoot."/".$paths[ 1 ]."/meta.info";
      if( -f $meta )
      {
        open( M, ">>$meta" ) or die "Cant open $meta\n";
        print M "genotype:fail\n";
        close( M );
      }
    }
  }
  close( IN );
}

=head2
dfd
=cut

sub symLinkHierarchyFiles
{
	croak "Usage: symLinkHierarchyFiles original_rootDir files_fofn release_root_directory hard|soft|cp" unless @_ == 4;
    
    my $original_root = shift;
    my $files_fofn = shift;
    my $release_root = shift;
	my $type = shift;
    
    croak "Cant find the lanes_fofn file: $files_fofn\n" unless -f $files_fofn;
    croak "Cant find root directory: $original_root\n" unless -d $original_root;
    croak "Cant find release root directory: $release_root\n" unless -d $release_root;
	croak "Must specify either hard or soft as link type" unless $type eq 'soft' || $type eq 'hard' || $type eq 'copy';
    
    my %studyToRoot;
    $studyToRoot{ 'LowCov' } = "SRP000031";
    $studyToRoot{ 'Trio' } = "SRP000032";
    $studyToRoot{ 'Exon' } = "SRP000033";
	
	open( LANES, $files_fofn ) or die "Cant open bam fofn file: $!\n";
    while( <LANES> )
    {
        chomp;
        
        my $file = $_;
        
        my $originalPath = qq[$original_root/$file];
        my $destinationPath = qq[$release_root/$file];
        my $destinationDir = dirname( $destinationPath );
        
        if( ! -f $originalPath )
        {
            print "Cant find original lane file: ".$originalPath."\n";
            next;
        }
        
        if( ! -d $destinationDir )
        {
            system( "mkdir -p $destinationDir" ) #or die "Failed to make release directory path: $release_root/$relativeDir\n";
        }
        
        if( ! -l $destinationPath && ! -f $destinationPath )
        {
            print "Linking file in: $destinationPath to $originalPath\n";
			if( $type eq 'soft' )
			{
				symlink( $originalPath, $destinationPath ) or die "Failed to link file in: $destinationPath to $originalPath: $!\n";
			}
			elsif( $type eq 'hard' )
			{
				link( $originalPath, $destinationPath ) or die "Failed to link file in: $destinationPath to $originalPath: $!\n";
			}
			else
			{
				copy( $originalPath, $destinationPath ) or die "Failed to copy file in: $destinationPath to $originalPath: $!\n";
			}
        }
    }
    close( LANES );
}

=head2 buildReleaseHierarchy

    Arg [1]    : original hierarchy root directory
    Arg [2]    : file of bam files (path relative to hierarchy root dir)
    Arg [3]    : root directory of release hierarchy
    Example    : buildReleaseHierarchy( '$G1K/MOUSE/MAPPING', 'bam.fofn', '$G1K/MOUSE/RELEASE-01');
    Description: Builds a sideways hierarchy and links in the bam files listed in the fofn
    Returntype : none
=cut
sub buildReleaseHierarchy
{
    croak "Usage: buildReleaseHierarchy original_rootDir files_fofn release_root_directory index_file" unless @_ == 4;
    
    my $original_root = shift;
    my $files_fofn = shift;
    my $release_root = shift;
    my $indexF = shift;
    
    croak "Cant find the lanes_fofn file: $files_fofn\n" unless -f $files_fofn;
    croak "Cant find root directory: $original_root\n" unless -d $original_root;
    croak "Cant find release root directory: $release_root\n" unless -d $release_root;
    croak "Cant find index file: $indexF\n" unless -f $indexF;
    
    my %studyToRoot;
    $studyToRoot{ 'LowCov' } = "SRP000031";
    $studyToRoot{ 'Trio' } = "SRP000032";
    $studyToRoot{ 'Exon' } = "SRP000033";
    
    my %indexFile;
    open( LANES, $indexF ) or die "Cannot open lanes fofn: $indexF\n";
    while( <LANES> )
    {
        chomp;
        
        my @s = split( /\t/, $_ );
        
        if( $s[ 20 ] == 0 )
        {
            $indexFile{ $s[ 2 ] } = \@s;
        }
        else
        {
            print "Skipping withdrawn lane: $s[ 2 ]\n";
        }
    }
    close( LANES );
    
    open( LANES, $files_fofn ) or die "Cant open bam fofn file: $!\n";
    while( <LANES> )
    {
        chomp;
        
        my $file = $_;
        
        my $originalPath = qq[$original_root/$file];
        my $destinationPath = qq[$release_root/$file];
        my $destinationDir = dirname( $destinationPath );
        
        if( ! -f $originalPath )
        {
            print "Cant find original lane file: ".$originalPath."\n";
            next;
        }
        
        #verify that the genotype is correct
        #if( $file =~ /.*\/(NA[0-9]+)\/.*\/([SRR|ERR][0-9]+)\/.*/ ){print "boo!";}exit;
        $file =~ /(LowCov|Trio|Exon)-.*\/(NA[0-9]+)\/.*\/(SRR[0-9]+|ERR[0-9]+).*/;
        my $ind = $2;
        my $acc = $3;
        my $laneStudyRoot = $1;
        
        if( ! defined( $indexFile{ $acc } ) )
        {
            print "No entry in index file for accession: $acc\n";
            print "$file\n";
            next;
        }
        elsif( $indexFile{ $acc }[ 9 ] ne $ind )
        {
            print "IND_SWAP: $_ -> @{ $indexFile{ $acc } }\n";
            next;
        }
        
        if( $studyToRoot{ $laneStudyRoot } ne $indexFile{ $acc }[ 3 ] )
        {
            print "STUDY_SWAP: $_ -> @{ $indexFile{ $acc } }\n";
            #next;
        }
        
        if( ! -d $destinationDir )
        {
            system( "mkdir -p $destinationDir" ) #or die "Failed to make release directory path: $release_root/$relativeDir\n";
        }
        
        if( ! -l $destinationPath && ! -f $destinationPath )
        {
            print "Linking file in: $destinationPath to $originalPath\n";
            symlink( $originalPath, $destinationPath ) or die "Failed to link file in: $destinationPath to $originalPath\n";
        }
    }
    close( LANES );
}

=head2 getFastqInfo

    Arg [1]    : lane directory
    Returntype : 2d array reference where each array consists of [ fastq name ][ read length ][ num reads ][ num bases ]
                where this information is derived from the meta.info file and fastqcheck files
=cut
sub getFastqInfo
{
	croak "Usage: getFastqInfo laneDir" unless @_ == 1;
	
	my $laneDir = shift;
	
	croak "Cant find lane directory: $laneDir\n" unless -d $laneDir;
	
	my $meta_info_file = File::Spec->catfile($laneDir, 'meta.info');
	
	croak "Cant find meta.info file in $laneDir\n" unless -f $meta_info_file;
	
	my $lane_read0 = '';
	my $lane_read1 = '';
	my $lane_read2 = '';
	if( -f $meta_info_file || -l $meta_info_file )
	{
		$lane_read0=`grep read0 $meta_info_file | head -1 | awk -F: '{print \$2}'`;
		$lane_read1=`grep read1 $meta_info_file | head -1 | awk -F: '{print \$2}'`;
		$lane_read2=`grep read2 $meta_info_file | head -1 | awk -F: '{print \$2}'`;
		chomp( $lane_read0 );
		chomp( $lane_read1 );
		chomp( $lane_read2 );
		$lane_read0 = File::Spec->catfile($laneDir, $lane_read0) if $lane_read0;
		$lane_read1 = File::Spec->catfile($laneDir, $lane_read1) if $lane_read1;
		$lane_read2 = File::Spec->catfile($laneDir, $lane_read2) if $lane_read2;
	}
	
	my @reads;
	
	#get the read lengths from the fastqcheck
	if( -f $lane_read0 )
	{
		my $length0 = `head -1 $lane_read0.fastqcheck | awk '{print \$6}'`;chomp( $length0 );
		my $num_reads = `head -1 $lane_read0.fastqcheck | awk '{print \$1}'`;chomp( $num_reads );
		my $num_bases = `head -1 $lane_read0.fastqcheck | awk '{print \$3}'`;chomp( $num_bases );
		$reads[ 0 ] = [ basename($lane_read0), $length0, $num_reads, $num_bases ];
	}
	else
	{
		$reads[ 0 ] = [ '', 0, 0, 0 ];
	}
	
	if( -f $lane_read1 )
	{
		my $length1 = `head -1 $lane_read1.fastqcheck | awk '{print \$6}'`;chomp( $length1 );
		my $num_reads = `head -1 $lane_read1.fastqcheck | awk '{print \$1}'`;chomp( $num_reads );
		my $num_bases = `head -1 $lane_read1.fastqcheck | awk '{print \$3}'`;chomp( $num_bases );
		$reads[ 1 ] = [ basename($lane_read1), $length1, $num_reads, $num_bases ];
	}
	else
	{
		$reads[ 1 ] = [ '', 0, 0, 0 ];
	}
	
	if( -f $lane_read2 )
	{
		my $length2 = `head -1 $lane_read2.fastqcheck | awk '{print \$6}'`;chomp( $length2 );
		my $num_reads = `head -1 $lane_read2.fastqcheck | awk '{print \$1}'`;chomp( $num_reads );
		my $num_bases = `head -1 $lane_read2.fastqcheck | awk '{print \$3}'`;chomp( $num_bases );
		$reads[ 2 ] = [ basename($lane_read2), $length2, $num_reads, $num_bases ];
	}
	else
	{
		$reads[ 2 ] = [ '', 0, 0, 0 ];
	}
	
	return \@reads;
}

1;
