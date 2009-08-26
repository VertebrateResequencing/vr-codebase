package HierarchyUtilities;
use strict;
use Carp;
use Utils;
use File::Basename;
use File::Spec; 
use File::Copy;
use Cwd qw(getcwd abs_path);
use Utils;

use AssemblyTools;

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
                    md5_ref             .. e.g. 28f4ff5cf14f5931d0d531a901236378
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
        $$info{'bwa_ref'}  = '/nfs/sf7/MOUSE/ref/NCBIM37_um';
        $$info{'fa_ref'}   = '/nfs/sf7/MOUSE/ref/NCBIM37_um.fa';
        $$info{'fai_ref'}  = '/nfs/sf7/MOUSE/ref/NCBIM37_um.fa.fai';
        $$info{'dict_ref'}  = '/lustre/scratch103/sanger/team145/mouse/ref/NCBIM37_um.dict';
        $$info{'ref_name'} = 'NCBIM37';
        $$info{'snps'}     = '/nfs/sf7/MOUSE/ref/mousehapmap.snps.bin';
        $$info{'md5_ref'}  = '36a352ec67f958c40f19a9cf6ceb0d1e';

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
            $$info{'bwa_ref'}  = '/nfs/sf8/G1K/ref/human_b36_female';
            $$info{'fa_ref'}   = '/nfs/sf8/G1K/ref/human_b36_female.fa';
            $$info{'fai_ref'}  = '/nfs/sf8/G1K/ref/human_b36_female.fa.fai';
            $$info{'dict_ref'}  = '/lustre/scratch103/sanger/team145/g1k/ref/human_b36_female.dict';
            $$info{'md5_ref'}  = '28f4ff5cf14f5931d0d531a901236378';
        }
        elsif ( $gender && $gender=~/\s+male/ )
        {
            $$info{'bwa_ref'}  = '/nfs/sf8/G1K/ref/human_b36_male';
            $$info{'fa_ref'}   = '/nfs/sf8/G1K/ref/human_b36_male.fa';
            $$info{'fai_ref'}  = '/nfs/sf8/G1K/ref/human_b36_male.fa.fai';
            $$info{'dict_ref'}  = '/lustre/scratch103/sanger/team145/g1k/ref/human_b36_male.dict';
            $$info{'md5_ref'}  = '7bcc140d6728a6c9fe6ed411ee633862';
        }
        else
        {
            Utils::error("FIXME: $ENV{'G1K'}/ref/genders.txt not accessible or no $$info{'sample'} in there?\n");
        }
        $$info{'ref_name'} = 'NCBI36';
        $$info{'snps'}     = '/nfs/sf8/G1K/ref/snps/hapmap3.snps.bin';
        $$info{'genotype'} = $$info{sample};
        $$info{'gtype_confidence'} = 1.4;
    }
    else
    {
        Utils::error("FIXME: no reference data for the project $$info{'project'} .. $lane\n");
    }

    return $info;
}

=head2 importExternalData

  Arg [1]	: sequence index file
  Arg [2]	: data hierarchy parent directory
  Arg [3]	: mapping hierarchy parent directory
  Arg [4]	: directory of fastq files
  Arg [5]	: mode (link or copy)
  Example	: importExternalData( 'sequence.index', '$G1K/MOUSE/DATA', '$G1K/MOUSE/MAPPING', '/lustre/data', 'link' );
  Description: Imports the files listed in the sequence index file into an existing hierarchy.
  Returntype : none

=cut

sub importExternalData
{
  croak "Usage: importExternalData seq_index hierarchy_parent_directory mapping_hierarchy_directory fastq_directory link|copy" unless @_ == 5;
  
  my $index_file = shift;
  my $data_hierarchy_dir = abs_path(shift);
  my $mapping_hierarchy_dir = abs_path(shift);
  my $fastqDir = abs_path(shift);
  my $method = shift;
  
  croak "Cant find index file" unless -f $index_file;
  croak "Cant find fastq directory" unless -d $fastqDir;
  croak "Cant find data hierarchy parent directory" unless -d $data_hierarchy_dir;
  croak "Cant find mapping hierarchy parent directory" unless -d $mapping_hierarchy_dir;
  croak "Method of import not specified correctly: " unless $method eq "link" || $method eq "copy";
  
  my %exp_read0;
  my %exp_read1;
  my %exp_read2;
  my %data_exp_path; #directory within data hierarchy
  my %mapping_exp_path; #mapping directory
  my %inserts; #insert sizes per accession
  
  my %srrCount;
  
  open( LCSV, $index_file ) or die "Failed to open lane index file: $!\n";
  while( <LCSV> )
  {
    chomp;
    my $original = $_;
    
    next if( $_ =~ /^\s+$/ || length( $_ ) == 0 || $_ =~ /^FASTQ_FILE/ );
    
    my $project_directory = getProjectDirectory( $_ );
    
    my @info = split( /\t/, $original );
    
    if( $info[ 18 ] ne "PAIRED" && $info[ 18 ] ne "SINGLE" )
    {
      print "Skipping line: $original\n";
      my $i = 0;
      foreach( @info )
      {
        print $i."->x".$info[ $i ]."x\n";
        $i ++;
      }
            
      exit;
      next;
    }
    
    #check SRR is unique
    if( defined( $srrCount{ $info[ 2 ] } ) )
    {
      if( ( $srrCount{ $info[ 2 ] } == 3 && $info[ 12 ] =~ /.*[roche|454].*/ ) || $info[ 19 ] eq 'SINGLE' )
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
    if( $info[ 18 ] eq 'SINGLE' || ( $info[ 18 ] eq 'PAIRED' && length( $info[ 19 ] ) == 0 ) )
    {
      $exp_read0{ $info[ 2 ] } = basename( $info[ 0 ] );
    }
    elsif( $info[ 18 ] eq 'PAIRED' )
    {
      if( ! defined( $exp_read1{ $info[ 2 ] } ) )
      {
        $exp_read1{ $info[ 2 ] } = basename( $info[ 0 ] );
      }
      else
      {
        $exp_read2{ $info[ 2 ] } = basename( $info[ 0 ] );
      }
    }
    
    my $dpath = $data_hierarchy_dir.'/'.$project_directory.'/'.$info[ 9 ];
    my $mpath = $mapping_hierarchy_dir.'/'.$project_directory.'/'.$info[ 9 ];
    
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
    
    #transform spaces to underscores
    $info[ 14 ] =~ s/\s+/_/g;
    
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
	if( $info[ 17 ] =~ /\d+/ )
	{
		$inserts{ $info[ 2 ] } = $info[ 17 ];
	}
  }
  close( LCSV );
  
  foreach( keys( %data_exp_path ) )
  {
    #check the fastq's exist
    if( defined( $exp_read1{ $_ } ) && ! defined( $exp_read2{ $_ } ) )
    {
      #turn it into an unpaired entry
      $exp_read0{ $_ } = $exp_read1{ $_ };
      delete( $exp_read1{ $_ } );
      print "No read2 entry so making read unpaired: $exp_read0{ $_ }\n";
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
    elsif( defined( $exp_read0{ $_ } ) && ! -f $fastqDir.'/'.$exp_read0{ $_ } )
    {
      print "WARNING: Failed on $_: Failed to read file $exp_read0{ $_ }\n";
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
    open( META, ">$pdirName/meta.info" ) or die "Cant open meta.info file in $pdirName\n";
    
    #copy the files into the directory
    if( defined( $exp_read0{ $_ } ) && ! -l $pdirName.'/'.$exp_read0{ $_ } && ! -f $pdirName.'/'.$exp_read0{ $_ } )
    {
        if( $method eq "copy" )
        {
            copy( $fastqDir.'/'.$exp_read0{ $_ }, $pdirName.'/'.$exp_read0{ $_ } ) or die "Failed to copy file: ".$exp_read0{ $_ };
            system( 'bsub -q yesterday -o fastqcheck.o -e fastqcheck.e "zcat '.$pdirName.'/'.$exp_read0{ $_ }.' | awk \'{print \$1}\' | fastqcheck > '.$pdirName.'/'.$exp_read0{ $_ }.'.fastqcheck; ln -s '.$pdirName.'/'.$exp_read0{ $_ }.'.fastqcheck '.$mdirName.'/'.$exp_read0{ $_ }.'.fastqcheck"' );
        }
        elsif( $method eq "link" )
        {
            symlink( $fastqDir.'/'.$exp_read0{ $_ }, $pdirName.'/'.$exp_read0{ $_ } ) or die "Failed to link file: ".$exp_read0{ $_ };
			system( 'bsub -q yesterday -o fastqcheck.o -e fastqcheck.e "zcat '.$pdirName.'/'.$exp_read0{ $_ }.' | awk \'{print \$1}\' | fastqcheck > '.$pdirName.'/'.$exp_read0{ $_ }.'.fastqcheck; ln -s '.$pdirName.'/'.$exp_read0{ $_ }.'.fastqcheck '.$mdirName.'/'.$exp_read0{ $_ }.'.fastqcheck"' );
        }

		#work out the clip points (if any)
		my $cmd = qq[bsub -q small -o import.o -e import.e perl -w -e "use AssemblyTools;AssemblyTools::writeClipPointMeta( '$pdirName/$exp_read0{ $_ }', '$pdirName/meta.info', '0' );"];
		system( $cmd );
		
		symlink( $pdirName.'/'.$exp_read0{ $_ }, $mdirName.'/'.$exp_read0{ $_ } );
        print META "read0:".$exp_read0{ $_ }."\n";
    }
    
    if( defined( $exp_read1{ $_ } ) && ! -l $pdirName.'/'.$exp_read1{ $_ } && ! -f $pdirName.'/'.$exp_read1{ $_ } )
    {
        if( $method eq "copy" )
        {
            copy( $fastqDir.'/'.$exp_read1{ $_ }, $pdirName.'/'.$exp_read1{ $_ } ) or die "Failed to copy file: ".$exp_read1{ $_ };
            system( 'bsub -q yesterday -o fastqcheck.o -e fastqcheck.e "zcat '.$pdirName.'/'.$exp_read1{ $_ }.' | awk \'{print \$1}\' | fastqcheck > '.$pdirName.'/'.$exp_read1{ $_ }.'.fastqcheck; ln -s '.$pdirName.'/'.$exp_read1{ $_ }.'.fastqcheck '.$mdirName.'/'.$exp_read1{ $_ }.'.fastqcheck"' );
        }
        elsif( $method eq "link" )
        {
            symlink( $fastqDir.'/'.$exp_read1{ $_ }, $pdirName.'/'.$exp_read1{ $_ } ) or die "Failed to link file: ".$exp_read1{ $_ };
			system( 'bsub -q yesterday -o fastqcheck.o -e fastqcheck.e "zcat '.$pdirName.'/'.$exp_read1{ $_ }.' | awk \'{print \$1}\' | fastqcheck > '.$pdirName.'/'.$exp_read1{ $_ }.'.fastqcheck; ln -s '.$pdirName.'/'.$exp_read1{ $_ }.'.fastqcheck '.$mdirName.'/'.$exp_read1{ $_ }.'.fastqcheck"' );
        }
        
		#work out the clip points (if any)
		my $cmd = qq[bsub -q small -o import.o -e import.e perl -w -e "use AssemblyTools;AssemblyTools::writeClipPointMeta( '$pdirName/$exp_read1{ $_ }', '$pdirName/meta.info', '1' );"];
		system( $cmd );

        symlink( $pdirName.'/'.$exp_read1{ $_ }, $mdirName.'/'.$exp_read1{ $_ } ) or die "Failed to sym link read into mapping directory: ".$mdirName.'/'.$exp_read1{ $_ };
        print META "read1:".$exp_read1{ $_ }."\n";
		
		#write the insert size information
		if( defined( $inserts{ $_ } ) )
		{
			print META "insert:$inserts{ $_ }\n";
		}
    }
    
    if( defined( $exp_read2{ $_ } ) && ! -l $pdirName.'/'.$exp_read2{ $_ } && ! -f $pdirName.'/'.$exp_read2{ $_ } )
    {
        if( $method eq "copy" )
        {
            copy( $fastqDir.'/'.$exp_read2{ $_ }, $pdirName.'/'.$exp_read2{ $_ } ) or die "Failed to copy file: ".$exp_read2{ $_ };
            system( 'bsub -q yesterday -o fastqcheck.o -e fastqcheck.e "zcat '.$pdirName.'/'.$exp_read2{ $_ }.' | awk \'{print \$1}\' | fastqcheck > '.$pdirName.'/'.$exp_read2{ $_ }.'.fastqcheck; ln -s '.$pdirName.'/'.$exp_read2{ $_ }.'.fastqcheck '.$mdirName.'/'.$exp_read2{ $_ }.'.fastqcheck"' );
        }
        elsif( $method eq "link" )
        {
            symlink( $fastqDir.'/'.$exp_read2{ $_ }, $pdirName.'/'.$exp_read2{ $_ } ) or die "Failed to symlink file: ".$exp_read2{ $_ };
			system( 'bsub -q yesterday -o fastqcheck.o -e fastqcheck.e "zcat '.$pdirName.'/'.$exp_read2{ $_ }.' | awk \'{print \$1}\' | fastqcheck > '.$pdirName.'/'.$exp_read2{ $_ }.'.fastqcheck; ln -s '.$pdirName.'/'.$exp_read2{ $_ }.'.fastqcheck '.$mdirName.'/'.$exp_read2{ $_ }.'.fastqcheck"' );
        }

		#work out the clip points (if any)
		my $cmd = qq[bsub -q small -o import.o -e import.e perl -w -e "use AssemblyTools;AssemblyTools::writeClipPointMeta( '$pdirName/$exp_read2{ $_ }', '$pdirName/meta.info', '2' );"];
		system( $cmd );
        
        symlink( $pdirName.'/'.$exp_read2{ $_ }, $mdirName.'/'.$exp_read2{ $_ } ) or die "Failed to sym link read into mapping directory: ".$mdirName.'/'.$exp_read2{ $_ };
        print META "read2:".$exp_read2{ $_ }."\n";
    }
    
    close( META );
    
    symlink( $pdirName.'/meta.info', $mdirName.'/meta.info' );
  }
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
		my $samples_csv = $ENV{ 'G1K' }.'/ref/G1K_samples.txt';
		
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

sub buildInternalHierarchy 
{
	my $projecthash = shift; # ref to hash of projectnames->samplenames->lists of fastq
	my $dhierarchyDir = shift;
	my $ahierarchyDir = shift;
	
	croak "Can't find data hierarchy directory\n" unless -d $dhierarchyDir;
	croak "Can't find analysis hierarchy directory\n" unless -d $ahierarchyDir;
	
	foreach my $project (keys %$projecthash)
	{
		print "Updating project: $project\n";
		
		( my $projectDirName = $project ) =~ s/\W+/_/g;
		
		my $projPath = $dhierarchyDir.'/'.$projectDirName;
		if( ! -d $projPath )
		{
			mkdir $projPath or die "Cannot create directory $projPath\n";
		}
		
		my $aprojPath = $ahierarchyDir.'/'.$projectDirName;
		if( ! -d $aprojPath )
		{
			mkdir $aprojPath or die "Cannot create directory $aprojPath\n";
		}
		
		my $numLibraries = 0;
		
		foreach my $library (keys %{$projecthash->{$project}})
		{
			print "Updating library: $library\n";
			
			#hack for G1K where sample starts with the individual
			my $individual = 1;
			if( $library =~ /^(NA\d+).*/ )
			{
				$individual = $1;
			}
			
			( my $libraryHierarchyName = $library ) =~ s/\W+/_/g;
			
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
			$path .= '/'.$libraryHierarchyName;
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
			$apath .= '/'.$libraryHierarchyName;
			
			if( ! -d $apath )
			{
				mkdir $apath or die "Cannot create directory $apath\n";
			}
			
			my @dirs;
			foreach my $fastq (@{$projecthash->{$project}{$library}})
			{
				$fastq = basename( $fastq );
				
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
        
        unless (getMpsaFastq($fastq)){
            print "Error retrieving $fastq with mpsa_download\n";
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
        my $cmd = qq[bsub -J split.$random -o import.o -e import.e -q $LSF_QUEUE perl -w -e "use AssemblyTools;AssemblyTools::sanger2SplitFastq( '$fastq', '$fastq1', '$fastq2');unlink '$fastq';"];
        system( $cmd );
        
        #work out the clip points (if any)
		$cmd = qq[bsub -J clip.$random.1 -w "done(split.$random)" -q small -o import.o -e import.e perl -w -e "use AssemblyTools;AssemblyTools::writeClipPointMeta( '$fastq1', 'meta.info', '1' );"];
		system( $cmd );

		$cmd = qq[bsub -J clip.$random.2 -w "done(split.$random)" -q small -o import.o -e import.e perl -w -e "use AssemblyTools;AssemblyTools::writeClipPointMeta( '$fastq2', 'meta.info', '2' );"];
        system( $cmd );
		
        #run fastqcheck
        $cmd = qq[bsub -J fastqcheck.$random.1 -o import.o -e import.e -q $LSF_QUEUE -w "done(split.$random)" "cat $fastq1 | fastqcheck > $fastq1.gz.fastqcheck;ln -fs $lPath/$fastq1.gz.fastqcheck $alPath/$fastq1.gz.fastqcheck"];
        system( $cmd );
        
        $cmd = qq[bsub -J fastqcheck.$random.2 -o import.o -e import.e -q $LSF_QUEUE -w "done(split.$random)" "cat $fastq2 | fastqcheck> $fastq2.gz.fastqcheck;ln -fs $lPath/$fastq2.gz.fastqcheck $alPath/$fastq2.gz.fastqcheck"];
        system( $cmd );
        
        #gzip the split fastq files
        $cmd = qq[bsub -q $LSF_QUEUE -w "done(clip.$random.*)&&done(fastqcheck.$random.*)" -o import.o -e import.e "gzip $fastq1 $fastq2; ln -fs $lPath/$fastq1.gz $alPath/$fastq1.gz;ln -fs $lPath/$fastq2.gz $alPath/$fastq2.gz"];
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
      
    unless (getMpsaFastq($fastq)){
        print "Error retrieving $fastq with mpsa_download\n";
        next;
    }
      
      #run fastqcheck
      my $cmd = qq[bsub -J fastqcheck.$random -o import.o -e import.e -q $LSF_QUEUE "cat $fastq | fastqcheck > $fastq.gz.fastqcheck; ln -fs $lPath/$fastq.gz.fastqcheck $alPath/$fastq.gz.fastqcheck"];
      system( $cmd );
      
      print "Writing meta.info for: $fastq.gz\n";
      
      #then write the meta info file
      open( META, ">>$lPath/meta.info" ) or die "Cannot create meta file in $lPath\n";
      
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
      $cmd = qq[bsub -J clip.$random -q small -o import.o -e import.e perl -w -e "use AssemblyTools;AssemblyTools::writeClipPointMeta( '$fastq', '$lPath/meta.info', '$readNum' );"];
	  
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

=head2 buildInternalHierarchyRecalibrated

    Arg [1]    : lane directory
    Returntype : 2d array reference where each array consists of [ fastq name ][ read length ][ num reads ][ num bases ]
                where this information is derived from the meta.info file and fastqcheck files
=cut
sub buildInternalHierarchyRecalibrated
{
    croak "Usage: buildInternalHierarchyRecalibrated fastq_files_fofn lane_dir_fofn new_root_dir copy|move" unless @_ == 4;
    
    my $recal_fofn = shift;
    my $old_lanes_fofn = shift;
    my $new_root = shift;
    my $mode = shift;
    
    croak "Cant find recal fofn: $recal_fofn\n" unless -f $recal_fofn;
    croak "Cant find old lanes fofn: $old_lanes_fofn\n" unless -f $old_lanes_fofn;
    croak "Cant find new root dir: $new_root\n" unless -d $new_root;
    
    croak "must specify copy or move as mode: $mode\n" unless $mode eq 'copy' || $mode eq 'move';
    
    mkdir( $new_root."/DATA" ) unless -d $new_root."/DATA";
    mkdir( $new_root."/MAPPING" ) unless -d $new_root."/MAPPING";
    my $dataRoot = $new_root."/DATA";
    my $mapRoot = $new_root."/MAPPING";
    
    my %lane_dirs;
    open( my $ldFh, $old_lanes_fofn ) or die $!;
    while( <$ldFh> )
    {
        chomp;
        
        $lane_dirs{ basename( $_ ) } = $_;
    }
    close( $ldFh );
    
    my %linked;
    open( my $recalFh, $recal_fofn ) or print $!;
    while( <$recalFh> )
    {
        chomp;
        my $filename = basename( $_ );
        
        croak "No fastq name found in fastq pathname: $_\n" unless $_ =~ /.*fastq.*/;
        
        next unless ! defined( $linked{ $filename } );
        
        $filename =~ /(\d+)_(\d+)[_\d+]*([\.recal]*\.fastq\.gz)/;
        
        croak "Malformed fastq filename: $filename\n" unless defined $1 && defined $2;
        
        my $laneDir = "$1_$2";
        
        if( defined( $lane_dirs{ $laneDir } ) )
        {
            #make the lane directory in the new root
            if( ! -d "$dataRoot/".$lane_dirs{ $laneDir } )
            {
                 system( "mkdir -p $dataRoot/".$lane_dirs{ $laneDir } ) == 0 or die "Failed to create data direcotry: $?\n";
                 system( "mkdir -p $mapRoot/".$lane_dirs{ $laneDir } ) == 0 or die "Failed to create mapping direcotry: $?\n";
            }
            
            chdir( "$dataRoot/".$lane_dirs{ $laneDir } );
            
            if( -s $_ )
            {
                print "Already imported fastq: $_\n";
                next;
            }
            
            if( $mode eq 'copy' )
            {
                #copy in the fastq to data
                copy( $_, "." ) or die "Failed to copy fastq file: $_ $!\n";
            }
            else
            {
                move( $_, "." ) or die "Failed to move fastq file: $_ $!\n";
            }
            
            #link the fastq into mapping
            symlink( "$dataRoot/".$lane_dirs{ $laneDir }."/".$filename, $mapRoot."/".$lane_dirs{ $laneDir }."/".$filename ) or die "failed to sym link fastq file: $!\n";
            
            #setup the reference links
            
            #run a fastqcheck
            my $cmd = qq[bsub -q small -o import.o -e import.e "zcat $filename | fastqcheck > $filename.fastqcheck; ln -s $dataRoot/$lane_dirs{ $laneDir }/$filename.fastqcheck $mapRoot/$lane_dirs{ $laneDir }"];
            system( $cmd );
            
            #setup the meta info file
            if( $filename =~ /^\d+_\d+[\.recal]*\.fastq\.gz$/ ) #if unpaired
            {
                system( qq[echo "read0:$filename" > meta.info] ) == 0 or die "Failed to write to meta.info file:$!\n";
                symlink( qq[$dataRoot/$lane_dirs{ $laneDir }/meta.info], "$mapRoot/$lane_dirs{ $laneDir }/meta.info" ) or die "Failed to sym link meta.info: $!";
                
                #run the clip point script
                $cmd = qq[bsub -q small -o import.o -e import.e perl -w -e "use AssemblyTools;AssemblyTools::writeClipPointMeta( '$filename', 'meta.info', '0' );"];
                system( $cmd );
            }
            else
            {
                if( $filename =~ /^\d+_\d+_2[\.recal]*\.fastq\.gz$/ ) #read 2
                {
                    #2nd read in pair
                    system( qq[echo "read2:$filename" >> meta.info] ) == 0 or die "Failed to write to meta.info file:$!\n";
                    symlink( qq[$dataRoot/$lane_dirs{ $laneDir }/meta.info], "$mapRoot/$lane_dirs{ $laneDir }/meta.info" ) or die "Failed to sym link meta.info: $!" unless -l "$mapRoot/$lane_dirs{ $laneDir }/meta.info";
                    
                    #run the clip point script
                    $cmd = qq[bsub -q small -o import.o -e import.e perl -w -e "use AssemblyTools;AssemblyTools::writeClipPointMeta( '$filename', 'meta.info', '2' );"];
                    system( $cmd );
                }
                elsif( $filename =~ /^\d+_\d+_1[\.recal]*\.fastq\.gz$/ ) #read 1
                {
                    system( qq[echo "read1:$filename" >> meta.info] ) == 0 or die "Failed to write to meta.info file:$!\n";
                    symlink( qq[$dataRoot/$lane_dirs{ $laneDir }/meta.info], "$mapRoot/$lane_dirs{ $laneDir }/meta.info" ) or die "Failed to sym link meta.info: $!" unless -l "$mapRoot/$lane_dirs{ $laneDir }/meta.info";
                    
                    #run the clip point script
                    $cmd = qq[bsub -q small -o import.o -e import.e perl -w -e "use AssemblyTools;AssemblyTools::writeClipPointMeta( '$filename', 'meta.info', '1' );"];
                    system( $cmd );
                }
            }
        }
        
        $linked{ $filename } = 1;
    }
    close( $recalFh );
}

1;
