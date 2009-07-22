package SamTools;
use strict;
use Carp;
use File::Spec;
use File::Basename;
use Cwd;
use Digest::MD5  qw(md5 md5_hex md5_base64);
use Utils;

use AssemblyTools;

=pod

=head1 DATA STRUCTURES

$FLAGS

=cut

our $FLAGS = 
{
    'paired_tech'    => 0x0001,
    'read_mapped'    => 0x0002,
    'unmapped'       => 0x0004,
    'mate_unmapped'  => 0x0008,
    'reverse_strand' => 0x0010,
    'mate_reverse'   => 0x0020,
    '1st_in_pair'    => 0x0040,
    '2nd_in_pair'    => 0x0080,
    'not_primary'    => 0x0100,
    'failed_qc'      => 0x0200,
    'duplicate'      => 0x0400,
};


=head1 METHODS

=head2 ssaha2samUnpaired

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
    $read_group = $read_group ? "\tRG:Z:".$read_group : '';
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
    my %cigar1Hash = %{ $ref };
    
    my $numReadsWritten = 0;
    my $numReadsDiscarded = 0;
    my $totalReadsHit = 0;
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
            $totalReadsHit ++;
            my $samCigar1 = ssahaCigar2SamCigar( $hits1[ 0 ], $seq1 );
            
            my $mapScore = detemineMappingScore( $cigar1Hash{ $name1 } );
            
            if( $mapScore > -1 )
            {
                my @s1 = split( /\s+/, $hits1[ 0 ] );
                
                my $flag = determineSamFlagUnpaired( $hits1[ 0 ] );
                
                $seq1 = AssemblyTools::revCompDNA( $seq1 ) if( $s1[ 4 ] eq '-' );
                $quals1 = reverse( $quals1 ) if( $s1[ 4 ] eq '-' );
                
                print SAM "$name1\t".$flag."\t".$s1[ 5 ]."\t".$s1[ 6 ]."\t".$mapScore."\t".$samCigar1."\t*\t0\t0\t".$seq1."\t".$quals1.$read_group."\n";
                $numReadsWritten ++;
            }
            else
            {
                $numReadsDiscarded ++;
            }
        }
    }
    close( READS1 );
    close( SAM );
    
    print "num_reads_written:$numReadsWritten\n";
    print "num_reads_discarded:$numReadsDiscarded\n";
    print "num_reads_cigar_hit:$totalReadsHit\n";
}

=head2 ssaha2samPaired

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
    $read_group = $read_group ? "\tRG:Z:".$read_group : '';
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
    my %cigar1Hash = %{ $ref };
    $ref = AssemblyTools::hash_ssaha_cigar_output( $cigar2 );
    my %cigar2Hash = %{ $ref };
    
    my $numReadsWritten = 0;
    my $numReadsDiscarded = 0;
    my $totalReadsHit = 0;
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
            $totalReadsHit += 2;
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
            my $mapScore2 = detemineMappingScore( \@hits2 );
            if( $flags[ 2 ] == 1 && $mapScore1 > -1 && $mapScore2 > -1 ) #if paired consistently - then double the mapping scores (RD suggest)
            {
                $mapScore1 = $mapScore1 * 2;
                $mapScore2 = $mapScore2 * 2;
                
                $mapScore1 = 254 if( $mapScore1 > 254 );
                $mapScore2 = 254 if( $mapScore2 > 254 );
            }
            
            $seq1 = AssemblyTools::revCompDNA( $seq1 ) if( $s1[ 4 ] eq '-' );
            $quals1 = reverse( $quals1 ) if( $s1[ 4 ] eq '-' );
            $seq2 = AssemblyTools::revCompDNA( $seq2 ) if( $s2[ 4 ] eq '-' );
            $quals2 = reverse( $quals2 ) if( $s2[ 4 ] eq '-' );
            
            if( $mapScore1 > -1 && $mapScore2 > -1 )
            {
                print SAM "$name1\t".$flags[ 0 ]."\t".$s1[ 5 ]."\t".$s1[ 6 ]."\t".$mapScore1."\t".$samCigar1."\t".($s2[ 5 ] eq $s1[ 5 ] ? "=" : $s2[ 5 ])."\t".$s2[ 6 ]."\t".( determineInsertSize( $s1[ 6], $s1[ 7 ], $s2[ 6 ], $s2[ 7 ] ) )."\t".$seq1."\t".$quals1.$read_group."\n";
                $numReadsWritten ++;
                print SAM "$name2\t".$flags[ 1 ]."\t".$s2[ 5 ]."\t".$s2[ 6 ]."\t".$mapScore2."\t".$samCigar2."\t".($s1[ 5 ] eq $s2[ 5 ] ? "=" : $s1[ 5 ])."\t".$s1[ 6 ]."\t".( determineInsertSize( $s1[ 6 ], $s1[ 7 ], $s2[ 6 ], $s2[ 7 ] ) )."\t".$seq2."\t".$quals2.$read_group."\n";
                $numReadsWritten ++;
                next;
            }
            elsif( $mapScore1 > -1 ) #just one read aligns
            {
                $flags[ 0 ] += hex( "0x0008" );
                print SAM "$name1\t".$flags[ 0 ]."\t".$s1[ 5 ]."\t".$s1[ 6 ]."\t".$mapScore1."\t".$samCigar1."\t*\t0\t0\t".$seq1."\t".$quals1.$read_group."\n";
                $numReadsWritten ++;
                $numReadsDiscarded ++;
                next;
            }
            elsif( $mapScore2 > -1 )
            {
                $flags[ 1 ] += hex( "0x0008" );
                print SAM "$name2\t".$flags[ 1 ]."\t".$s2[ 5 ]."\t".$s2[ 6 ]."\t".$mapScore2."\t".$samCigar2."\t*\t0\t0\t".$seq2."\t".$quals2.$read_group."\n";
                $numReadsWritten ++;
                $numReadsDiscarded ++;
                next;
            }
            else
            {
                $numReadsDiscarded += 2;
            }
        }
        elsif( @hits1 > 0 )
        {
            $totalReadsHit ++;
            my $samCigar1 = ssahaCigar2SamCigar( $hits1[ 0 ], $seq1 );
            my @s1 = split( /\s+/, $hits1[ 0 ] );
            
            my $flag = determineSamFlagUnpaired( $hits1[ 0 ] );
            
            my $mapScore = detemineMappingScore( \@hits1 );
            
            if( $mapScore > -1 )
            {
                $seq1 = AssemblyTools::revCompDNA( $seq1 ) if( $s1[ 4 ] eq '-' );
                $quals1 = reverse( $quals1 ) if( $s1[ 4 ] eq '-' );
                
                print SAM "$name1\t".$flag."\t".$s1[ 5 ]."\t".$s1[ 6 ]."\t".$mapScore."\t".$samCigar1."\t*\t0\t0\t".$seq1."\t".$quals1.$read_group."\n";
                $numReadsWritten ++;
            }
            else
            {
                $numReadsDiscarded ++;
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
            $totalReadsHit ++;
            my $samCigar2 = ssahaCigar2SamCigar( $hits2[ 0 ], $seq2 );
            my @s2 = split( /\s+/, $hits2[ 0 ] );
            
            my $flag = determineSamFlagUnpaired( $hits2[ 0 ] );
            
            my $mapScore = detemineMappingScore( \@hits2 );
            
            if( $mapScore > -1 )
            {
                $seq2 = AssemblyTools::revCompDNA( $seq2 ) if( $s2[ 4 ] eq '-' );
                $quals2 = reverse( $quals2 ) if( $s2[ 4 ] eq '-' );
                
                print SAM "$name2\t".$flag."\t".$s2[ 5 ]."\t".$s2[ 6 ]."\t".$mapScore."\t".$samCigar2."\t*\t0\t0\t".$seq2."\t".$quals2.$read_group."\n";
                $numReadsWritten ++;
            }
            else
            {
                $numReadsDiscarded ++;
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
    
    print "num_reads_written:$numReadsWritten\n";
    print "num_reads_discarded:$numReadsDiscarded\n";
    print "num_reads_cigar_hit:$totalReadsHit\n";
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

=head2 determineHits

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

=head2 detemineMappingScore

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
        return 254;
    }
    else
    {
        my @hit2 = split( /\s+/, $hits[ 1 ] );
        
        my $d = ( $hit1[ 9 ] - $hit2[ 9 ] ) - ( 2 * scalar( @hits ) );
        
        return $d < 255 ? $d : 254;
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
        if( $s2[ 4 ] eq '-' )
        {
            $flag1 += hex( '0x0020' );
        }
        
        $flag2 += hex( '0x0010' ) if( $s2[ 4 ] eq '-' );
        if( $s1[ 4 ] eq '-' )
        {
            $flag2 += hex( '0x0020' );
        }
        
        $flag1 += hex( '0x0040' );
        $flag2 += hex( '0x0080' );
        
        $consistent = 1;
    }
    elsif( length( $cigar1 ) > 0 )
    {
        $flag1 += hex( '0x0008' );
        $flag1 += hex( '0x0010' ) if( $s1[ 4 ] eq '-' );
    }
    elsif( length( $cigar2 ) > 0 )
    {
        $flag2 += hex( '0x0008' );
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

sub bam_stat_fofn
{
    croak "Usage: bam_stat bam_fofn\n" unless @_ == 1;
    my $bam_fofn = shift;
    
    croak "Cant find bam_fofn file: $bam_fofn\n" unless -f $bam_fofn;
    
    open( F, $bam_fofn ) or die "cant open bam fofn: $bam_fofn\n";
    while( <F> )
    {
        chomp;
        if( -f $_ )
        {
            my $statF = $_."bamstat";
            my $cmd = qq[bsub -q normal -o $$.bam.o -e $$.bam.e perl -w -e "use SamTools;SamTools::bam_stat( '$_', '$statF');"];
            #print $cmd."\n";
            system( $cmd );
        }
        else
        {
            print "Cant find BAM file: $_\n";
        }
    }
    close( F );
}

sub bam_stat
{
    croak "Usage: bam_stat bam output_file\n" unless @_ == 2;
    my $bam = shift;
    my $output = shift;
    
    my $numReads = 0;
    my $numBases = 0;
    my $insert_mean = 0;
    open( BAM, "samtools view $bam |" ) or die "Cannot open bam file: $bam\n";
    while( <BAM> )
    {
        chomp;
        my @s = split( /\s+/, $_ );
        # only consider reads that actually mapped
        next if ($s[1] & $FLAGS->{unmapped});
        $numBases += length( $s[ 9 ] );
        my $ins = abs( $s[ 8 ] );
        $insert_mean = ( ( $insert_mean * $numReads ) + $ins ) / ( $numReads + 1 ) unless $ins > 10000 && $ins != 0;
        $numReads ++;
    }
    close( BAM );
    
    #do a second pass to get the std deviations (uses too much memory if store all inserts in memory)
    my $sumSqVariances = 0;
    open( BAM, "samtools view $bam |" ) or die "Cannot open bam file: $bam\n";
    while( <BAM> )
    {
        chomp;
        my @s = split( /\s+/, $_ );
        my $ins = abs( $s[ 8 ] );
        $sumSqVariances += ( ( $insert_mean - $ins ) * ( $insert_mean - $ins ) ) unless $ins > 10000 && $ins != 0;
    }
    close( BAM );
    
    my $std = $sumSqVariances ? sqrt( ( $sumSqVariances / $numReads ) ) : 0;
    
    open( O, ">$output" ) or die "Cannot create output: $!\n";
    print O "num_bases:$numBases\n";
    print O "num_reads:$numReads\n";
    print O "avg_insert:$insert_mean\n";
    print O "std_dev:$std\n";
    close( O );
}

sub parse_bam_header_line
{
    my ($line) = @_;
    my $out = {};

    # Header line:
    #   @RG     ID:ERR001773    PL:ILLUMINA     PU:IL23_337_6   LB:g1k-sc-NA12878-CEU-2 PI:200  SM:NA12878      CN:SC
    my @items = split /\t/, $line;
    shift(@items);
    for my $pair (@items)
    {
        my ($key,$value) = split /:/,$pair;
        $$out{$key} = $value;
    }
    return $out;
}

sub parse_bam_line
{
    my ($line) = @_;
    my $out = {};

    #   IL14_1902:3:83:1158:1446        89      1       23069154        37      54M     =       23069154        0       TGCAC ... RG:Z:ERR001720
    my @items = split /\t/, $line;

    $$out{'flag'}  = $items[1];
    $$out{'chrom'} = $items[2];
    $$out{'pos'}   = $items[3];
    $$out{'cigar'} = $items[5];
    $$out{'isize'} = $items[8];
    $$out{'seq'}   = $items[9];

    my $nitems = @items;
    for (my $i=11; $i<$nitems; $i++)
    {
        my ($key,@vals) = split /:/,$items[$i];

        # e.g. RG:Z:ERR001720
        if ( $key eq 'RG' )
        {
            $$out{'RG'} = $vals[1];
        }

        # e.g. NM:i:0
        elsif( $key eq 'NM' )
        {
            $$out{'NM'} = $vals[1];
        }
    }
    return $out;
}


=head2 collect_detailed_bam_stats

        Description: Reads the given bam file and runs detailed statistics, such as 
                        histograms of insert size; QC content of mapped and both mapped
                        and unmapped sequences; reads distribution with respect to
                        chromosomes; duplication rate.
        Arg [1]    : The .bam file (sorted and indexed).
        Arg [2]    : The .fai file (to determine chromosome lengths). (Can be NULL with do_chrm=>0.)
        Arg [3]    : Options (hash) [Optional]
        Options    : (see the code for default values)
                        do_chrm          .. should collect the chromosome distrib. stats.
                        do_gc_content    .. should we collect the gc_content (default is 1)
                        do_rmdup         .. default is 1 (calculate the rmdup); alternatively supply the filename of a pre-calculated rmdup

                        insert_size_bin  .. the length of the distribution intervals for the insert size frequencies
                        gc_content_bin   

        Returntype : Hash with the following entries, each statistics type is a hash as well. See also Graphs::plot_stats.
                        insert_size =>
                            yvals      .. array of insert size frequencies 
                            xvals      .. 
                            max => x   .. maximum values
                            max => y
                            average    .. average value (estimated from histogram, the binsize may influence the accuracy)
                        gc_content_forward  .. gc content of sequences with the 0x0040 flag set
                        gc_content_reverse  .. gc content of seqs with the 0x0080 flag set
                        reads_chrm_distrib  .. read distribution with respect to chromosomes
                        reads_total
                        reads_paired
                        reads_mapped        .. paired + unpaired
                        reads_unpaired
                        reads_unmapped
                        bases_total
                        bases_mapped
                        duplication
                        num_mismatches (if defined in bam file, otherwise 0)
=cut

sub collect_detailed_bam_stats
{
    my ($bam_file,$fai_file,$options) = @_;
    if ( !$bam_file ) { Utils::error("Expected .bam file as a parameter.\n") }
    
    $options = {} unless $options;
    my $insert_size_bin = exists($$options{'insert_size_bin'}) ? $$options{'insert_size_bin'} : 1;
    my $gc_content_bin  = exists($$options{'gc_content_bin'}) ? $$options{'gc_content_bin'} : 1;

    my $do_chrm  = exists($$options{'do_chrm'}) ?  $$options{'do_chrm'} : 1;
    my $do_gc    = exists($$options{'do_gc_content'}) ? $$options{'do_gc_content'} : 1;
    my $do_rmdup = exists($$options{'do_rmdup'}) ? $$options{'do_rmdup'} : 1;

    my $chrm_lengths = $do_chrm ? Utils::fai_chromosome_lengths($fai_file) : {};

    # Use hashes, not arrays - the data can be broken and we might end up allocating insanely big arrays 

    #   reads_unmapped  ..   That is, not aligned to the ref sequence (unmapped flag)
    #   reads_paired    ..      both mates are mapped and form a pair (read_mapped flag)
    #   reads_unpaired  ..      both mates mapped, but do not form a pair (neither unmapped nor read_mapped flag set)
    #   bases_total     ..   The total number of bases determined as \sum_{seq} length(seq)
    #   bases_mapped    ..      number of bases with 'M' in cigar

    my $raw_stats = {};

    # Collect the statistics - always collect the total statistics for all lanes ('total') and
    #   when the RG information is present in the header, collect also statistics individually
    #   for each ID (@RG ID:xyz). The statistics are collected into raw_stats hash, individual
    #   keys are 'total' and IDs.
    #
    my $i=0;
    open(my $fh, "samtools view -h $bam_file |") or Utils::error("samtools view -h $bam_file |: $!");
    while (my $line=<$fh>)
    {
        # Header line:
        #   @RG     ID:ERR001773    PL:ILLUMINA     PU:IL23_337_6   LB:g1k-sc-NA12878-CEU-2 PI:200  SM:NA12878      CN:SC
        #
        # Data line:
        #   IL14_1902:3:83:1158:1446        89      1       23069154        37      54M     =       23069154        0       TGCAC ... RG:Z:ERR001720
        #
        if ( $line=~/^\@/ )
        {
            if ( !($line=~/^\@RG/) ) { next }

            my $header = parse_bam_header_line($line);
            if ( !exists($$header{'ID'}) ) { Utils::error("No ID in the header line? $line\n"); }

            $$raw_stats{$$header{'ID'}}{'header'} = $header;
            next;
        }
		
        # The @stats array is a convenient way how to reuse the same code for the total and individual
        #   statistics - we will always add the same numbers to 'total' and the ID.
        #
        my @stats = ('total');
        my $data  = parse_bam_line($line);
        if ( exists($$data{'RG'}) ) { push @stats, $$data{'RG'}; }

        my $flag  = $$data{'flag'};
        my $chrom = $$data{'chrom'};
        my $pos   = $$data{'pos'};
        my $cigar = $$data{'cigar'};
        my $isize = $$data{'isize'};
        my $seq   = $$data{'seq'};
        
        my $seq_len  = length($seq);
        my $mismatch = exists($$data{'NM'}) ? $$data{'NM'} : 0;
        for my $stat (@stats)
        {
            $$raw_stats{$stat}{'reads_total'}++;
            $$raw_stats{$stat}{'bases_total'} += $seq_len;
            $$raw_stats{$stat}{'num_mismatches'} += $mismatch;
        }

        my $paired = ($flag & $$FLAGS{'read_mapped'}) && ($flag & $$FLAGS{'paired_tech'});
        if ( $paired || !($flag & $$FLAGS{'paired_tech'}) )
        {
            if ( $paired ) 
            { 
                for my $stat (@stats) { $$raw_stats{$stat}{'reads_paired'}++; }

                # Insert Size Frequencies
                #
                my $bin = abs(int( $isize / $insert_size_bin ));
                for my $stat (@stats) { $$raw_stats{$stat}{'insert_size_freqs'}{$bin}++; }
            }

            my $cigar_info = cigar_stats($cigar);
            if ( exists($$cigar_info{'M'}) )
            {
                for my $stat (@stats) { $$raw_stats{$stat}{'bases_mapped'} += $$cigar_info{'M'}; }
            }

            # Chromosome distribution
            #
            if ( $do_chrm && $chrom=~/^(?:\d+|X|Y)$/i ) 
            {
                for my $stat (@stats) { $$raw_stats{$stat}{'chrm_distrib_freqs'}{$chrom}++; }
            }
        }
        elsif ( $flag & $$FLAGS{'unmapped'} ) 
        { 
            for my $stat (@stats) { $$raw_stats{$stat}{'reads_unmapped'}++;  }
        }
        else 
        { 
            for my $stat (@stats) { $$raw_stats{$stat}{'reads_unpaired'}++; }

            my $cigar_info = cigar_stats($cigar);
            if ( exists($$cigar_info{'M'}) )
            {
                for my $stat (@stats) { $$raw_stats{$stat}{'bases_mapped'} += $$cigar_info{'M'}; }
            }
        }

        # GC Content Frequencies - collect stats for both pairs separately
        if ( $do_gc )
        {
            my $gc_count = 0;
            for (my $ipos=0; $ipos<$seq_len; $ipos++)
            {
                my $nuc  = substr($seq, $ipos, 1);
                if ( $nuc eq 'g' || $nuc eq 'G' || $nuc eq 'c' || $nuc eq 'C' ) { $gc_count++; }
            }
            $gc_count = $gc_count*100./$seq_len;
            my $bin = abs(int( $gc_count / $gc_content_bin ));
            if ( $flag & $$FLAGS{'1st_in_pair'} )
            {
                for my $stat (@stats) { $$raw_stats{$stat}{'gc_content_fwd_freqs'}{$bin}++; }
            }
            elsif ( $flag & $$FLAGS{'2nd_in_pair'} )
            {
                for my $stat (@stats) { $$raw_stats{$stat}{'gc_content_rev_freqs'}{$bin}++; }
            }
            elsif ( !($flag & $$FLAGS{'paired_tech'}) )  # Not a paired-read technology
            { 
                # Either it is a non-paired-read technology, or the 1st_in_pair and
                #   and 2nd_in_pair flags got lost in the process. (The specs allows this.)
                for my $stat (@stats) { $$raw_stats{$stat}{'gc_content_fwd_freqs'}{$bin}++; }
            }
        }

        #if ( $i++>50000 ) { last }
    }
    close $fh;

    # This sanity check could be used for paired reads only.
    #
    #   if ( ($flag & $$FLAGS{'paired_tech'}) && $reads_total != $reads_paired + $reads_unmapped + $reads_unpaired )
    #   {
    #       Utils::error("FIXME: Incorrect assumption: paired + unpaired + unmapped != total ($reads_paired+$reads_unpaired+$reads_unmapped!=$reads_total)\n");
    #   }
    #
    # This calculation worked only for paired reads.
    #   my $reads_mapped = $reads_unpaired + $reads_paired;
    #
    for my $stat (keys %$raw_stats) 
    { 
        $$raw_stats{$stat}{'reads_unmapped'} = 0 unless exists($$raw_stats{$stat}{'reads_unmapped'});
        $$raw_stats{$stat}{'reads_unpaired'} = 0 unless exists($$raw_stats{$stat}{'reads_unpaired'});
        $$raw_stats{$stat}{'reads_paired'}   = 0 unless exists($$raw_stats{$stat}{'reads_paired'});
        $$raw_stats{$stat}{'reads_total'}    = 0 unless exists($$raw_stats{$stat}{'reads_total'});

        $$raw_stats{$stat}{'reads_mapped'} = $$raw_stats{$stat}{'reads_total'} - $$raw_stats{$stat}{'reads_unmapped'}; 
    }

    # Find out the duplication rate
    if ( $do_rmdup )
    {
    my $rmdup_reads_total;
    if (-f $do_rmdup && -s $do_rmdup) {
        chomp(($rmdup_reads_total) = Utils::CMD("wc -l $do_rmdup"));
    }
    else {
        chomp(($rmdup_reads_total) = Utils::CMD("samtools rmdup $bam_file - 2>/dev/null | samtools view - | wc -l"));
    }
        $$raw_stats{'total'}{'rmdup_reads_total'} = $rmdup_reads_total;
    }

    # Now process the reults. The raw_stats hash now contains the total statistics (the key 'total')
    #   and possibly also separate statistics for individual libraries (other keys, e.g. 'ERR001773').
    #
    my $stats = {};
    for my $stat (keys %$raw_stats)
    {
        $$stats{$stat} = 
        {
            'reads_total'       => $$raw_stats{$stat}{'reads_total'},
            'reads_paired'      => $$raw_stats{$stat}{'reads_paired'},
            'reads_mapped'      => $$raw_stats{$stat}{'reads_mapped'},
            'reads_unpaired'    => $$raw_stats{$stat}{'reads_unpaired'},
            'reads_unmapped'    => $$raw_stats{$stat}{'reads_unmapped'},
            'bases_total'       => $$raw_stats{$stat}{'bases_total'},
            'bases_mapped'      => $$raw_stats{$stat}{'bases_mapped'},
            'num_mismatches'    => $$raw_stats{$stat}{'num_mismatches'},

            'insert_size' => 
            {
                'data'       => $$raw_stats{$stat}{'insert_size_freqs'},
                'bin_size'   => $insert_size_bin,
            },
        };
        if ( exists($$raw_stats{$stat}{'header'}) )
        {
            $$stats{$stat}{'header'} = $$raw_stats{$stat}{'header'};
        }
        if ( exists($$raw_stats{$stat}{'chrm_distrib_freqs'}) )
        {
            $$stats{$stat}{'chrm_distrib_freqs'} = $$raw_stats{$stat}{'chrm_distrib_freqs'};
        }
        if ( exists($$raw_stats{$stat}{'rmdup_reads_total'}) )
        {
            $$stats{$stat}{'rmdup_reads_total'} = $$raw_stats{$stat}{'rmdup_reads_total'};
            $$stats{$stat}{'duplication'}       = $$raw_stats{$stat}{'rmdup_reads_total'}/$$raw_stats{$stat}{'reads_total'};
        }
        if ( exists($$raw_stats{$stat}{'gc_content_fwd_freqs'}) )
        {
            $$stats{$stat}{'gc_content_forward'} = 
            {
                'data'       => $$raw_stats{$stat}{'gc_content_fwd_freqs'},
                'bin_size'   => $gc_content_bin,
            };
        }
        if ( exists($$raw_stats{$stat}{'gc_content_rev_freqs'}) )
        {
            $$stats{$stat}{'gc_content_reverse'} =
            {
                'data'       => $$raw_stats{$stat}{'gc_content_rev_freqs'},
                'bin_size'   => $gc_content_bin,
            };
        }
    }

    # Convert the hashes into arrays and find extreme values - this is for convenience only,
    #   TrackQC uses this. Although the code is lengthy, the histograms are small and take
    #   no time to process.
    #
    for my $stat_name (keys %$stats)
    {
        my $stat = $$stats{$stat_name};

        for my $key (keys %$stat)
        {
            if ( ref($$stat{$key}) ne 'HASH' ) { next }
            if ( !exists($$stat{$key}{'data'}) ) { next }
            if ( !exists($$stat{$key}{'bin_size'}) ) { next }

            my @yvals = ();
            my @xvals = ();
            my ($ymax,$xmax);

            my $avg  = 0;
            my $navg = 0;

            my $data = $$stat{$key}{'data'};
            for my $ibin (sort {$a<=>$b} keys %$data)
            {
                my $bin = $ibin * $$stat{$key}{'bin_size'};
                if ( $$stat{$key}{'bin_size'}>1 ) { $bin += $$stat{$key}{'bin_size'}*0.5; }
                push @xvals, $bin;

                my $yval = $$data{$ibin};
                push @yvals, $yval;

                $avg  += $yval*$bin;    # yval is the count and bin the value
                $navg += $yval;

                if ( !defined $ymax || $yval>$ymax ) { $xmax=$bin; $ymax=$yval }
            }
            $$stat{$key}{'xvals'} = scalar @xvals ? \@xvals : [0];  # Yes, this can happen,e.g. AKR_J_SLX_200_NOPCR_1/1902_3
            $$stat{$key}{'yvals'} = scalar @yvals ? \@yvals : [0];
            $$stat{$key}{'max'}{'x'} = defined $xmax ? $xmax : 0;
            $$stat{$key}{'max'}{'y'} = defined $ymax ? $ymax : 0;
            $$stat{$key}{'average'}  = $avg/$navg;
        }

        # Chromosome distribution histograms (reads_chrm_distrib) are different - xvalues are not numeric
        #
        if ( exists($$stat{'chrm_distrib_freqs'}) )
        {
            my @yvals = ();
            my @xvals = ();
            for my $key (sort Utils::cmp_mixed keys %{$$stat{'chrm_distrib_freqs'}})
            {
                if ( !exists($$chrm_lengths{$key}) ) { Utils::error("The chromosome \"$key\" not in $fai_file.\n") }
                push @yvals, $$stat{'chrm_distrib_freqs'}{$key} / $$chrm_lengths{$key};
                push @xvals, $key;
            }   
            $$stat{'reads_chrm_distrib'}{'yvals'} = scalar @yvals ? \@yvals : [0];
            $$stat{'reads_chrm_distrib'}{'xvals'} = scalar @xvals ? \@xvals : [0];
            $$stat{'reads_chrm_distrib'}{'scaled_dev'} = scaled_chrm_dev(\@xvals,\@yvals);
        }
    }
    
    return $stats;
}


# Calculates some sort of standard devitation, except that the
#   the difference from mean is scaled by mean (so that the formula
#   works for samples of different size) and the instead of mean
#   the maximum is used (we want to detect one oversampled chromosome): 
#       sqrt[(1/N)*sum_i((y_i-max)/max)**2]
#
sub scaled_chrm_dev
{
    my ($xvals,$yvals) = @_;
    my $max   = 0;
    my $ndata = 0;
    for (my $i=0; $i<scalar @$xvals; $i++)
    {
        # We do not know if both chromosome X and Y should be counted - we don't know the sex
        if ( !($$xvals[$i]=~/^\d+$/) ) { next } 

        if ( !defined($max) || $max<$$yvals[$i] ) { $max=$$yvals[$i]; }
        $ndata++;
    }
    if ( !$ndata ) { return 0 }

    my $stddev = 0;
    for (my $i=0; $i<scalar @$xvals; $i++)
    {
        if ( !($$xvals[$i]=~/^\d+$/) ) { next } 

        $stddev += (($$yvals[$i] - $max)/$max)**2;
    }
    return sqrt($stddev/$ndata);
}


=head2 cigar_stats

        Description: Parse the SAM cigar and return some info.
        Arg [1]    : The cigar (e.g. 36M, 5M1I30M, 21M15S, etc.)
        Returntype : The hash e.g. { 'M'=>35, 'I'=>1 } 

=cut

sub cigar_stats
{
    my ($cigar) = @_;

    my $stats = {};
    while ($cigar)
    {
        if ( !($cigar=~/^(\d+)(\D)/) ) { last }
        $$stats{$2} += $1;
        $cigar = $';
    }

    return $stats;
}


=head2 print_flags

        Description: For debugging purposes, prints flags for all sequences in human readable form.
        Arg [1]    : The .bam file (sorted and indexed).
        Returntype : None

=cut

sub print_flags
{
    my ($bam_file,$options) = @_;

    my $fh = \*STDIN;
    if ( $bam_file ) { open($fh, "samtools view $bam_file |") or Utils::error("samtools view $bam_file |: $!"); }
    while (my $line=<$fh>)
    {
        # IL14_1902:3:83:1158:1446        89      1       23069154        37      54M     =       23069154        0       TGCAC
        my @items = split /\t/, $line;
        my $flag  = $items[1];
        my $chrom = $items[2];
        my $pos   = $items[3];
        my $isize = $items[8];
        my $seq   = $items[9];

        print $line;
        print "\t insert_size=$isize flag=$flag\n\t --\n";
        print debug_flag($flag);
        print "\n";
    }
    close $fh;
    return;
}

sub debug_flag
{
    my ($flag) = @_;

    my $out = '';
    for my $key (sort keys %$FLAGS)
    {
        if ( $flag & $$FLAGS{$key} ) { $out .= "\t $key\n" }
    }
    return $out;
}


=head2 determineUnmappedFlag

    Arg [1]    : 1/0 flag to say whether the lane is a paired lane
    Arg [2]    : -1/0/1 flag of whether the read is the first read of a pair
    Example    : determineUnmappedFlag( 0, 1 );
    Description: works out the SAM flag value for an unmapped read
    Returntype : int
=cut
sub determineUnmappedFlag
{
    croak "Usage: determineUnmappedFlag paired(0/1) read1(-1/0/1)\n" unless @_ == 2;
    
    my $paired = shift;
    my $read1 = shift;
    
    croak "Paired flag must be 0 or 1\n" unless $paired == 0 || $paired == 1;
    
    croak "Read1 flag must be -1, 0 or 1\n" unless $read1 == 0 || $read1 == 1 || $read1 == -1;
    
    my $total = $paired == 1 ? hex( '0x0001' ) : 0;
    $total += hex( '0x0004' );
    $total += hex( '0x0008' );
    
    if( $paired == 1 )
    {
        #if the read is a fragment from an paired run - then it wont be either end of a pair (e.g. 454)
        if( $read1 != -1 )
        {
            $total += $read1 == 1 ? hex( '0x0040' ) : hex( '0x0080' );
        }
    }
    
    return $total;
}

=head2 pileup2Intervals

    Arg [1]    : samtools pileup file
    Arg [2]    : sam header file
    Arg [3]    : output intervals file
    Example    : pileup2Intervals
    Description: creates an intervals file for the broad recalibrator which excludes all positions +/- 20bp around all indels
    Returntype : none
=cut

sub pileup2Intervals
{
    croak "Usage: pileup2Intervals samtools_pileup_file sam_header_file output_file\n" unless @_ == 3;
    
    my $pileup = shift;
    my $sam_header = shift;
    my $output = shift;
    
    croak "Cant find pileup file: $pileup" unless -f $pileup;
    croak "Cant find sam header file: $sam_header" unless -f $sam_header;
    
    my $pfh;
    if( $pileup =~ /\.gz$/ )
    {
        open( $pfh, "gunzip -c $pileup |" ) or die "Cannot open pileup file: $!\n";
    }
    else
    {
        open( $pfh, $pileup ) or die "Cannot open pileup file: $!\n";
    }
    
    system( qq[cat $sam_header > $output] );
    
    open( my $out, ">>$output" ) or die "Cannot create output file: $!\n";
    my $currentStartPos = 1;
    my $c = 0;
    
    my $currentChr = 'none';
    while( <$pfh> )
    {
        chomp;
        if( $_ =~ /^.+\t\d+\t\*\t.*/ )
        {
            my @s = split( /\t/, $_ );
            
            if( $currentChr eq 'none' || $currentChr ne $s[ 0 ] )
            {
                $currentChr = $s[ 0 ];
                $currentStartPos = 1;
            }
            
            my $stop = $s[ 1 ] - 20;
            
            if( $stop > $currentStartPos )
            {
                print $out qq/$s[ 0 ]\t$currentStartPos\t$stop\t+\ttarget_$c\n/;
                $currentStartPos = $s[ 1 ] + 20;
                $c ++;
            }
        }
    }
    close( $pfh );
    close( $out );
}

=head2 makeBamStat

	Arg [1]    : bam file
	Arg [2]    : output file
	Example    : makeBamStat
	Description: creates a bam stat file in the following tab-delimited format
				Study	sample	platform	RG	mapped_bases	total_reads	mapped_reads	paired_reads	properly_paired_reads	num_mismatches	avg_mapped_bases	avg_insert	sd_insert
	Returntype : none
=cut
sub makeBamStat
{
	croak "Usage: makeBamStat bam_file fai_file seq_index output_file\n" unless @_ == 4;
	
	my $bam_file = shift;
	my $fai_file = shift;
	my $output = shift;
	
	my %opts = ( 'do_chrm' => 0, 'do_gc_content' => 0, 'do_rmdup' => 0 );
	#print "Collecting bam stats.....\n";
	#my $stats = collect_detailed_bam_stats( $bam_file, $fai_file, \%opts );
	
	#get the study id
	#my $study = `grep `;
	print "Printing bam stats per RG....\n"
	open( my $bfh, "samtools view -H $bam_file|" ) or die $!;
	my $laneMeta = {};
	while( <$bfh> )
	{
		chomp;
		if( $_ =~ /^\@ZG/ )
		{
			$_ =~ /\t\@RG:/;
			print "Reading Header: $_\n";
			my $ref = parse_bam_header_line( $_ );
			#$$laneMeta{ $$ref{ 'RG' } } = $ref;
		}
	}
	close( $bfh );
=pod	
	foreach( keys( %$laneMeta ) )
	{
		print "$_\n";
		
	}

	my $output = '';
	foreach( my $id ( keys( %$stats ) ) )
	{
		$output .= qq[];
	}
	
	open( my $fh, ">$output" ) or die "Cant create output $!";
	print $fh qq[];
	close( $fh );
=cut
}

1;
