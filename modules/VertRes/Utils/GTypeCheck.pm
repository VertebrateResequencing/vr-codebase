package VertRes::Utils::GTypeCheck;
use base qw(VertRes::Base);

use strict;
use warnings;
use Carp;
use VertRes::LSF;
use Utils;

sub new 
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    return $self;
}


sub check_genotype
{
    my ($self) = @_;

    if ( !exists($$self{'bam'}) ) { $self->throw("Expected the \"bam\" option.\n") }
    if ( !exists($$self{'fa_ref'}) ) { $self->throw("Expected the \"fa_ref\" option.\n") }
    if ( !exists($$self{'glf'}) ) { $self->throw("Expected the \"glf\" option.\n") }
    if ( !exists($$self{'snps'}) ) { $self->throw("Expected the \"snps\" option.\n") }
    if ( !exists($$self{'samtools'}) ) { $self->throw("Expected the \"samtools\" option.\n") }
    if ( !exists($$self{'min_glf_ratio'}) ) { $self->throw("Expected the \"min_glf_ratio\" option.\n") }
    if ( !exists($$self{'lock_file'}) ) { $self->throw("Expected the \"lock_file\" option.\n") }

    my $bam         = $$self{'bam'};
    my $fa_ref      = $$self{'fa_ref'};
    my $glf         = $$self{'glf'};
    my $snps        = $$self{'snps'};
    my $samtools    = $$self{'samtools'};
    my $genotype    = exists($$self{'genotype'}) ? $$self{'genotype'} : '';
    my $ratio       = $$self{'min_glf_ratio'};
    my $lock_file   = $$self{'lock_file'};
    my $prefix      = exists($$self{'prefix'}) ? $$self{'prefix'} : '_';

    # Splits e.g. '/path/to/file.bam' to '/path/to' and 'file'
    my ($dir,$name,$suff) = Utils::basename($bam);
    if ( !$dir ) { $dir = '.'; }

    my $pileup_out = "$name.bcf";
    my $pileup_cmd = "bin2hapmap -l $snps > $name.tmp.sites && $samtools mpileup -ugDI -d 1000 -l $name.tmp.sites -f $fa_ref $bam > $name.tmp.pileup && mv $name.tmp.pileup $pileup_out";
    my $checkGenotype_cmd = "$glf checkGenotype -s - $snps $pileup_out > $name.gtypey";

    if (exists $self->{'snp_sites'} ) {
        my $snp_sites_string;
        open my $f, $self->{'snp_sites'} or $self->throw("Error opening file ". $self->{'snp_sites'});

        while (<$f>) {
            chomp;
            my @a = split /\t/;
            $snp_sites_string .= "$a[0]:$a[1]-$a[1] ";
        }

        close $f;
        my $gtype_bam = "$name.gtype.tmp.sort";
        $pileup_out = "$name.bcf";
        $pileup_cmd = "$samtools view -bh $bam $snp_sites_string > $name.gtype.tmp.bam && $samtools sort $name.gtype.tmp.bam $gtype_bam && $samtools index $gtype_bam.bam && $samtools mpileup -ugDI -d 1000 -l " . $self->{'snp_sites'} . " -f $fa_ref $gtype_bam.bam > $name.gtype.tmp.pileup && mv $name.gtype.tmp.pileup $pileup_out && rm $name.gtype.tmp.bam";
        $checkGenotype_cmd = "$glf checkGenotype -s - $snps $pileup_out > $name.gtypey";
    }

    # Dynamic script to be run by LSF.
    open(my $fh, '>', "$dir/${prefix}genotype.pl") or $self->throw("$dir/${prefix}genotype.pl: $!");
    print $fh
    qq{
use Utils;
use VertRes::Utils::GTypeCheck;

\$base = VertRes::Base->new();

if ( ! -e "$bam.bai" || Utils::file_newer("$bam","$bam.bai") ) 
{ 
    Utils::CMD("$samtools index $bam"); 
}
if ( ! -e "$pileup_out" || Utils::file_newer("$bam","$pileup_out") )
{
    Utils::CMD(q[$pileup_cmd]);
}
if ( ! -e "$name.gtypex" || Utils::file_newer(q[$pileup_out],"$name.gtypex") )
{
    Utils::CMD(q[$checkGenotype_cmd]);
    if ( ! -s "$name.gtypey" ) 
    { 
        \$base->throw("FIXME: this should not happen, what's wrong with:\n\t$checkGenotype_cmd\n???\n");
    }
    rename("$name.gtypey", "$name.gtypex") or Utils::CMD("rename $name.gtypey $name.gtypex: \$!");
    if ( ! -s "$name.gtypex" ) 
    { 
        \$base->throw("FIXME: this should not happen, what's wrong with:\n\trename $name.gtypey $name.gtypex\n???\n");
    }
}
if ( -s "$name.gtypex" )
{
    my \$gtc = VertRes::Utils::GTypeCheck->new();

    open (my \$fh,'>',"$name.gtype") or Utils::error("$name.gtype: \$!");
    print \$fh \$gtc->is_genotype_ok("$name.gtypex",'$genotype',$ratio);
    close \$fh;
}
};
    close $fh;

    VertRes::LSF::run($lock_file,$dir,"$prefix${name}_glf",$self, qq{perl -w ${prefix}genotype.pl});

    return $$self{'No'};
}


sub is_genotype_ok
{
    my ($self,$gtype_file,$expected,$min_ratio) = @_;

    open (my $fh,'<',$gtype_file) or $self->throw("$gtype_file; $!");

    # Reading the first two lines would be enough, but we must see if there were
    #   data for the expected genotype.
    #
    #   entropy 16.7
    #   sample C57BL_6J likelihood 9863 over 2583 sites, avg depth 1.05
    #   sample C57BLKS_J likelihood 16203 over 2469 sites, avg depth 1.05

	# We also need to take account of examples where the top 2 records have an equal likelihood
	# and one of them is the expected sample, but as it is listed second it is recorded as unconfirmed
	# e.g. sample UK10K_CIL5062115 lane 6436_7#5:
	#	entropy 0.1, hets included 1
	#	sample UK10K_CIL5002407 likelihood 0 over 26 sites, score 0.000000, avg depth 30.192308
	#	sample UK10K_CIL5062115 likelihood 0 over 26 sites, score 0.000000, avg depth 30.192308

    my $has_data = 0;
    my ($hit1,$hit2,$gtype1,$lhood1,$gtype2,$lhood2);

    my $entrp = <$fh>;
    while (my $line=<$fh>)
    {
        if ( $has_data && defined($hit1) && defined($hit2) ) { last; }

        # The regex matches both lines:
        #   sample NA12717 likelihood 104917 over 13788 sites, score 1.222, avg depth 0.026506
        #   sample xx/NA12717.snp likelihood 104917 over 13788 sites, score 1.222, avg depth 0.026506
        #
        if ( !($line =~ m{sample\s+(?:.*/)?(\S+?)(?:\.snp)?\s+likelihood \d+ over \d+ sites, score (\S+),} ) ) 
        { 
            $self->throw("Could not parse $gtype_file: $hit1") 
        }

        if ( $expected && $1 eq $expected ) { $has_data = 1; } 

        if ( !defined $hit1 ) 
        { 
            $hit1   = $line; 
            $gtype1 = $1;
            $lhood1 = $2;
        }
        elsif ( !defined $hit2 ) 
        { 
            $hit2   = $line; 
            $gtype2 = $1;
            $lhood2 = $2;
        }
    }
    close $fh;

	# need to take account of special circumstance above where the 'correct' sample is listed second, but has 
	# the same likelihood as that reported first->  $expected eq $gtype2 and $lhood1==$lhood2
    my $expected_gtype1 = ($expected eq $gtype1 && $lhood1 == $lhood2) ? 1 : 0;
    my $expected_gtype2 = ($expected eq $gtype2 && $lhood1 == $lhood2) ? 1 : 0; 
     
    if ( $expected && !$has_data ) { $expected = 0; }

    my $ratio = $lhood1!=0 ? $lhood2/$lhood1 : $lhood2/1e-6;

	if ( $expected_gtype1 ) { return "status=confirmed expected=$expected found=$gtype1 ratio=$ratio\n"; }
	if ( $expected_gtype2 ) { return "status=confirmed expected=$expected found=$gtype2 ratio=$ratio\n"; }

    if ( $ratio<$min_ratio ) 
    { 
        if ( $expected ) { return "status=unconfirmed expected=$expected found=$gtype1 ratio=$ratio\n"; }
        return "status=unknown expected=none found=$gtype1 ratio=$ratio\n";
    }
    if ( !$expected ) { return "status=candidate expected=none found=$gtype1 ratio=$ratio\n"; }
    if ( $expected eq $gtype1 ) { return "status=confirmed expected=$expected found=$gtype1 ratio=$ratio\n"; }
    return "status=wrong expected=$expected found=$gtype1 ratio=$ratio\n";
}


sub get_status
{
    my ($gtype_file) = @_;
    my %info = ( status=>'', expected=>'', found=>'', ratio=>'' );

    open(my $fh,'<',$gtype_file) or croak "$gtype_file: $!";
    my (@lines) = <$fh>;
    close($fh) or croak "$gtype_file: $!";
    if ( !scalar @lines ) { croak "Could not read $gtype_file\n"; }

    # status=confirmed expected=NA19107 found=NA19107 ratio=1.28724761001092
    if ( !($lines[0]=~/status=(\S+)\s+expected=(\S+)\s+found=(\S+)\s+ratio=(\S+)/) ) { croak "Could not parse $gtype_file: $lines[0]"; }
    $info{status}   = $1;
    $info{expected} = $2;
    $info{found}    = $3;
    $info{ratio}    = $4;

    return \%info;
}


1;

