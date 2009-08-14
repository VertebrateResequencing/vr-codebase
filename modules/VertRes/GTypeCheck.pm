package VertRes::GTypeCheck;
use base qw(VertRes::Base);

use strict;
use warnings;
use LSF;
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

    # Dynamic script to be run by LSF.
    open(my $fh, '>', "$dir/${prefix}genotype.pl") or $self->throw("$dir/${prefix}genotype.pl: $!");
    print $fh
    qq{
use Utils;
use VertRes::GTypeCheck;

\$base = VertRes::Base->new();

if ( ! -e "$bam.bai" || Utils::file_newer("$bam","$bam.bai") ) 
{ 
    Utils::CMD("$samtools index $bam"); 
}
if ( ! -e "$name.glf" || Utils::file_newer("$bam","$name.glf") )
{
    Utils::CMD("$samtools pileup -g -f $fa_ref $bam > $name.glfx");
    rename("$name.glfx", "$name.glf") or Utils::CMD("rename $name.glfx $name.glf: \$!");
}
if ( ! -e "$name.gtypex" || Utils::file_newer("$name.glf","$name.gtypex") )
{
    Utils::CMD("$glf checkGenotype $snps $name.glf > $name.gtypey");
    if ( ! -s "$name.gtypey" ) 
    { 
        \$base->throw("FIXME: this should not happen, what's wrong with:\n\t$glf checkGenotype $snps $name.glf > $name.gtypey\n???\n");
    }
    rename("$name.gtypey", "$name.gtypex") or Utils::CMD("rename $name.gtypey $name.gtypex: \$!");
    if ( ! -s "$name.gtypex" ) 
    { 
        \$base->throw("FIXME: this should not happen, what's wrong with:\n\trename $name.gtypey $name.gtypex\n???\n");
    }
}
if ( -s "$name.gtypex" )
{
    my \$gtc = VertRes::GTypeCheck->new();

    open (my \$fh,'>',"$name.gtype") or Utils::error("$name.gtype: \$!");
    print \$fh \$gtc->is_genotype_ok("$name.gtypex",'$genotype',$ratio);
    close \$fh;
}
};
    close $fh;

    LSF::run($lock_file,$dir,"$prefix${name}_glf",$self, qq{perl -w ${prefix}genotype.pl});

    return $$self{'No'};
}


sub is_genotype_ok
{
    my ($self,$gtype_file,$expected,$min_ratio) = @_;

    open (my $fh,'<',$gtype_file) or $self->throw("$gtype_file; $!");

    # entropy 16.7
    # sample C57BL_6J likelihood 9863 over 2583 sites, avg depth 1.05
    # sample C57BLKS_J likelihood 16203 over 2469 sites, avg depth 1.05

    my $entrp = <$fh>;
    my $hit1  = <$fh>;
    my $hit2  = <$fh>;
    close $fh;

    # The regex matches both lines:
    #   sample NA12717 likelihood 104917 over 13788 sites, avg depth 0.026506
    #   sample xx/NA12717.snp likelihood 104917 over 13788 sites, avg depth 0.026506
    #
    if ( !($hit1 =~ m{sample\s+(?:.*/)?(\S+?)(?:\.snp)?\s+likelihood (\d+)} ) ) { $self->throw("Could not parse $gtype_file: $hit1") }
    my $gtype1 = $1;
    my $lhood1 = $2;
    if ( !($hit2 =~ m{sample\s+(?:.*/)?(\S+?)(?:\.snp)?\s+likelihood (\d+)} ) ) { $self->throw("Could not parse $gtype_file: $hit2") }
    my $gtype2 = $1;
    my $lhood2 = $2;

    my $ratio = $lhood1!=0 ? $lhood2 / $lhood1 : 0;

    if ( $ratio<$min_ratio ) 
    { 
        if ( $expected ) { return "status=unconfirmed expected=$expected found=none ratio=$ratio\n"; }
        return "status=unknown expected=none found=none ratio=$ratio\n";
    }
    if ( !$expected ) { return "status=candidate expected=none found=$gtype1 ratio=$ratio\n" }
    if ( $expected eq $gtype1 ) { return "status=confirmed expected=$expected found=$gtype1 ratio=$ratio\n" }
    return "status=wrong expected=$expected found=$gtype1 ratio=$ratio\n";
}

#---------- Debugging and error reporting -----------------

#   sub warn
#   {
#       my ($self,@msg) = @_;
#       if ($self->verbose > 0) 
#       {
#           print STDERR join('',@msg);
#       }
#       $self->log(join('',@msg));
#   }
#   
#   sub debug
#   {
#       my ($self,@msg) = @_;
#       if ($self->verbose > 0) 
#       {
#           print STDERR join('',@msg);
#           $self->log(join('',@msg));
#       }
#   }
#   
#   sub throw
#   {
#       my ($self,@msg) = @_;
#       Utils::error(@msg);
#   }


1;

