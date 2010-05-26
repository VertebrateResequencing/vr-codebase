package FaSlice;

use strict;
use warnings;
use Carp;

sub new
{
    my ($class,@args) = @_;
    my $self = @args ? {@args} : {};
    bless $self, ref($class) || $class;
    if ( !$$self{file} ) { $self->throw("Missing the parameter file\n"); }
    $$self{chr}  = undef;
    $$self{from} = undef;
    $$self{to}   = undef;
    if ( !$$self{size} ) { $$self{size}=1_000_000; }
    return $self;
}

sub throw
{
    my ($self,@msg) = @_;
    confess(@msg);
}

sub cmd
{
    my ($self,$cmd) = @_;
    my @out = `$cmd`;
    if ( $? )
    {
        my @msg = ();
        push @msg, qq[The command "$cmd" returned non-zero status $?];
        if ( $! )
        {
            push @msg, ": $!\n";
        }
        else
        {
            push @msg, ".\n";
        }
        if ( scalar @out )
        {
            push @msg, @out;
        }
        $self->throw(@msg);
    }
    return (@out);
}

sub read_chunk
{
    my ($self,$chr,$pos) = @_;
    my $to = $pos + $$self{size};
    my $cmd = "samtools faidx $$self{file} $chr:$pos-$to";
    my @out = $self->cmd($cmd) or $self->throw("$cmd: $!");
    my $line = shift(@out);
    if ( !($line=~/^>$chr:(\d+)-(\d+)/) ) { $self->throw("Could not parse: $line"); }
    $$self{chr}  = $chr;
    $$self{from} = $1;
    my $chunk = '';
    while ($line=shift(@out))
    {
        chomp($line);
        $chunk .= $line;
    }
    $$self{to} = $$self{from} + length($chunk) - 1;
    $$self{chunk} = $chunk;
    return;
}

sub get_base
{
    my ($self,$chr,$pos) = @_;
    if ( !$$self{chr} || $chr ne $$self{chr} || $pos<$$self{from} || $pos>$$self{to} )
    {
        $self->read_chunk($chr,$pos);
    }
    my $idx = $pos - $$self{from};
    if ( $$self{from}>$$self{to} ) { print STDERR "No such site $chr:$pos\n"; return ''; }
    return substr($$self{chunk},$idx,1);
}


sub get_slice
{
    my ($self,$chr,$from,$to) = @_;
    if ( $to-$from >= $$self{size} ) { $self->throw("Too big region requested, $from-$to >= $$self{size}\n"); }
    if ( !$$self{chr} || $chr ne $$self{chr} || $from<$$self{from} || $to>$$self{to} )
    {
        $self->read_chunk($chr,$from);
    }
    my $idx = $from - $$self{from};
    if ( $$self{from}>$$self{to} || $to>$$self{to} ) { print STDERR "No such region $chr:$from-$to\n"; return ''; }
    return substr($$self{chunk},$idx,$to-$from+1);
}


# http://www.illumina.com/documents/products/technotes/technote_topbot.pdf
sub illumina_alleles_TOP_to_ref
{
    my ($self,$a1,$a2,$chr,$pos,$ref) = @_;
    my %map = (A=>'T', C=>'G', G=>'C', T=>'A');
    my %top = ( 
            A=>{A=>-2,C=> 1,G=> 1,T=>-1}, 
            C=>{A=> 1,C=>-2,G=>-1,T=> 0}, 
            G=>{A=> 1,C=>-1,G=>-2,T=> 0}, 
            T=>{A=>-1,C=> 0,G=> 0,T=>-2} ); 

    my $stat = $top{$a1}{$a2};
    if ( $stat==-2 ) { $self->throw("Expected two different bases, got $a1 and $a2.\n"); }
    if ( $stat==-1 )
    {
        # Now we should do the sequence walking to see if the reference is TOP or BOT,
        #   but we do not this in ill-to-vcf: C/G would become G/C and A/T would become T/A.
        return ($a1,$a2);
    }
    if ( $stat==0 ) { $self->throw("Expected Illumina TOP, got $a1 and $a2.\n"); }
    if ( $ref eq $a1 or $ref eq $a2 ) { return ($a1,$a2); }
    return ($map{$a1},$map{$a2});
}


1;

