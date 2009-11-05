package Vcf;

# http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcfv3.2
#
# Authors: lh3, pd3
# for VCF v3.2

=head1 NAME

Vcf

=head1 SYNOPSIS

From the command line:
    perl -MVcf -e validate example.vcf
    perl -MVcf -ne '$x=vcf1;print if $x->{NA00001}{GT} eq "0/1"' example.vcf

From a script:
    use Vcf;

    # Do some simple parsing 
    my $vcf = Vcf->new(file=>'example.vcf');
    $vcf->parse_header();
    while (my $x=$vcf->parse_next_data_line()) 
    { 
        for my $gt (keys %{$$x{gtypes}})
        {
            print "\t$gt:", $$x{gtypes}{$gt}{GT};
        }
        print "\n";
    }

    # Output a subset of the original file, only the columns NA00001, NA00002 and NA00003
    my $vcf = Vcf->new(file=>'example.vcf');
    $vcf->parse_header();
    my %columns = map { $_ => $i++ } qw(NA00001 NA00002 NA00003);
    print $vcf->sprint_header(\%columns);
    while (my $x=$vcf->parse_next_data_line())
    {
        print $vcf->sprint_line($x,\%columns);     # this will recalculate AC and AN counts
    }

=cut


use strict;
use warnings;
use Carp;
use Exporter;

use vars qw/@ISA @EXPORT %s2i @i2s/;
@ISA = qw/Exporter/;
@EXPORT = qw/vcf1 validate/;

sub vcf1 {
  my %x = ();
  my $i;
  $_ = "$_\n" unless (/\n/);
  if (/^#/) {
	if (/^#CHR/) {
	  @i2s = split;
	  $i2s[0] =~ s/^#//;
	  for (my $i = 0; $i < @i2s; ++$i) {
		$s2i{$i2s[$i]} = $i;
	  }
	}
  } else {
	if (@i2s == 0) {
	  warn("[Vcf::vcf2] fail to find VCF header\n");
	  return;
	}
	my @t = split;
	# mandatory fields
	for my $i (0 .. 6, 8) {
	  $x{$i2s[$i]} = $t[$i];
	}
	# INFO
	my $z = \%{$x{$i2s[7]}};
	$t[7] =~ s/([^\s=;]+)(=([^\s;]+))?/$z->{$1} = defined($3)? $3 : 1/eg;
	# FORMAT
	my @s = split(':', $t[8]);
	# samples
	for $i (9 .. $#t) {
	  my $j = 0;
	  $z = \%{$x{$i2s[$i]}};
	  $t[$i] =~ s/([^\s:]+)/$z->{$s[$j++]}=$1;/eg;
	}
  }
  return \%x;
}

sub validate
{
    my $vcf = $ARGV[0] ? Vcf->new(file=>$ARGV[0]) : Vcf->new();
    $vcf->parse_header();
    if ( !$vcf->{header} ) { $vcf->warn("No VCF header found.\n"); }
    if ( !$vcf->{columns} ) { $vcf->warn("No column descriptions found.\n"); }

    my $warn_sorted=1;
    my ($prev_chrm,$prev_pos);
    while (my $x=$vcf->parse_next_data_line()) 
    {
        if ( $warn_sorted )
        {
            if ( $prev_chrm && $prev_chrm eq $$x{CHROM} && $prev_pos > $$x{POS} ) 
            { 
                $vcf->warn("The file not sorted, seee e.g. $prev_chrm:$prev_pos and $$x{CHROM}:$$x{POS}\n");
                $warn_sorted = 0;
            }
            $prev_chrm = $$x{CHROM};
            $prev_pos  = $$x{POS};
        }
        if ( exists($$x{INFO}{AN}) || exists($$x{INFO}{AC}) )
        {
            my ($an,$ac) = $vcf->calc_an_ac($$x{gtypes});
            if ( exists($$x{INFO}{AN}) && $an ne $$x{INFO}{AN} ) 
            { 
                $vcf->throw("$$x{CHROM}:$$x{POS} .. AN is $$x{INFO}{AN}, should be $an\n"); 
            }
            if ( exists($$x{INFO}{AC}) && $ac ne $$x{INFO}{AC} ) 
            { 
                $vcf->throw("$$x{CHROM}:$$x{POS} .. AC is $$x{INFO}{AC}, should be $ac\n"); 
            }
        }
    }
}


sub new
{
    my ($class,@args) = @_;
    my $self = {@args};
    bless($self);
    if ( !exists($$self{fh}) ) 
    { 
        if ( exists($$self{file}) ) 
        { 
            open($$self{fh},'<',$$self{file}) or $self->throw("$$self{file}: $!"); 
        }
        else { $$self{fh} = *STDIN; }
    }
    $$self{buffer}  = [];       # buffer stores the lines in the reverse order
    $$self{header}  = undef;
    $$self{columns} = undef;    # column names 
    $$self{mandatory} = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']; 
    return $self;
}

sub sprint_header
{
    my ($self,$columns) = @_;
    my $out;
    if ( $$self{header} ) 
    {
        while (my ($key,$value) = each %{$$self{header}})
        {
            $out .= "##" . $key .'='. $value ."\n";
        }
    }
    if ( $$self{columns} )
    {
        my @out ;
        if ( $columns )
        {
            my $icol=0;
            for my $col (@{$$self{columns}})
            {
                if ( $icol++>8 && !exists($$columns{$col}) ) { next; }
                push @out, $col;
            }
        }
        else 
        { 
            @out = @{$$self{columns}}; 
        }
        $out .= "#". join("\t", @out). "\n";
    }
    return $out;
}

sub sprint_line
{
    my ($self,$record,$columns) = @_;
    # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA00001 NA00002 NA00003
    my $out = $$record{CHROM} ."\t". $$record{POS} ."\t". $$record{ID} ."\t". $$record{REF};
    $out .= "\t". join(',',@{$$record{ALT}});
    $out .= "\t". $$record{QUAL} ."\t". $$record{FILTER};

    my %gtypes;
    if ( $columns )
    {
        %gtypes = map { $_ => $$record{gtypes}{$_} } sort { $$columns{$a} <=> $$columns{$b} } keys %$columns;
    }
    else
    {
        %gtypes = %{$$record{gtypes}};
    }
    my ($an,$ac) = $self->calc_an_ac(\%gtypes);
    my @info = ("AN=$an","AC=$ac");
    while (my ($key,$value) = each %{$$record{INFO}})
    {
        if ( $key eq 'AN' ) { next; }
        if ( $key eq 'AC' ) { next; }
        push @info, (defined $value ? "$key=$value" : $key);
    }
    $out .= "\t". join(';', @info);
    $out .= "\t". join(':', @{$$record{FORMAT}});

    for my $gt (values %gtypes)
    {
        $out .= "\t" . join(':', map { exists($$gt{$_}) ? $$gt{$_} : '' } @{$$record{FORMAT}});
    }

    $out .= "\n";
    return $out;
}

sub parse_header
{
    my ($self) = @_;

    # First come the header lines prefixed by ##
    while ($self->next_header_line()) { ; }

    # Now comes the column names line prefixed by #
    $self->read_column_names();
}


# This contains quite thorough validation. If too slow, a simpler parser can be added.
sub parse_next_data_line
{
    my ($self) = @_;
    my $line = $self->get_line();
    if ( !$line ) { return undef; }
    chomp($line);
    my @cols = split(/\t/, $line);
    if ( $$self{columns} && scalar @{$$self{columns}} != scalar @cols )
    {
        my $msg = "The number of columns (".scalar @cols
            .") does not match the number of column names (" .scalar @{$$self{columns}}. ").\n"
            ."The offendig line was [$line].\n";
        $self->throw($msg);
    }
    my %out = (); 
    my $ncols = scalar @cols;
    my ($format,$alts);
    if ( $ncols<8 ) { $self->throw("Missing some mandatory columns on line [$line].\n"); } 
    for (my $i=0; $i<$ncols; $i++)
    {
        my $cname = $$self{columns} ? $$self{columns}->[$i] : $i;
        my $value = $cols[$i];
        if ( $i==1 ) 
        {
            if ( !($value=~/^\d+$/) ) 
            { 
                $self->throw("Error parsing the POS field [$value] on line [$line]\n"); 
            }
        }
        elsif ( $i==3 ) 
        {
            if ( !($value=~/^[ACGTN]$/) ) 
            { 
                $self->throw("Error parsing the REF field [$value] on line [$line]\n"); 
            }
        }
        elsif ( $i==4 ) 
        {
            $value = $self->parse_alt_field($value);
            $alts  = $value;
            if ( !defined($value) )
            {
                $self->throw("Error parsing the ALT field [$cols[$i]] on line [$line].\n");
            }
        }
        elsif ( $i==7 ) 
        { 
            $value = $self->parse_info_field($value); 
            if ( !defined($value) )
            {
                $self->throw("Error parsing the INFO field [$cols[$i]] on line [$line].\n");
            }
        }
        elsif ( $i==8 ) 
        { 
            $value  = $self->parse_format_field($value); 
            $format = $value;
            if ( !defined($value) )
            {
                $self->throw("Error parsing the FORMAT field [$cols[$i]] on line [$line].\n");
            }
        }
        elsif ( $i>8 )
        {
            $value = $self->parse_gtype_field($value, $format, $alts);
            if ( !defined($value) )
            {
                $self->throw("Error parsing the column $cname [$cols[$i]] on line [$line].\n");
            }
            $out{gtypes}{$cname} = $value;
            next;
        }
        $out{$cname} = $value; 
    }
    return \%out;
}

sub parse_alt_field
{
    my ($self,$value) = @_;

    my @items = split(/,/,$value);
    if ( scalar @items>1 && $value=~/\./ ) 
    { 
        $self->warn("[parse_alt_field] A mixture of multiple alleles and '.'\n");
        return undef; 
    }

    for my $item (@items)
    {
        if ( $item eq '.' ) { next; }
        elsif ( $item=~/^[ACTGN]$/ ) { next; }
        elsif ( $item=~/^I[ACTGN]+$/ ) { next; }
        elsif ( $item=~/^D\d+$/ ) { next; }
        $self->warn("[parse_al_field] Unrecognised format for ALT.\n");
        return undef;
    }
    return \@items;
}


sub parse_info_field
{
    my ($self,$value) = @_;
    my %out;
    my %infos = ( AA=>1, AC=>1, AN=>1, AF=>2, DP=>1, MQ=>1, NS=>1, BQ=>1, SB=>1, DB=>0, H2=>0 );
    my @items = split(/;/,$value);
    for my $item (@items)
    {
        my ($key,$value) = split(/=/,$item);
        if ( !exists($infos{$key}) ) 
        { 
            $self->warn("[parse_info_field] No such key for INFO: [$key]\n");
            return undef; 
        }
        if ( $infos{$key}>0 && !defined($value) ) 
        { 
            $self->warn("[parse_info_field] Expected a value for [$key]\n");
            return undef; 
        }
        $out{$key} = $value;
    }
    return \%out;
}


sub parse_format_field
{
    my ($self,$value) = @_;
    my %fields = ( GT=>1, GQ=>1, DP=>1, HQ=>1, FT=>1 );
    my @items  = split(/:/, $value);
    for my $item (@items)
    {
        if ( !exists($fields{$item}) ) 
        { 
            $self->warn("[parse_info_field] No such key for FORMAT: [$item]\n");
            return undef; 
        }
    }
    if ( $items[0] ne 'GT' ) 
    { 
        $self->warn("[parse_info_field] Wrong order, GT must be first\n");
        return undef; 
    }
    return \@items;
}

sub parse_gtype_field
{
    my ($self, $value, $format, $alts) = @_;
    my %out;
    my @items = split(/:/, $value);
    if ( scalar @items > scalar @$format ) 
    { 
        $self->warn("[parse_gtype_field] Incomplete FORMAT.\n");
        return undef; 
    }
    for (my $i=0; $i<@items; $i++)
    {
        if ( !defined($items[$i]) ) 
        { 
            # The field can be empty
            $out{$$format[$i]} = $items[$i];
            next;
        }

        if ( $$format[$i] eq 'GT' )
        {
            if ( $items[$i] eq '.' ) # This is only to convert existing malformatted vcf's
            { 
                $items[$i] = './.';
            }
            elsif ( !($items[$i]=~m{^(?:\.|(\d+))[\|/](?:\.|(\d+))$}) ) 
            { 
                # Dots or two numbers separated by either of \|/
                $self->warn("[parse_gtype_field] Could not parse GT [$items[$i]]\n");
                return undef; 
            }
            if ( (defined($1) && $1>@$alts) || (defined($2) && $2>@$alts) )
            {
                $self->warn("[parse_gtype_field] The gtype value too big in [$items[$i]]\n");
                return undef;
            }
        }
        elsif ( $$format[$i] eq 'GQ' )
        {
            if ( !($items[$i]=~/^\d+$/) )
            {
                $self->warn("[parse_gtype_field] Could not parse GQ [$items[$i]], expected number\n");
                return undef;
            }
            if ( $items[$i]>99 )
            {
                $self->warn("[parse_gtype_field] Could not parse GQ [$items[$i]], value bigger than 99\n");
                return undef;
            }
        }
        elsif ( $$format[$i] eq 'DP' && !($items[$i]=~/^\d+$/) )
        {
            $self->warn("[parse_gtype_field] Could not parse DP [$items[$i]], expected number\n");
            return undef;
        }
        elsif ( $$format[$i] eq 'HQ' && !($items[$i]=~/^\d+,\d+$/) )
        {
            $self->warn("[parse_gtype_field] Could not parse HQ [$items[$i]], expected number,number\n");
            return undef;
        }
        elsif ( $$format[$i] eq 'FT' && !($items[$i]=~/^[,0-9]+$/) )
        {
            $self->warn("[parse_gtype_field] Could not parse FT [$items[$i]], expected colon separated numbers\n");
            return undef;
        }
        $out{$$format[$i]} = $items[$i];
    }
    return \%out;
}

sub calc_an_ac
{
    my ($self,$gtypes) = @_;
    my ($an,%ac_counts);
    $an = 0;
    for my $gt (keys %$gtypes)
    {
        my $value = $$gtypes{$gt}{GT};
        if ( $value eq './.' ) { next; }
        my ($al1,$al2) = split(m{[\|/]},$value);
        $an += 2;
        $ac_counts{$al1}++ unless $al1 eq '0';
        $ac_counts{$al2}++ unless $al2 eq '0';
    }
    my @ac = map { $ac_counts{$_} } sort { $a <=> $b } keys %ac_counts;
    if ( !@ac ) { @ac = '0'; }
    return ($an,join(',',@ac));
}


# Stores the header lines in a hash
sub next_header_line
{
    my ($self) = @_;
    my $line = $self->get_line();
    if ( substr($line,0,2) ne '##' )
    {
        $self->unread_line($line);
        return undef;
    }
    if ( !($line=~/^\#\#(\S+)=(.*)$/) ) # allowing white space 
    { 
        chomp($line);
        $self->throw("Could not parse the header line: [$line]\n"); 
    }
    if ( exists($$self{header}{$1}) )
    {
        chomp($line);
        if ( $$self{header}{$1} eq $2 ) 
        { 
            warn("Duplicate header line [$line].\n"); 
        }
        $self->throw("Conflicting header line [$1=$$self{header}{$1}] vs [$1=$2].");
    }
    $$self{header}{$1} = $2;
    return $line;
}

sub read_column_names
{
    my ($self) = @_;
    my $line = $self->get_line();
    if ( substr($line,0,1) ne '#' || substr($line,1,1) eq '#' )
    {
        $self->unread_line($line);
        return undef;
    }
    chomp($line);
    my @cols  = split(/\t/, substr($line,1));
    my $ncols = scalar @cols;
    if ( $ncols == 1 )
    {
        # If there is only one name, it can be space-seprated instead of tab separated
        my @items = split(/\s+/, $cols[0]);
        if ( scalar @items >= 8 )
        {
            $self->throw(
                "Too few column names in [$line].\n"
                . "Hint: the column names must be tab separated.\n");
        }
    }

    # Check the names of the mandatory columns
    if ( $ncols < 8 ) { $self->throw("Mssing mandatory columns in [$line].\n"); }
    my $fields = $$self{mandatory};
    for (my $i=0; $i<$ncols; $i++)
    {
        if ( @$fields <= $i ) { last; }
        if ( $cols[$i] ne $$fields[$i] ) 
        { 
            $self->throw("Expected mandatory column [$$fields[$i]], got [$cols[$i]]\n"); 
        }
    }
    $$self{columns} = \@cols;
    return $$self{columns};
}

sub get_line
{
    my ($self) = @_;
    if ( @{$$self{buffer}} ) { return shift(@{$$self{buffer}}); }
    my $fh = $$self{fh};
    my $line = <$fh>;
    return $line;
}

sub unread_line
{
    my ($self,$line) = @_;
    unshift @{$$self{buffer}}, $line;
    return;
}

sub throw
{
    my ($self,@msg) = @_;
    croak @msg;
}

sub warn
{
    my ($self,@msg) = @_;
    warn @msg;
}

1;

