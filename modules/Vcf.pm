package Vcf;

# http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcfv3.2
# http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf3.3
#
# Authors: lh3, pd3
# for VCF v3.2 and v3.3

=head1 NAME

Vcf

=head1 SYNOPSIS

From the command line:
    perl -MVcf -e validate example.vcf

From a script:
    use Vcf;

    my $vcf = Vcf->new(file=>'example.vcf.gz');
    $vcf->parse_header();

    # Do some simple parsing. Most thorough but slowest way how to get the data.
    while (my $x=$vcf->next_data_hash()) 
    { 
        for my $gt (keys %{$$x{gtypes}})
        {
            print "\t$gt:", $$x{gtypes}{$gt}{GT};
        }
        print "\n";
    }

    # This will split the fields and print a list of CHR:POS
    while (my $x=$vcf->next_data_array()) 
    {
        print "$$x[0]:$$x[1]\n";
    }

    # This will return the lines as they were read, including the newline at the end
    while (my $x=$vcf->next_line()) 
    { 
        print $x;
    }

    # Output a subset of the original file in the v3.2 format. Only the columns 
    #   NA00001, NA00002 and NA00003 will be printed.
    my @columns = qw(NA00001 NA00002 NA00003);
    $vcf->output_version('3.2');
    print $vcf->format_header(\@columns);
    while (my $x=$vcf->next_data_array())
    {
        # this will recalculate AC and AN counts, unless $vcf->recalc_ac_an was set to 0
        print $vcf->format_line($x,\@columns); 
    }

=cut


use strict;
use warnings;
use Carp;
use Exporter;

use vars qw/@ISA @EXPORT %s2i @i2s/;
@ISA = qw/Exporter/;
@EXPORT = qw/vcf1 validate/;


# This is the original code, not used, left only for backward compatibility.
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



=head2 validate

    About   : Validates the VCF file.
    Usage   : perl -MVcf -e validate example.vcf.gz     # (from the command line)
              validate('example.vcf.gz');               # (from a script)
              validate(\*STDIN);
    Args    : File name or file handle. When no argument given, the first command line
              argument is interpreted as the file name.

=cut

sub validate
{
    my ($fh) = @_;

    if ( !$fh && @ARGV && -e $ARGV[0] ) { $fh = $ARGV[0]; }

    my $vcf;
    if ( $fh ) { $vcf = fileno($fh) ? Vcf->new(fh=>$fh) : Vcf->new(file=>$fh); }
    else { $vcf = Vcf->new(); }

    $vcf->parse_header();
    if ( !exists($$vcf{header_lc}{fileformat}) ) 
    { 
        $vcf->warn(qq[The "fileformat" field not present in the header, assuming VCFv$$vcf{version}\n]);
    }
    if ( !$vcf->{header} ) { $vcf->warn("No VCF header found.\n"); }
    if ( !$vcf->{columns} ) { $vcf->warn("No column descriptions found.\n"); }

    my $warn_sorted=1;
    my ($prev_chrm,$prev_pos);
    while (my $x=$vcf->next_data_hash()) 
    {
        # Is the position numeric?
        if ( !($$x{POS}=~/^\d+$/) ) { $vcf->warn("Expected integer for the position at $$x{CHROM}:$$x{POS}\n"); }

        # Is the file sorted?
        if ( $warn_sorted )
        {
            if ( $prev_chrm && $prev_chrm eq $$x{CHROM} && $prev_pos > $$x{POS} ) 
            { 
                $vcf->warn("The file is not sorted, for example $$x{CHROM}:$$x{POS} comes after $prev_chrm:$prev_pos\n");
                $warn_sorted = 0;
            }
            $prev_chrm = $$x{CHROM};
            $prev_pos  = $$x{POS};
        }

        # Is the ID non-empty?
        if ( $$x{ID}=~/^$/ ) { $vcf->warn("The ID should be set to . if unknown at $$x{CHROM}:$$x{POS}\n"); }

        # The reference base: one of A,C,G,T,N, non-empty.
        if ( !($$x{REF}=~/^[ACGTN]$/) ) { $vcf->warn("Expected one of A,C,G,T,N for the reference base at $$x{CHROM}:$$x{POS}, got [$$x{REF}]\n"); }
        
        # The ALT field (alternate non-reference base)
        my $err = $vcf->validate_alt_field($$x{ALT});
        if ( $err ) { $vcf->warn("$$x{CHROM}:$$x{POS} .. $err\n"); }

        # The QUAL field
        my $ret = $vcf->validate_float($$x{QUAL},-1);
        if ( $ret ) { $vcf->warn("QUAL field at $$x{CHROM}:$$x{POS} .. $ret\n"); }

        # The FILTER field
        $err = $vcf->validate_filter_field($$x{FILTER});
        if ( $err ) { $vcf->warn("FILTER field at $$x{CHROM}:$$x{POS} .. $err\n"); }

        # The INFO field
        $err = $vcf->validate_info_field($$x{INFO});
        if ( $err ) { $vcf->warn("INFO field at $$x{CHROM}:$$x{POS} .. $err\n"); } 

        while (my ($gt,$data) = each %{$$x{gtypes}})
        {
            $err = $vcf->validate_gtype_field($data,$$x{ALT});
            if ( $err ) { $vcf->warn("$gt column at $$x{CHROM}:$$x{POS} .. $err\n"); }
        }

        if ( exists($$x{INFO}{AN}) || exists($$x{INFO}{AC}) )
        {
            my ($an,$ac) = $vcf->calc_an_ac($$x{gtypes});
            if ( exists($$x{INFO}{AN}) && $an ne $$x{INFO}{AN} ) 
            { 
                $vcf->warn("$$x{CHROM}:$$x{POS} .. AN is $$x{INFO}{AN}, should be $an\n"); 
            }
            if ( exists($$x{INFO}{AC}) && $ac ne $$x{INFO}{AC} ) 
            { 
                $vcf->warn("$$x{CHROM}:$$x{POS} .. AC is $$x{INFO}{AC}, should be $ac\n"); 
            }
        }
    }
}

=head2 new

    About   : Creates new VCF reader/writer. 
    Usage   : my $vcf = Vcf->new(file=>'my.vcf', version=>'3.2');
    Args    : 
                file    .. The file name. If not given, STDIN is assumed.
                silent  .. Unless set to 0, warning messages may be printed.
                strict  .. Unless set to 0, the reader will die when the file violates the specification.
                version .. If not given, '3.2' is assumed.

=cut

sub new
{
    my ($class,@args) = @_;
    my $self = {@args};
    bless($self);
    if ( !exists($$self{fh}) ) 
    { 
        if ( exists($$self{file}) ) 
        {
            if ( $$self{file}=~/\.gz$/i )
            {
                open($$self{fh},"zcat $$self{file} |") or $self->throw("$$self{file}: $!");
            }
            else
            {
                open($$self{fh},'<',$$self{file}) or $self->throw("$$self{file}: $!"); 
            }
        }
        else { $$self{fh} = *STDIN; }
    }
    $$self{version}   = '3.3';
    $$self{silent}    = 0 unless exists($$self{silent});
    $$self{strict}    = 0 unless exists($$self{strict});
    $$self{buffer}    = [];       # buffer stores the lines in the reverse order
    $$self{header_lc} = {};
    $$self{columns}   = undef;    # column names 
    $$self{mandatory} = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] unless exists($$self{mandatory}); 
    $$self{recalc_ac_an} = 1;
    return $self;
}

=head2 next_line

    About   : Reads next VCF line.
    Usage   : my $vcf = Vcf->new(); 
              my $x   = $vcf->next_line();
    Args    : none

=cut

sub next_line
{
    my ($self) = @_;
    if ( @{$$self{buffer}} ) { return shift(@{$$self{buffer}}); }
    return readline($$self{fh});
}

=head2 next_data_array

    About   : Reads next VCF line and splits it into an array. The last element is chomped.
    Usage   : my $vcf = Vcf->new(); 
              $vcf->parse_header(); 
              my $x = $vcf->next_data_array();
    Args    : none

=cut

sub next_data_array
{
    my ($self) = @_;
    my $line;
    if ( @{$$self{buffer}} ) { $line = shift(@{$$self{buffer}}); }
    else { $line = readline($$self{fh}); }
    if ( !$line ) { return undef; }
    my @items = split(/\t/,$line);
    chomp($items[-1]);
    return \@items;
}

=head2 next_data_hash

    About   : Reads next VCF line and splits it into a hash. This is the slowest way to obtain the data.
    Usage   : my $vcf = Vcf->new(); 
              $vcf->parse_header(); 
              my $x = $vcf->next_data_hash();
    Args    : none

=cut

sub next_data_hash
{
    my ($self) = @_;
    my $line;
    if ( @{$$self{buffer}} ) { $line = shift(@{$$self{buffer}}); }
    else { $line = readline($$self{fh}); }
    if ( !$line ) { return undef; }
    my @items = split(/\t/,$line);
    chomp($items[-1]);

    if ( !$$self{columns} ) { $self->_fake_column_names(scalar @items); }
    my $cols = $$self{columns};
    my %out;

    # Mandatory fields
    for my $i (0..3,5) { $out{$$cols[$i]} = $items[$i]; }

    # ALT field
    $out{$$cols[4]} = [ split(/,/,$items[4]) ];

    # FILTER field
    $out{$$cols[6]} = [ split(/;/,$items[6]) ];

    # Info, e.g. NS=58;DP=258;AF=0.786;DB;H2
    if ( $items[7] )
    {
        my %hash;
        for my $info (split(/;/,$items[7]))
        {
            my ($key,$val) = split(/=/,$info);
            $hash{$key} = $val || $$self{header}{$$cols[7]}{$key}{default};
        }
        $out{$$cols[7]} = \%hash;
    }

    # Format, e.g. GT:GQ:DP:HQ
    my $format = $out{$$cols[8]} = [ split(/:/,$items[8]) ];

    # Genotype fields
    my %gtypes;
    my $check_nformat = $$self{version} < 3.3 ? 0 : 1;
    for (my $icol=9; $icol<@items; $icol++)
    {
        my @fields = split(/:/, $items[$icol]);
        if ( $check_nformat && @fields != @$format ) 
        { 
            $self->warn("Different number of fields in the format spec and the column $$cols[$icol] at $items[0]:$items[1] ("
                .scalar @fields." vs ".scalar @$format.")\n"); 
        }
        my %hash;
        for (my $ifield=0; $ifield<@fields; $ifield++)
        {
            $hash{$$format[$ifield]} = $fields[$ifield];
        }
        $gtypes{$$cols[$icol]} = \%hash;
    }
    $out{gtypes} = \%gtypes;

    return \%out;
}

sub _unread_line
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
    if ( $$self{silent} ) { return; }
    if ( $$self{strict} ) { $self->throw(@msg); }
    warn @msg;
}


=head2 parse_header

    About   : Reads (and stores) the VCF header.
    Usage   : my $vcf = Vcf->new(); $vcf->parse_header();
    Args    : none

=cut

sub parse_header
{
    my ($self) = @_;

    # First come the header lines prefixed by ##
    while ($self->_next_header_line()) { ; }

    # Now comes the column names line prefixed by #
    $self->_read_column_names();
}


=head2 _next_header_line

    About   : Stores the header lines and meta information, such fields types, etc.
    Args    : none

=cut

sub _next_header_line
{
    my ($self) = @_;
    my $line = $self->next_line();
    if ( !defined $line ) { return undef; }
    if ( substr($line,0,2) ne '##' )
    {
        $self->_unread_line($line);
        return undef;
    }
    push @{$$self{header_lines}}, $line;

    if ( !($line=~/^\#\#(\S+)=(.*)$/) )
    { 
        chomp($line);
        $self->warn("Could not parse the header line: [$line]\n"); 
    }
    if ( $1 eq 'INFO' ) { $self->_add_field($1,$2); }
    elsif ( $1 eq 'FILTER' ) { $self->_add_filter_field($2); }
    elsif ( $1 eq 'FORMAT' ) { $self->_add_field($1,$2); }
    else 
    {
        $$self{header}{$1} = $2; 
        my $key = lc($1);
        $$self{header_lc}{$key} = $2;

        # If the header line contains the file format version info, parse it and save.
        if ( $key eq 'format' ) { $key='fileformat'; }
        if ( $key eq 'fileformat' ) 
        { 
            my $value = $2;
            if ( !($value=~/(\d+(?:\.\d+)?)$/) ) 
            { 
                $self->warn("Could not parse the fileformat version string [$value], assuming VCFv3.3\n");
                $$self{version} = '3.3';
            }
            else { $$self{version} = $1; }
        }
    }

    return $line;
}


=head2 _add_field

    About   : Stores the field types of INFO or FORMAT
    Usage   : $vcf->_add_field('INFO', q[AA,1,String,"Ancestral Allele"]);
              $vcf->_add_field('FORMAT', q[GT,1,String,"Genotype"]);
    Args    : name, value 

=cut

sub _add_field
{
    my ($self,$field,$string) = @_;
    my @values = split(/,/,$string);
    if ( @values < 3 ) { $self->throw("Could not parse [$field=$string].\n"); }
    elsif ( @values < 4 ) { $self->warn("No description in [$field=$string].\n"); } 

    if ( !($values[1]=~/^-?\d+$/) ) { $self->throw("Expected number of arguments in [$field=$string].\n"); }

    my $missing;
    my $handler;
    if ( $values[2] eq 'Integer' ) { $missing=-1; $handler = \&validate_int; }
    elsif ( $values[2] eq 'Float' ) { $missing=-1; $handler = \&validate_float; }
    elsif ( $values[2] eq 'Character' ) { $missing='.'; $handler = \&validate_char; }
    elsif ( $values[2] eq 'String' ) { $missing='.'; }
    elsif ( $values[2] eq 'Flag' ) { }
    else { $self->throw("Unknown field type [$field=$string].\n"); }

    unshift(@values,$missing); # prepend with the default (missing) value

    if ( exists($$self{header}{$field}{$values[1]}) ) { $self->warn("The field specified twice [$field=$string].\n"); }
    $$self{header}{$field}{$values[1]} = 
    {
        default => $values[0],
        name    => $values[1],
        nparams => $values[2],
        type    => $values[3],
        desc    => $values[4],
        handler => $handler,
    };
}


=head2 _add_filter_field

    About   : Stores the field types of FILTER
    Usage   : $vcf->_add_filter_field(q[q10,"Quality below 10"]);
    Args    : name, value 

=cut

sub _add_filter_field
{
    my ($self,$string) = @_;
    my @values = split(/,/,$string);
    if ( @values < 1 ) { $self->throw("Could not parse [FILTER=$string].\n"); }
    elsif ( @values < 2 ) { $self->warn("No description in [FILTER=$string].\n"); } 

    if ( exists($$self{header}{'FILTER'}{$values[0]}) ) { $self->warn("The field specified twice [FILTER=$string].\n"); }
    $$self{header}{'FILTER'}{$values[0]} = 
    {
        name    => $values[1],
        desc    => $values[4],
    };
}


=head2 _read_column_names

    About   : Stores the columns names as array $$self{columns} and hash $$self{has_column}{COL_NAME}=index.
              The indexes goes from 1.
    Usage   : $vcf->_read_column_names();
    Args    : none

=cut

sub _read_column_names
{
    my ($self) = @_;
    my $line = $self->next_line();
    if ( !defined $line ) { return undef; }
    if ( substr($line,0,1) ne '#' || substr($line,1,1) eq '#' )
    {
        $self->_unread_line($line);
        return undef;
    }
    chomp($line);
    my @cols  = split(/\t/, substr($line,1));
    my $ncols = scalar @cols;
    if ( $ncols == 1 )
    {
        # If there is only one name, it can be space-seprated instead of tab separated
        @cols  = split(/\s+/, $cols[0]);
        $ncols = scalar @cols;
        chomp($line);
        if ( $ncols <= 1 ) { $self->warn("Could not parse the column names. [$line]\n"); return; }
        $self->warn("The column names not tab-separated? [$line]\n");
    }

    my $fields  = $$self{mandatory};
    my $nfields = scalar @$fields;

    # Check the names of the mandatory columns
    if ( $ncols < $nfields ) { chomp($line); $self->warn("Missing mandatory column names. [$line].\n"); return; }

    for (my $i=0; $i<$ncols; $i++)
    {
        if ( $i<$nfields && $cols[$i] ne $$fields[$i] ) 
        { 
            $self->warn("Expected mandatory column [$$fields[$i]], got [$cols[$i]]\n"); 
        }
        $$self{has_column}{$cols[$i]} = $i+1;
    }
    $$self{columns} = \@cols;
    return;
}


=head2 format_header

    About   : When no header is present, fake column names as the default mandatory ones + numbers
    Usage   : print $vcf->format_header();
    Args    : The number of columns total (i.e. including the mandatory columns)

=cut

sub _fake_column_names
{
    my ($self,$ncols) = @_;

    $$self{columns} = [ @{$$self{mandatory}} ];
    my $i = scalar @{$$self{columns}};
    while ($i<$ncols) { push @{$$self{columns}}, $i-8; $i++; }
}


=head2 format_header

    About   : Returns the header.
    Usage   : print $vcf->format_header();
    Args    : none

=cut

sub format_header
{
    my ($self,$columns) = @_;

    my $out;
    if ( $$self{header_lines} ) 
    {
        for my $line (@{$$self{header_lines}}) { $out .= $line; }
    }

    if ( !$$self{columns} ) { return $out; }

    my @out_cols;
    if ( $columns )
    {
        @out_cols = @{$$self{columns}}[0..8];
        for my $col (@$columns)
        {
            if ( exists($$self{has_column}{$col}) ) { push @out_cols, $col; }
        }
    }
    else 
    { 
        @out_cols = @{$$self{columns}}; 
    }
    $out .= "#". join("\t", @out_cols). "\n";
    
    return $out;
}


=head2 format_line

    About   : Returns the header.
    Usage   : $x = $vcf->next_data_hash(); print $vcf->format_line($x);
              $x = $vcf->next_data_array(); print $vcf->format_line($x);
    Args 1  : The columns or hash in the format returned by next_data_hash or next_data_array.
         2  : The columns to include [optional]

=cut

sub format_line
{
    my ($self,$record,$columns) = @_;

    if ( ref($record) eq 'HASH' ) { return $self->_format_line_hash($record,$columns); }
    $self->throw("FIXME: todo\n");
}


=head2 recalc_ac_an

    About   : Control if the AC and AN values should be updated.
    Usage   : $vcf->recalc_ac_an(1); $x = $vcf->next_data_hash(); print $vcf->format_line($x);
    Args 1  : 0 .. never recalculate
              1 .. recalculate if present
              2 .. recalculate if present and add if missing

=cut

sub recalc_ac_an
{
    my ($self,$value) = @_;
    if ( $value eq '0' || $value eq '1' || $value eq '2' ) { $$self{recalc_ac_an} = $value; }
    return;
}



sub _format_line_hash
{
    my ($self,$record,$columns) = @_;

    my $ngtypes = scalar keys %{$$record{gtypes}};
    if ( !$$self{columns} ) { $self->_fake_column_names(9 + $ngtypes); }
    my $cols = $$self{columns};

    # CHROM  POS     ID      REF
    my $out;
    for my $i (0..3) { $out .= $$record{$$cols[$i]} ."\t";  }

    # ALT
    $out .= join(',',@{$$record{$$cols[4]}});

    # QUAL
    $out .= "\t". $$record{$$cols[5]};

    # FILTER
    $out .= "\t". join(';',@{$$record{$$cols[6]}});

    # Collect the gtypes of interest
    my $gtypes;
    if ( $columns )
    {
        # Select only those gtypes keys with a corresponding key in columns.
        for my $col (@$columns) { $$gtypes{$col} = $$record{gtypes}{$col}; }
    }
    else
    {
        $gtypes = $$record{gtypes};
    }

    # INFO
    # .. calculate AN and AC, but only if recalc_ac_an is set
    my $needs_an_ac;
    my @info;
    while (my ($key,$value) = each %{$$record{$$cols[7]}})
    {
        if ( $$self{recalc_ac_an}>0 )
        {
            if ( $key eq 'AN' ) { $needs_an_ac=1; next; }
            if ( $key eq 'AC' ) { $needs_an_ac=1; next; }
        }
        push @info, (defined $value ? "$key=$value" : $key);
    }
    if ( $needs_an_ac )
    {
        my ($an,$ac) = $self->calc_an_ac($gtypes);
        push @info, "AN=$an","AC=$ac";
    }
    $out .= "\t". join(';', sort @info);

    # FORMAT
    $out .= "\t". join(':',@{$$record{$$cols[8]}});

    # genotypes
    if ( $columns )
    {
        for my $col (@$columns)
        {
            my $gt = $$gtypes{$col};
            $out .= "\t" . join(':', map { exists($$gt{$_}) ? $$gt{$_} : '' } @{$$record{FORMAT}});
        }
    }
    else
    {
        for my $i (1..$ngtypes)
        {
            my $gt = $$gtypes{$$cols[$i+8]};
            # Get all the fields specified in FORMAT. If not available, print empty string instead.
            $out .= "\t" . join(':', map { exists($$gt{$_}) ? $$gt{$_} : '' } @{$$record{FORMAT}});
        }
    }

    $out .= "\n";
    return $out;

}

sub calc_an_ac
{
    my ($self,$gtypes) = @_;
    my ($an,%ac_counts);
    $an = 0;
    for my $gt (keys %$gtypes)
    {
        my $value = $$gtypes{$gt}{GT};
        if ( $value eq '.' || $value eq './.' ) { next; }
        my ($al1,$al2) = split(m{[\|/]},$value);
        if ( defined($al1) )
        {
            $an++;
            if ( $al1 ne '0' ) { $ac_counts{$al1}++; }
        }
        if ( defined($al2) )
        {
            $an++;
            if ( $al2 ne '0' ) { $ac_counts{$al2}++; }
        }
    }
    my @ac;
    for my $ac ( sort { $a <=> $b } values %ac_counts) { push @ac, $ac; }
    if ( !@ac ) { @ac = ('0'); }
    return ($an,join(',',@ac));
}


=head2 validate_alt_field

    Usage   : my $x = $vcf->next_data_hash(); $vcf->validate_alt_field($$x{ALT});
    Args    : The ALT arrayref
    Returns : Error message in case of an error.

=cut

sub validate_alt_field
{
    my ($self,$values) = @_;

    if ( @$values == 1 && $$values[0] eq '.' ) { return undef; }
    
    my @err;
    for my $item (@$values)
    {
        if ( $item=~/^[ACTGN]$/ ) { next; }
        elsif ( $item=~/^I[ACTGN]+$/ ) { next; }
        elsif ( $item=~/^D\d+$/ ) { next; }

        push @err, $item;
    }
    if ( !@err ) { return undef; }
    return 'Could not parse the allele(s) [' .join(',',@err). ']';
}


=head2 validate_filter_field

    Usage   : my $x = $vcf->next_data_hash(); $vcf->validate_filter_field($$x{FILTER});
    Args    : The FILTER arrayref
    Returns : Error message in case of an error.

=cut

sub validate_filter_field
{
    my ($self,$values) = @_;

    if ( @$values == 1 && $$values[0] eq '.' ) { return undef; }
    
    my @errs;
    for my $item (@$values)
    {
        if ( $item eq '0' ) { next; }
        if ( exists($$self{header}{FILTER}{$item}) ) { next; }
        push @errs, $item;
        $self->_add_filter_field("$item,No description");
    }
    if ( !@errs ) { return undef; }
    return 'The filter(s) [' . join(',',@errs) . '] not listed in the header.';
}


=head2 validate_info_field

    Usage   : my $x = $vcf->next_data_hash(); $vcf->validate_info_field($$x{INFO});
    Args    : The INFO hashref
    Returns : Error message in case of an error.

=cut

sub validate_info_field
{
    my ($self,$values) = @_;

    my @errs;
    while (my ($key,$value) = each %$values)
    {
        if ( !exists($$self{header}{INFO}{$key}) )
        {
            push @errs, "INFO tag [$key] not listed in the header";
            my $nargs = defined $value ? -1 : 0;
            $self->_add_field('INFO',qq[$key,$nargs,String,"No description"]);
            next;
        }
        my $type = $$self{header}{INFO}{$key};
        if ( $$type{nparams}==0 ) 
        { 
            if ( defined($value) ) { push @errs, "INFO tag [$key] did not expect any parameters, got [$value]"; }
            next; 
        }

        my @vals = split(/,/, $value);
        if ( $$type{nparams}!=-1 && @vals!=$$type{nparams} )
        {
            push @errs, "INFO tag [$key=$value] expected different number of values ($$type{nparams})";
        }
        if ( !$$type{handler} ) { next; }
        for my $val (@vals)
        {
            my $err = &{$$type{handler}}($self,$val,$$type{missing});
            if ( $err ) { push @errs, $err; }
        }
    }
    if ( !@errs ) { return undef; }
    return join(',',@errs);
}



=head2 validate_gtype_field

    Usage   : my $x = $vcf->next_data_hash(); $vcf->validate_gtype_field($$x{FORMAT},$$x{gtypes}{NA00001});
    Args    : The genotype data hashref
              The ALT arrayref
    Returns : Error message in case of an error.

=cut

sub validate_gtype_field
{
    my ($self,$data,$alts) = @_;

    my @errs;
    while (my ($key,$value) = each %$data)
    {
        if ( !exists($$self{header}{FORMAT}{$key}) )
        {
            push @errs, "FORMAT tag [$key] not listed in the header";
            $self->_add_field('FORMAT',qq[$key,-1,String,"No description"]);
            next;
        }
        my $type = $$self{header}{FORMAT}{$key};

        my @vals = split(/,/, $value);
        if ( $$type{nparams}!=-1 && @vals!=$$type{nparams} )
        {
            push @errs, "FORMAT tag [$key=$value] expected different number of values ($$type{nparams})";
        }
        if ( !$$type{handler} ) { next; }
        for my $val (@vals)
        {
            my $err = &{$$type{handler}}($self,$val,$$type{missing});
            if ( $err ) { push @errs, $err; }
        }
    }
    if ( !exists($$data{GT}) ) { push @errs, "The mandatory tag GT not present."; }
    elsif ( !($$data{GT} =~ m{^(\.|\d+)(?:[\|/](\.|\d+))?$}) ) { push @errs, "Unable to parse the GT field [$$data{GT}]."; } 
    else
    {
        my $nalts = @$alts==1 && $$alts[0] eq '.' ? 0 : @$alts;

        my $a = $1;
        my $b = $2;
        my $err = $self->validate_int($a,'.');
        if ( $err ) { push @errs,$err; }
        elsif ( $a<0 || $a>$nalts ) { push @errs, "Bad ALT value in the GT field [$$data{GT}]."; }
        if ( defined($b) )
        {
            $err = $self->validate_int($b,'.');
            if ( $err ) { push @errs,$err; }
            elsif ( $b<0 || $b>$nalts ) { push @errs, "Bad ALT value in the GT field [$$data{GT}]."; }
        }
    }
    if ( !@errs ) { return undef; }
    return join(',',@errs);
}


sub validate_int
{
    my ($self,$value,$default) = @_;

    if ( defined($default) && $value eq $default ) { return undef; }
    if ( $value =~ /^-?\d+$/ ) { return undef; }
    return "Could not validate the int [$value]";
}

sub validate_float
{
    my ($self,$value,$default) = @_;
    if ( defined($default) && $value eq $default ) { return undef; }
    if ( $value =~ /^-?\d*(?:.\d*)$/ ) { return undef; }
    return "Could not validate the float [$value]";
}

sub validate_char
{
    my ($self,$value,$default) = @_;

    if ( defined($default) && $value eq $default ) { return undef; }
    if ( length($value)==1) { return undef; }
    return "Could not validate the char value [$value]";
}


1;

