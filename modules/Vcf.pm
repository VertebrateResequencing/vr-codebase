package Vcf;

# http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:variant_call_format
# http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf4.0
# http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf3.3
# http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcfv3.2
#
# Authors: pd3
# for VCF v3.2, v3.3, v4.0

=head1 NAME

Vcf.pm.  Module for validation, parsing and creating VCF files. 
         Supported versions: 3.2, 3.3, 4.0

=head1 SYNOPSIS

From the command line:
    perl -MVcf -e validate example.vcf
    perl -I/path/to/the/module/ -MVcf -e validate_v32 example.vcf

From a script:
    use Vcf;

    my $vcf = Vcf->new(file=>'example.vcf.gz');
    $vcf->parse_header();

    # Do some simple parsing. Most thorough but slowest way how to get the data.
    while (my $x=$vcf->next_data_hash()) 
    { 
        for my $gt (keys %{$$x{gtypes}})
        {
            my ($al1,$sep,$al2) = $vcf->parse_alleles($x,$gt);
            print "\t$gt: $al1$sep$al2\n";
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

    # Only the columns NA00001, NA00002 and NA00003 will be printed.
    my @columns = qw(NA00001 NA00002 NA00003);
    print $vcf->format_header(\@columns);
    while (my $x=$vcf->next_data_array())
    {
        # this will recalculate AC and AN counts, unless $vcf->recalc_ac_an was set to 0
        print $vcf->format_line($x,\@columns); 
    }

    $vcf->close();

=cut


use strict;
use warnings;
use Carp;
use Exporter;
use Data::Dumper;

use vars qw/@ISA @EXPORT/;
@ISA = qw/Exporter/;
@EXPORT = qw/validate validate_v32/;

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

    $vcf->run_validation();
}


=head2 validate_v32

    About   : Same as validate, but assumes v3.2 VCF version.
    Usage   : perl -MVcf -e validate_v32 example.vcf.gz     # (from the command line)
    Args    : File name or file handle. When no argument given, the first command line
              argument is interpreted as the file name.

=cut

sub validate_v32
{
    my ($fh) = @_;

    if ( !$fh && @ARGV && -e $ARGV[0] ) { $fh = $ARGV[0]; }

    my %params = ( version=>'3.2' );

    my $vcf;
    if ( $fh ) { $vcf = fileno($fh) ? Vcf->new(fh=>$fh, %params) : Vcf->new(file=>$fh, %params); }
    else { $vcf = Vcf->new(%params); }

    $vcf->run_validation();
}


=head2 new

    About   : Creates new VCF reader/writer. 
    Usage   : my $vcf = Vcf->new(file=>'my.vcf', version=>'3.2');
    Args    : 
                file    .. The file name. If not given, STDIN is assumed.
                silent  .. Unless set to 0, warning messages may be printed.
                strict  .. Unless set to 0, the reader will die when the file violates the specification.
                version .. If not given, '4.0' is assumed. The header information overrides this setting.

=cut

sub new
{
    my ($class,@args) = @_;
    my $self = {@args};
    bless $self, ref($class) || $class;
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
    $$self{silent}    = 0 unless exists($$self{silent});
    $$self{strict}    = 0 unless exists($$self{strict});
    $$self{buffer}    = [];       # buffer stores the lines in the reverse order
    $$self{columns}   = undef;    # column names 
    $$self{mandatory} = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] unless exists($$self{mandatory}); 
    $$self{recalc_ac_an} = 1;
    $$self{has_header} = 0;
    $$self{default_version} = '4.0';

    $self->_set_version();
    return $self;
}

sub _set_version
{
    my ($self) = @_;
    if ( !$$self{version} ) { $$self{version}=$$self{default_version}; }

    if ( $$self{version} eq '3.2' ) { Vcf3_2->renew($self); }
    elsif ( $$self{version} eq '3.3' ) { Vcf3_3->renew($self); } 
    elsif ( $$self{version} eq '4.0' ) { Vcf4_0->renew($self); }
    else
    {
        $self->warn(qq[The version "$$self{version}" not supported, assuming VCFv4.0\n]);
    }
}


=head2 close

    About   : Close the filehandle
    Usage   : $vcf->close();
    Args    : none

=cut

sub close
{
    my ($self) = @_;
    if ( !$$self{fh} ) { return; }
    close($$self{fh});
    delete($$self{fh});
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
    Args    : Optional line to parse

=cut

sub next_data_array
{
    my ($self,$line) = @_;
    if ( !$line )
    {
        if ( @{$$self{buffer}} ) { $line = shift(@{$$self{buffer}}); }
        else { $line = readline($$self{fh}); }
    }
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

              # Or having a VCF data line $line
              my $x = $vcf->next_data_hash($line);

    Args    : Optional line to parse.

=cut

sub next_data_hash
{
    my ($self,$line) = @_;
    if ( !$line )
    {
        if ( @{$$self{buffer}} ) { $line = shift(@{$$self{buffer}}); }
        else { $line = readline($$self{fh}); }
    }
    if ( !$line ) { return undef; }
    my @items = split(/\t/,$line);
    chomp($items[-1]);

    my $cols = $$self{columns};
    if ( !$$self{columns} ) 
    { 
        $self->_fake_column_names(scalar @items); 
        $cols = $$self{columns};
    }
    else
    {
        # Check the number of columns
        if ( scalar @items != scalar @$cols )  
        { 
            $self->warn("Different number of columns at $items[0]:$items[1] (expected ".scalar @$cols.", got ".scalar @items.")\n");
            while ( $items[-1] eq '' ) { pop(@items); }
            if ( scalar @items != scalar @$cols ) 
            {
                my @test = split(/\s+/,$line);
                if ( scalar @test == scalar @$cols ) { $self->warn("(Have spaces been used instead of tabs?)\n\n"); }
                else { $self->throw("Error not recoverable, exiting.\n"); }

                @items = @test;
            }
            else { $self->warn("(Trailing tabs?)\n\n"); }
        }
    }
    my %out;

    # Mandatory fields
    for my $i (0..3,5) { $out{$$cols[$i]} = $items[$i]; }

    # ALT field
    $out{$$cols[4]} = [ split(/,/,$items[4]) ];

    # FILTER field
    $out{$$cols[6]} = [ split(/;/,$items[6]) ];

    # Info, e.g. NS=58;DP=258;AF=0.786;DB;H2
    if ( defined $items[7] )
    {
        my %hash;
        for my $info (split(/;/,$items[7]))
        {
            my ($key,$val) = split(/=/,$info);
            if ( defined $val )
            {
                $hash{$key} = $val;
            }
            elsif ( exists($$self{header}{$$cols[7]}{$key}) )
            {
                $hash{$key} = $$self{header}{$$cols[7]}{$key}{default};
            }
            else
            {
                $hash{$key} = undef;
            }
        }
        $out{$$cols[7]} = \%hash;
    }

    # The format field may not be present. GT:GQ:DP:HQ
    my $format;
    if ( $$cols[8] || $items[8] )
    {
        $format = $out{$$cols[8]} = [ split(/:/,$items[8]) ];
        if ( !$$format[0] || $$format[0] ne 'GT' ) { $self->warn("Expected GT as the first genotype field at $items[0]:$items[1]\n"); } 
    }

    # Genotype fields
    my %gtypes;
    my $check_nformat = $$self{drop_trailings} ? 0 : 1;
    for (my $icol=9; $icol<@items; $icol++)
    {
        if ( $items[$icol] eq '' ) { $self->warn("Empty column $$cols[$icol] at $items[0]:$items[1]\n"); next; }

        my @fields = split(/:/, $items[$icol]);
        if ( $check_nformat && @fields != @$format ) 
        {
            $self->warn("Different number of fields in the format and the column $$cols[$icol] at $items[0]:$items[1] ("
                .scalar @fields." vs ".scalar @$format.": [",join(',',@fields),"] vs [",join(',',@$format),"])\n");
        }
        my %hash;
        for (my $ifield=0; $ifield<@fields; $ifield++)
        {
            $hash{$$format[$ifield]} = $fields[$ifield];
        }
        $gtypes{$$cols[$icol]} = \%hash;
    }
    $out{gtypes} = \%gtypes;
    # Was this required? Does not seem so...
    #   $out{format} = $format;

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
    confess @msg;
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

    About   : Stores the header lines and meta information, such as fields types, etc.
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

    my $rec = $self->parse_header_line($line);
    if ( $rec ) { $self->add_header_line($rec); }

    return $rec;
}

=head2 add_header_line

    Usage   : $vcf->add_header_line({key=>'INFO', ID=>'AC',Number=>-1,Type=>'Integer',Description=>'Allele count in genotypes'})
              $vcf->add_header_line({key=>'reference',value=>'1000GenomesPilot-NCBI36'})
    Args    : 
    Returns : 

=cut

sub add_header_line
{
    my ($self,$rec) = @_;

    my $key = $$rec{key};
    if ( !$key ) { $self->throw("Missing the key\n"); }

    if ( exists($$rec{Type}) )
    {
        if ( !exists($$rec{default}) )
        {
            my $type = $$rec{Type};
            if ( exists($$self{defaults}{$type}) ) { $$rec{default}=$$self{defaults}{$type}; }
            else { $$rec{default}=$$self{defaults}{default}; }
        }
        if ( !exists($$rec{handler}) )
        {
            my $type = $$rec{Type};
            if ( exists($$self{handlers}{$type}) ) { $$rec{handler}=$$self{handlers}{$type}; }
        }
    }

    if ( $key eq 'INFO' or $key eq 'FILTER' or $key eq 'FORMAT' )
    {
        my $id = $$rec{ID};
        if ( !defined $id ) { $self->throw("Missing ID for the key $key: ",Dumper($rec)); }
        if ( exists($$self{header}{$key}{$id}) ) 
        {
            $self->warn("The header tag $key:$id already exists, ignoring.\n");
            return;
        }
        $$self{header}{$key}{$id} = $rec;
        push @{$$self{header_lines}}, $rec;
        return;
    }

    if ( $key eq 'format' ) { $key='fileformat'; }
    if ( $key eq 'fileformat' )
    {
        my $value = $$rec{value};
        if ( !($value=~/(\d+(?:\.\d+)?)$/) )
        {
            $self->warn("Could not parse the fileformat version string [$value], assuming VCFv$$self{default_version}\n");
            $$self{version} = $$self{default_version};
        }
        else { $$self{version} = $1; }
        $self->_set_version();

        $$rec{key} = $key;
    }

    if ( exists($$self{header}{$key}) ) 
    {
        $self->warn("The header tag $key already exists, ignoring.\n");
        return;
    }

    $$self{header}{$key} = $rec;
    push @{$$self{header_lines}}, $rec;
}

=head2 parse_header_line

    Usage   : $vcf->parse_header_line(q[##reference=1000GenomesPilot-NCBI36])
              $vcf->parse_header_line(q[##INFO=NS,1,Integer,"Number of Samples With Data"])
    Args    : 
    Returns : 

=cut

sub parse_header_line
{
    my ($self,$line) = @_;
    
    chomp($line);
    $line =~ s/^##//;
    
    if ( !($line=~/^([^=]+)=/) ) { return { key=>$line, value=>'' }; }
    my $key   = $1;
    my $value = $';

    my $desc;
    if ( $value=~/,\"([^\"]+)\"$/ ) { $desc=$1; $value=$`; }

    if ( !$desc ) { return { key=>$key, value=>$value }; }

    if ( $key eq 'INFO' or $key eq 'FORMAT' )
    {
        my ($id,$number,$type,@rest) = split(/,/,$value);
        if ( !$type or scalar @rest ) { $self->throw("Could not parse the header line: $line\n"); }
        return { key=>$key, ID=>$id, Number=>$number, Type=>$type, Description=>$desc };
    }
    if ( $key eq 'FILTER' )
    {
        my ($id,@rest) = split(/,/,$value);
        if ( !$id or scalar @rest ) { $self->throw("Could not parse the header line: $line\n"); }
        return { key=>$key, ID=>$id, Description=>$desc };
    }
    $self->throw("Could not parse the header line: $line\n"); 
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
    $$self{column_line} = $line;

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
    my $nfields = scalar @$fields - 1;  # The FORMAT field is in fact not mandatory

    # Check the names of the mandatory columns
    if ( $ncols < $nfields ) { chomp($line); $self->warn("Missing mandatory column names. [$line].\n"); return; }

    for (my $i=0; $i<$ncols; $i++)
    {
        if ( $cols[$i]=~/^\s+/ or $cols[$i]=~/\s+$/ ) 
        {
            $self->warn("The column name contains leading/trailing spaces, removing: '$cols[$i]'\n");
            $cols[$i] =~ s/^\s+//;
            $cols[$i] =~ s/\s+$//;
        }
        if ( $i<$nfields && $cols[$i] ne $$fields[$i] ) 
        { 
            $self->warn("Expected mandatory column [$$fields[$i]], got [$cols[$i]]\n"); 
            $cols[$i] = $$fields[$i];
        }
        $$self{has_column}{$cols[$i]} = $i+1;
    }
    $$self{columns} = \@cols;
    return;
}


=head2 _fake_column_names

    About   : When no header is present, fake column names as the default mandatory ones + numbers
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
    Args    : The columns to include on output [optional]

=cut

sub format_header
{
    my ($self,$columns) = @_;

    my $out = '';
    for my $line (@{$$self{header_lines}}) { $out .= $self->format_header_line($line); }

    # This is required when using the API for writing new VCF files and the caller does not add the line explicitly
    if ( $$self{header_lines}[0]{key} ne 'fileformat' ) { $out = "##fileformat=VCFv$$self{version}\n" .$out; }
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
    for my $i (0..3) { $out .= ($$record{$$cols[$i]} ? $$record{$$cols[$i]} : '.') ."\t";  }

    # ALT
    $out .= join(',',@{$$record{$$cols[4]}} ? @{$$record{$$cols[4]}} : '.');

    # QUAL
    $out .= "\t". $$record{$$cols[5]};

    # FILTER
    $out .= "\t". join(';',$$record{$$cols[6]} ? @{$$record{$$cols[6]}} : '.');

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
    # .. calculate NS, AN and AC, but only if recalc_ac_an is set
    my $needs_an_ac = $$self{recalc_ac_an}==2 ? 1 : 0;
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
    if ( !@info ) { push @info, '.'; }
    $out .= "\t". join(';', sort @info);

    # FORMAT, the column is not required, it may not be present when there are no genotypes
    if ( exists($$cols[8]) )
    {
        $out .= "\t". join(':',@{$$record{$$cols[8]}});
    }

    # Genotypes: output all columns or only a selection?
    my @col_names = $columns ? @$columns : @$cols[9..@$cols-1];
    my $nformat = @{$$record{FORMAT}};
    for my $col (@col_names)
    {
        my $gt = $$gtypes{$col};
        my $can_drop = $$self{drop_trailings};
        my @gtype;
        for (my $i=$nformat-1; $i>=0; $i--)
        {
            my $field = $$record{FORMAT}[$i];
            if ( $i==0 ) { $can_drop=0; }

            if ( exists($$gt{$field}) ) { unshift @gtype,$$gt{$field}; $can_drop=0; }
            elsif ( $can_drop ) { next; }
            elsif ( exists($$self{header}{FORMAT}{$field}{default}) ) { unshift @gtype,$$self{header}{FORMAT}{$field}{default}; $can_drop=0; }
            else { $self->throw(qq[No value for the field "$field" and no default available, $col at $$record{CHROM}:$$record{POS}.\n]); }
        }
        $out .= "\t" . join(':',@gtype);
    }

    $out .= "\n";
    return $out;

}

sub calc_an_ac
{
    my ($self,$gtypes) = @_;
    my $sep_re = $$self{regex_gtsep};
    my ($an,%ac_counts);
    $an = 0;
    for my $gt (keys %$gtypes)
    {
        my $value = $$gtypes{$gt}{GT};
        if ( $value eq '.' || $value eq './.' ) { next; }
        my ($al1,$al2) = split($sep_re,$value);
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
    for my $ac ( sort { $a <=> $b } keys %ac_counts) { push @ac, $ac_counts{$ac}; }
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

=head2 event_type

    Usage   :   my $x = $vcf->next_data_hash(); 
                my ($al1,$sep,$al2) = $vcf->parse_alleles($x,'NA00001'); 
                my ($type,$len,$ht) = $vcf->event_type($x,$al1);
    Args    : VCF data line parsed by next_data_hash
            : Allele
    Returns :   's' for SNP and number of SNPs in the record
                'i' for indel and a positive (resp. negative) number for the length of insertion (resp. deletion)
                'r' identical to the reference, length 0
                'o' for other (complex events) and the number of affected bases

=cut

sub event_type
{
    my ($self,$rec,$allele) = @_;
    if ( exists($$rec{_cached_events}{$allele}) ) { return (@{$$rec{_cached_events}{$allele}}); }

    my ($type,$len,$ht);
    if ( $allele eq $$rec{REF} or $allele eq '.' ) { $len=0; $type='r'; $ht=$$rec{REF}; }
    elsif ( $allele=~/^[ACGT]$/ ) { $len=1; $type='s'; $ht=$allele; }
    elsif ( $allele=~/^I/ ) { $len=length($allele)-1; $type='i'; $ht=$'; }
    elsif ( $allele=~/^D(\d+)/ ) { $len=-$1; $type='i'; $ht=''; }
    else { $self->throw("Eh?: $$rec{CHROM}:$$rec{POS} .. $$rec{REF} $allele\n"); }

    $$rec{_cached_events}{$allele} = [$type,$len,$ht];
    return ($type,$len,$ht);
}


=head2 parse_alleles

    Usage   : my $x = $vcf->next_data_hash(); my ($al1,$sep,$al2) = $vcf->parse_alleles($x,'NA00001');
    Args    : VCF data line parsed by next_data_hash
            : The genotype column name
    Returns : Alleles and the separator. If only one allele is present, $sep and $al2 will be an empty string.

=cut

sub parse_alleles
{
    my ($self,$rec,$column) = @_;
    if ( !exists($$rec{gtypes}) || !exists($$rec{gtypes}{$column}) ) { $self->throw("The column not present: '$column'\n"); }

    my $gtype = $$rec{gtypes}{$column}{GT};
    if ( !($gtype=~$$self{regex_gt}) ) { $self->throw("Could not parse gtype string [$gtype]\n"); }
    my $al1 = $1;
    my $sep = $2;
    my $al2 = $3;

    if ( !$al1 ) { $al1 = $$rec{REF}; }
    elsif ( $al1 ne '.' ) 
    { 
        if ( !($al1=~/^\d+$/) ) { $self->throw("Uh, what is this? [$al1] $$rec{CHROM}:$$rec{POS}\n"); } 
        $al1 = $$rec{ALT}[$al1-1]; 
    }

    if ( !defined $al2 or $al2 eq '' )
    {
        $sep = '';
        $al2 = '';
    }
    else
    {
        if ( !$al2 ) { $al2 = $$rec{REF}; }
        elsif ( $al2 ne '.' ) { $al2 = $$rec{ALT}[$al2-1]; }
    }
    return ($al1,$sep,$al2);
}


=head2 format_genotype_strings

    Usage   : my $x = { REF=>'A', gtypes=>{'NA00001'=>'A/C'}, FORMAT=>['GT'], CHROM=>1, POS=>1, FILTER=>['.'], QUAL=>-1 };
              $vcf->format_genotype_strings($x); 
              print $vcf->format_line($x);
    Args    : VCF data line in the format as if parsed by next_data_hash with alleles written as letters.
    Returns : Modifies the ALT array and the genotypes so that ref alleles become 0 and non-ref alleles 
                numbers starting from 1.

=cut

sub format_genotype_strings
{
    my ($self,$rec) = @_;

    if ( !exists($$rec{gtypes}) ) { return; }

    my $ref = $$rec{REF};
    my $nalts = 0;
    my %alts  = ();
    my $gt_re = $$self{regex_gt2};

    for my $key (keys %{$$rec{gtypes}})
    {
        my $gtype = $$rec{gtypes}{$key}{GT};
        if ( !($gtype=~$gt_re) ) { $self->throw("Could not parse gtype string [$gtype]\n"); }
        my $al1 = $1;
        my $sep = $2;
        my $al2 = $3;

        if ( $al1 eq $ref || $al1 eq '0' || $al1 eq '*' ) { $al1 = 0; }
        else
        {
            if ( $al1=~/^\d+$/ ) { $al1 = $$rec{ALT}[$al1-1]; }

            if ( exists($alts{$al1}) ) { $al1 = $alts{$al1} }
            elsif ( $al1=~$$self{regex_snp} or $al1=~$$self{regex_ins} or $al1=~$$self{regex_del} )
            {
                $alts{$al1} = ++$nalts;
                $al1 = $nalts;
            }
            elsif ( $al1 ne '.' )
            {
                $self->throw("Could not parse the genotype string [$gtype]\n");
            }
        }

        if ( defined $al2 && $al2 ne '' )
        {
            if ( $al2 eq $ref || $al2 eq '0' || $al2 eq '*' ) { $al2 = 0; }
            else
            {
                if ( $al2=~/^\d+$/ ) { $al2 = $$rec{ALT}[$al2-1]; }

                if ( exists($alts{$al2}) ) { $al2 = $alts{$al2} }
                elsif ( $al2=~$$self{regex_snp} or $al2=~$$self{regex_ins} or $al2=~$$self{regex_del} ) 
                {
                    $alts{$al2} = ++$nalts;
                    $al2 = $nalts;
                }
                elsif ( $al2 ne '.' )
                {
                    $self->throw("Could not parse the genotype string [$gtype]\n");
                }
            }
        }
        else
        {
            $al2 = '';
            $sep = '';
        }

        $$rec{gtypes}{$key}{GT} = $al1.$sep.$al2;
    }

    $$rec{ALT} = [ sort { $alts{$a}<=>$alts{$b} } keys %alts ];
}

sub fill_ref_alt_mapping
{
    my ($self,$map) = @_;
    
    my $new_ref;
    for my $ref (keys %$map)
    {
        $new_ref = $ref;
        if ( $ref ne $new_ref ) { $self->throw("The reference prefixes do not agree: $ref vs $new_ref\n"); }
        for my $alt (keys %{$$map{$ref}})
        {
            $$map{$ref}{$alt} = $alt;
        }
    }
    $$map{$new_ref}{$new_ref} = $new_ref;
    return $new_ref;
}


=head2 format_header_line

    Usage   : $vcf->format_header_line({key=>'INFO', ID=>'AC',Number=>-1,Type=>'Integer',Description=>'Allele count in genotypes'})
    Args    : 
    Returns : 

=cut

sub format_header_line
{
    my ($self,$rec) = @_;
    my $line = "##$$rec{key}";
    $line .= "=$$rec{value}" unless !exists($$rec{value});
    $line .= "=$$rec{ID}" unless !exists($$rec{ID});
    $line .= ",$$rec{Number}" unless !exists($$rec{Number});
    $line .= ",$$rec{Type}" unless !exists($$rec{Type});
    $line .= qq[,"$$rec{Description}"] unless !exists($$rec{Description});
    $line .= "\n";
    return $line;
}

=head2 add_columns

    Usage   : $vcf->add_columns('NA001','NA0002');
    Args    : 
    Returns : 

=cut

sub add_columns
{
    my ($self,@columns) = @_;
    if ( !$$self{columns} ) 
    { 
        $$self{columns} = [ @{$$self{mandatory}} ]; 
        for my $col (@{$$self{columns}}) { $$self{has_column}{$col}=1; }
    }
    my $ncols = @{$$self{columns}};
    for my $col (@columns)
    {
        if ( $$self{has_column}{$col} ) { next; }
        $ncols++;
        push @{$$self{columns}}, $col;
    }
}

=head2 add_format_field

    Usage   : $x=$vcf->next_data_hash(); $vcf->add_format_field($x,'FOO'); $$x{gtypes}{NA0001}{FOO}='Bar'; print $vcf->format_line($x);
    Args    : The record obtained by next_data_hash
            : The field name
    Returns : 

=cut

sub add_format_field
{
    my ($self,$rec,$field) = @_;

    if ( !$$rec{FORMAT} ) { $$rec{FORMAT}=[]; }

    for my $key (@{$$rec{FORMAT}})
    {
        if ( $key eq $field ) { return; } # already there
    }
    push @{$$rec{FORMAT}}, $field;
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
    my @missing;
    for my $item (@$values)
    {
        if ( $item eq $$self{filter_passed} ) { next; }
        if ( $item=~/,/ ) { push @errs,"Expected semicolon as a separator."; }
        if ( exists($$self{header}{FILTER}{$item}) ) { next; }
        push @missing, $item;
        $self->add_header_line({key=>'FILTER',ID=>$item,Description=>'No description'});
    }
    if ( !@errs && !@missing ) { return undef; }
    if ( $$self{version}<3.3 ) { return undef; }
    return join(',',@errs) .' '. 'The filter(s) [' . join(',',@missing) . '] not listed in the header.';
}


sub _add_unknown_field
{
    my ($self,$field,$key,$nargs) = @_;
    $self->add_header_line({key=>$field,ID=>$key,Number=>$nargs,Type=>'String',Description=>'No description'});
}


=head2 validate_info_field

    Usage   : my $x = $vcf->next_data_hash(); $vcf->validate_info_field($$x{INFO});
    Args    : The INFO hashref
    Returns : Error message in case of an error.

=cut

sub validate_info_field
{
    my ($self,$values) = @_;

    if ( !defined $values ) { return 'Empty INFO field.'; }

    # First handle the empty INFO field (.)
    if ( scalar keys %$values == 1 && exists($$values{'.'}) ) { return undef; }

    my @errs;
    while (my ($key,$value) = each %$values)
    {
        if ( !exists($$self{header}{INFO}{$key}) )
        {
            push @errs, "INFO tag [$key] not listed in the header" unless $$self{version}<3.3;
            my $nargs = defined $value ? -1 : 0;
            $self->_add_unknown_field('INFO',$key,$nargs);
            next;
        }
        my $type = $$self{header}{INFO}{$key};
        if ( $$type{Number}==0 ) 
        {
            if ( defined($value) ) { push @errs, "INFO tag [$key] did not expect any parameters, got [$value]"; }
            next; 
        }

        my @vals = split(/,/, $value);
        if ( $$type{Number}!=-1 && @vals!=$$type{Number} )
        {
            push @errs, "INFO tag [$key=$value] expected different number of values ($$type{Number})";
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
    my ($self,$data,$alts,$format) = @_;

    my @errs;
    while (my ($key,$value) = each %$data)
    {
        if ( !exists($$self{header}{FORMAT}{$key}) )
        {
            push @errs, "FORMAT tag [$key] not listed in the header" unless $$self{version}<3.3;
            $self->_add_unknown_field('FORMAT',$key,-1);
            next;
        }
        my $type = $$self{header}{FORMAT}{$key};

        my @vals = split(/,/, $value);
        if ( $$type{Number}!=-1 && @vals!=$$type{Number} )
        {
            push @errs, "FORMAT tag [$key] expected different number of values ($$type{Number})";
        }
        if ( !$$type{handler} ) { next; }
        for my $val (@vals)
        {
            my $err = &{$$type{handler}}($self,$val,$$type{missing});
            if ( $err ) { push @errs, $err; }
        }
    }
    if ( !exists($$data{GT}) ) { push @errs, "The mandatory tag GT not present."; }
    elsif ( !($$data{GT} =~ $$self{regex_gt}) ) { push @errs, "Unable to parse the GT field [$$data{GT}]."; } 
    else
    {
        my $nalts = @$alts==1 && $$alts[0] eq '.' ? 0 : @$alts;

        my $a = $1;
        my $b = $3;
        my $err = $self->validate_int($a,'.');
        if ( $err ) { push @errs,$err; }
        elsif ( $a ne '.' && ($a<0 || $a>$nalts) ) { push @errs, "Bad ALT value in the GT field [$$data{GT}]."; }
        if ( $b ne '' )
        {
            $err = $self->validate_int($b,'.');
            if ( $err ) { push @errs,$err; }
            elsif ( $b ne '.' && ($b<0 || $b>$nalts) ) { push @errs, "Bad ALT value in the GT field [$$data{GT}]."; }
        }
    }
    if ( !@errs ) { return undef; }
    return join(',',@errs);
}


sub validate_ref_field
{
    my ($self,$ref) = @_;
    if ( !($ref=~/^[ACGTN]$/) ) { return "Expected one of A,C,G,T,N, got [$ref]\n"; }
    return undef;
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
    if ( $value =~ /^-?\d+(?:\.\d*)$/ ) { return undef; }
    if ( $value =~ /^-?\d*(?:\.\d+)$/ ) { return undef; }
    if ( $value =~ /^-?\d+$/ ) { return undef; }
    return "Could not validate the float [$value]";
}

sub validate_char
{
    my ($self,$value,$default) = @_;

    if ( defined($default) && $value eq $default ) { return undef; }
    if ( length($value)==1) { return undef; }
    return "Could not validate the char value [$value]";
}


=head2 run_validation

    About   : Validates the VCF file.
    Usage   : my $vcf = Vcf->new(file=>'file.vcf'); $vcf->run_validation('example.vcf.gz');
    Args    : File name or file handle.

=cut

sub run_validation
{
    my ($self) = @_;

    $self->parse_header();
    if ( !exists($$self{header}) ) { $self->warn(qq[The header not present.\n]); }
    elsif ( !exists($$self{header}{fileformat}) ) 
    {
        $self->warn(qq[The "fileformat" field not present in the header, assuming VCFv$$self{version}\n]);
    }
    elsif ( $$self{header_lines}[0]{key} ne 'fileformat' ) 
    {
        $self->warn(qq[The "fileformat" not the first line in the header\n]);
    }
    if ( !exists($$self{columns}) ) { $self->warn("No column descriptions found.\n"); }

    my $default_qual = $$self{defaults}{QUAL};
    my $warn_sorted=1;
    my ($prev_chrm,$prev_pos);
    while (my $x=$self->next_data_hash()) 
    {
        # Is the position numeric?
        if ( !($$x{POS}=~/^\d+$/) ) { $self->warn("Expected integer for the position at $$x{CHROM}:$$x{POS}\n"); }

        if ( $prev_chrm && $prev_chrm eq $$x{CHROM} && $prev_pos eq $$x{POS} )
        {
            $self->warn("Duplicate entry $$x{CHROM}:$$x{POS}\n");
        }

        # Is the file sorted?
        if ( $warn_sorted )
        {
            if ( $prev_chrm && $prev_chrm eq $$x{CHROM} && $prev_pos > $$x{POS} ) 
            { 
                $self->warn("The file is not sorted, for example $$x{CHROM}:$$x{POS} comes after $prev_chrm:$prev_pos\n");
                $warn_sorted = 0;
            }
            $prev_chrm = $$x{CHROM};
            $prev_pos  = $$x{POS};
        }

        # Is the ID non-empty?
        if ( $$x{ID}=~/^$/ ) { $self->warn("The ID should be set to . if unknown at $$x{CHROM}:$$x{POS}\n"); }

        # The reference base: one of A,C,G,T,N, non-empty.
        my $err = $self->validate_ref_field($$x{REF});
        if ( $err ) { $self->warn("$$x{CHROM}:$$x{POS} .. $err\n"); }
        
        # The ALT field (alternate non-reference base)
        $err = $self->validate_alt_field($$x{ALT},$$x{REF});
        if ( $err ) { $self->warn("$$x{CHROM}:$$x{POS} .. $err\n"); }

        # The QUAL field
        my $ret = $self->validate_float($$x{QUAL},$default_qual);
        if ( $ret ) { $self->warn("QUAL field at $$x{CHROM}:$$x{POS} .. $ret\n"); }
        elsif ( $$x{QUAL}=~/^-?\d+$/ && $$x{QUAL}<-1 ) { $self->warn("QUAL field at $$x{CHROM}:$$x{POS} is negative .. $$x{QUAL}\n"); }

        # The FILTER field
        $err = $self->validate_filter_field($$x{FILTER});
        if ( $err ) { $self->warn("FILTER field at $$x{CHROM}:$$x{POS} .. $err\n"); }

        # The INFO field
        $err = $self->validate_info_field($$x{INFO});
        if ( $err ) { $self->warn("INFO field at $$x{CHROM}:$$x{POS} .. $err\n"); } 

        while (my ($gt,$data) = each %{$$x{gtypes}})
        {
            $err = $self->validate_gtype_field($data,$$x{ALT},$$x{FORMAT});
            if ( $err ) { $self->warn("$gt column at $$x{CHROM}:$$x{POS} .. $err\n"); }
        }

        if ( exists($$x{INFO}{AN}) || exists($$x{INFO}{AC}) )
        {
            my ($an,$ac) = $self->calc_an_ac($$x{gtypes});
            if ( exists($$x{INFO}{AN}) && $an ne $$x{INFO}{AN} ) 
            { 
                $self->warn("$$x{CHROM}:$$x{POS} .. AN is $$x{INFO}{AN}, should be $an\n"); 
            }
            if ( exists($$x{INFO}{AC}) && $ac ne $$x{INFO}{AC} ) 
            { 
                $self->warn("$$x{CHROM}:$$x{POS} .. AC is $$x{INFO}{AC}, should be $ac\n"); 
            }
        }
    }
}


=head2 get_chromosomes

    About   : Get list of chromosomes from the VCF file. Must be bgzipped and tabix indexed.
    Usage   : my $vcf = Vcf->new(); $vcf->get_chromosomes();
    Args    : none

=cut

sub get_chromosomes
{
    my ($self) = @_;
    if ( !$$self{file} ) { $self->throw(qq[The parameter "file" not set.\n]); }
    my (@out) = `tabix $$self{file} -l`;
    if ( $? ) 
    { 
        $self->throw(qq[The command "tabix $$self{file} -l" exited with an error. Is the file tabix indexed?\n]); 
    }
    for (my $i=0; $i<@out; $i++) { chomp($out[$i]); }
    return \@out;
}


=head2 open_tabix

    About   : Open fh for reading from tabix
    Usage   : my $vcf = Vcf->new(); $vcf->open_tabix('1:100-1000'); while (my $line=$vcf->next_line()) { ... }
    Args    : Tabix command line argument

=cut

sub open_tabix
{
    my ($self,$arg) = @_;
    if ( !$$self{file} ) { $self->throw(qq[The parameter "file" not set.\n]); }
    $self->close();
    open($$self{fh},"tabix $$self{file} $arg 2>/dev/null|");
}


#------------------------------------------------
# Version 3.2 specific functions

package Vcf3_2;
use base qw(Vcf);

sub renew
{
    my ($class,$ref) = @_;
    my $self = $ref;
    bless $self, ref($class) || $class;

    $$self{version} = '3.2';
    $$self{drop_trailings} = 1;
    $$self{filter_passed}  = 0;

    $$self{defaults}{QUAL}    = '-1';
    $$self{defaults}{default} = '.';
    $$self{defaults}{Flag}    = undef;
    $$self{defaults}{GT}      = '.';

    $$self{handlers}{Integer}   = \&Vcf::validate_int;
    $$self{handlers}{Float}     = \&Vcf::validate_float;
    $$self{handlers}{Character} = \&Vcf::validate_char;

    $$self{regex_snp}   = qr/^[ACGTN]$/i;
    $$self{regex_ins}   = qr/^I[ACGTN]+$/;
    $$self{regex_del}   = qr/^D\d+$/;
    $$self{regex_gtsep} = qr{[\\|/]};
    $$self{regex_gt}    = qr{^(\.|\d+)([\\|/]?)(\.?|\d*)$};
    $$self{regex_gt2}   = qr{^(\.|[0-9ACGTNIDacgtn]+)([\\|/]?)([0-9ACGTNIDacgtn]*)$}; 

    return $self;
}


#------------------------------------------------
# Version 3.3 specific functions

package Vcf3_3;
use base qw(Vcf);

sub renew
{
    my ($class,$ref) = @_;
    my $self = $ref;
    bless $self, ref($class) || $class;

    $$self{version} = '3.3';
    $$self{drop_trailings} = 0;
    $$self{filter_passed}  = 0;

    $$self{defaults}{QUAL}      = '-1';
    $$self{defaults}{Integer}   = '-1';
    $$self{defaults}{Float}     = '-1';
    $$self{defaults}{Character} = '.';
    $$self{defaults}{String}    = '.';
    $$self{defaults}{Flag}      = undef;
    $$self{defaults}{GT}        = './.';
    $$self{defaults}{default}   = '.';

    $$self{handlers}{Integer}   = \&Vcf::validate_int;
    $$self{handlers}{Float}     = \&Vcf::validate_float;
    $$self{handlers}{Character} = \&Vcf::validate_char;

    $$self{regex_snp}   = qr/^[ACGTN]$/i;
    $$self{regex_ins}   = qr/^I[ACGTN]+$/;
    $$self{regex_del}   = qr/^D\d+$/;
    $$self{regex_gtsep} = qr{[\\|/]};
    $$self{regex_gt}    = qr{^(\.|\d+)([\\|/]?)(\.?|\d*)$};
    $$self{regex_gt2}   = qr{^(\.|[0-9ACGTNIDacgtn]+)([\\|/]?)((?:\.|[0-9ACGTNIDacgtn]+)?)$}; # . 0/1 0|1 A/A A|A D4/IACGT

    return $self;
}


#------------------------------------------------
# Version 4.0 specific functions

=head1 SYNOPSIS

VCFv4.0 specific functions

=cut

package Vcf4_0;
use base qw(Vcf);

sub renew
{
    my ($class,$ref) = @_;
    my $self = $ref;
    bless $self, ref($class) || $class;

    $$self{version} = '4.0';
    $$self{drop_trailings} = 1;
    $$self{filter_passed}  = 'PASS';

    $$self{defaults}{QUAL}    = '.';
    $$self{defaults}{Flag}    = undef;
    $$self{defaults}{GT}      = '.';
    $$self{defaults}{default} = '.';

    $$self{handlers}{Integer}   = \&Vcf::validate_int;
    $$self{handlers}{Float}     = \&Vcf::validate_float;
    $$self{handlers}{Character} = \&Vcf::validate_char;

    $$self{regex_snp}   = qr/^[ACGTN]$|^<[^<>\s]+>$/i;
    $$self{regex_ins}   = qr/^[ACGTN]+$/;
    $$self{regex_del}   = qr/^[ACGTN]+$/;
    $$self{regex_gtsep} = qr{[|/]};                     # | /
    $$self{regex_gt}    = qr{^(\.|\d+)([|/]?)(\d*)$};   # . 0/1 0|1
    $$self{regex_gt2}   = qr{^(\.|[0-9ACGTNacgtn]+)([|/]?)([0-9ACGTNacgtn]*)$};   # . 0/1 0|1 A/A A|A

    return $self;
}

sub Vcf4_0::format_header_line
{
    my ($self,$rec) = @_;

    my $number = exists($$rec{Number}) && $$rec{Number}==-1 ? '.' : $$rec{Number};

    my $line = "##$$rec{key}=";
    $line .= $$rec{value} unless !exists($$rec{value});
    $line .= '<' unless !exists($$rec{ID});
    $line .= "ID=$$rec{ID}" unless !exists($$rec{ID});
    $line .= ",Number=$number" unless !defined $number;
    $line .= ",Type=$$rec{Type}" unless !exists($$rec{Type});
    $line .= ",Description=\"$$rec{Description}\"" unless !exists($$rec{Description});
    $line .= ">" unless !exists($$rec{ID});
    $line .= "\n";
    return $line;
}

=head2 parse_header_line

    Usage   : $vcf->parse_header_line(q[##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">])
              $vcf->parse_header_line(q[reference=1000GenomesPilot-NCBI36])
    Args    : 
    Returns : 

=cut

sub Vcf4_0::parse_header_line
{
    my ($self,$line) = @_;
    
    chomp($line);
    $line =~ s/^##//;
    
    if ( !($line=~/^([^=]+)=/) ) { $self->throw("Expected key=value pair in the header: $line\n"); }
    my $key   = $1;
    my $value = $';

    if ( !($value=~/^<(.+)>$/) ) { return { key=>$key, value=>$value }; }

    my $rec = { key=>$key };
    my $tmp = $1;
    while ($tmp)
    {
        my ($key,$value);
        if ( $tmp=~/^([^=]+)="([^\"]+)"/ ) { $key=$1; $value=$2; }
        elsif ( $tmp=~/^([^=]+)=([^,]+)/ ) { $key=$1; $value=$2; }
        else { $self->throw(qq[Could not parse $value, got stuck here: "$tmp"\n]); }

        $$rec{$key} = $value;

        $tmp = $';
        if ( $tmp=~/^,/ ) { $tmp = $'; }
    }

    if ( !exists($$rec{ID}) ) { $self->throw("Missing the ID tag in the $value\n"); }
    if ( exists($$rec{Number}) && $$rec{Number} eq '.' ) { $$rec{Number}=-1; }
    return $rec;
}

sub Vcf4_0::validate_ref_field
{
    my ($self,$ref) = @_;
    if ( !($ref=~/^[ACGTN]+$/) ) { return "Expected combination of A,C,G,T,N, got [$ref]\n"; }
    return undef;
}

sub Vcf4_0::validate_alt_field
{
    my ($self,$values,$ref) = @_;

    if ( @$values == 1 && $$values[0] eq '.' ) { return undef; }

    my $ref_len = length($ref);
    my $ref1 = substr($ref,0,1);

    my @err;
    my $msg = '';
    for my $item (@$values)
    {
        if ( !($item=~/^[ACTGN]+$|^<[^<>\s]+>$/) ) { push @err,$item; next; }
        if ( $ref_len==1 && length($item)==1 ) { next; }
        if ( $item=~/^<[^<>\s]+>$/ ) { next; }
        if ( substr($item,0,1) ne $ref1 ) { $msg=', first base does not match the reference.'; push @err,$item; next; }
    }
    if ( !@err ) { return undef; }
    return 'Could not parse the allele(s) [' .join(',',@err). ']' . $msg;
}


=head2 fill_ref_alt_mapping

    About   : A tool for merging VCFv4.0 records. The subroutine unifies the REFs and creates a mapping
                from the original haplotypes to the haplotypes based on the new REF. Consider the following
                example:
                    REF ALT
                    G    GA
                    GT   G
                    GT   GA
                    GT   GAA
                    GTC  G
                my $map={G=>{GA=>1},GT=>{G=>1,GA=>1,GAA=>1},GTC=>{G=>}};   
                my $new_ref=$vcf->fill_ref_alt_mapping($map);
                
              The call returns GTC and $map is now
                    G    GA     ->      GTC  GATC
                    GT   G      ->      GTC  GC
                    GT   GA     ->      GTC  GAC
                    GT   GAA    ->      GTC  GAAC
                    GTC  G      ->      GTC  G
    Args    : 
    Returns : New REF string and fills the hash with appropriate ALT.

=cut

sub Vcf4_0::fill_ref_alt_mapping
{
    my ($self,$map) = @_;
    
    my $max_len = 0;
    my $new_ref;
    for my $ref (keys %$map)
    {
        my $len = length($ref);
        if ( $max_len<$len ) 
        { 
            $max_len = $len; 
            $new_ref = $ref;
        }
        $$map{$ref}{$ref} = 1;
    }
    for my $ref (keys %$map)
    {
        my $rlen = length($ref);
        if ( substr($new_ref,0,$rlen) ne $ref ) { $self->throw("The reference prefixes do not agree: $ref vs $new_ref\n"); }
        for my $alt (keys %{$$map{$ref}})
        {
            my $new = $alt;
            if ( $rlen<$max_len ) { $new .= substr($new_ref,$rlen); }
            $$map{$ref}{$alt} = $new;
        }
    }
    return $new_ref;
}



sub Vcf4_0::event_type
{
    my ($self,$rec,$allele) = @_;
    if ( exists($$rec{_cached_events}{$allele}) ) { return (@{$$rec{_cached_events}{$allele}}); }

    if ( $allele=~/^<([^>]+)>$/ ) 
    { 
        $$rec{_cached_events}{$allele} = ['u',0,$1];
        return ('u',0,$1); 
    }

    my $ref = $$rec{REF};
    my $reflen = length($ref);
    my $len = length($allele);
    
    my $ht;
    my $type;
    if ( $len==$reflen )
    {
        # This can be a reference, a SNP, or multiple SNPs
        my $mism = 0;
        for (my $i=0; $i<$len; $i++)
        {
            if ( substr($ref,$i,1) ne substr($allele,$i,1) ) { $mism++; }
        }
        if ( $mism==0 ) { $type='r'; $len=0; }
        else { $type='s'; $len=$mism; }
    }
    else
    {
        ($len,$ht)=is_indel($ref,$allele);
        if ( $len )
        {
            # Indel
            $type = 'i';
            $allele = $ht;
        }
        else 
        {
            $type = 'o'; $len = $len>$reflen ? $len-1 : $reflen-1;
        }
    }

    $$rec{_cached_events}{$allele} = [$type,$len,$allele];
    return ($type,$len,$allele);
}


sub is_indel
{
    my ($seq1,$seq2) = @_;

    my $len1 = length($seq1);
    my $len2 = length($seq2);
    if ( $len1 eq $len2 ) { return (0,''); }

    my $ht;
    my $len = $len1<$len2 ? $len1 : $len2;
    my $match  = 1;
    for (my $i=0; $i<$len; $i++)
    {
        if ( substr($seq1,$i,1) ne substr($seq2,$i,1) ) { $match=0; last; }
    }
    if ( $match ) 
    { 
        if ( $len1<$len2 ) { $ht=substr($seq2,$len); }
        else { $ht=substr($seq1,$len); }
        return ($len2-$len1,$ht); 
    }
    $match = 1;
    for (my $i=1; $i<=$len; $i++)
    {
        if ( substr($seq1,$len1-$i,1) ne substr($seq2,$len2-$i,1) ) { $match=0; last; }
    }
    if ( $match ) 
    { 
        if ( $len1<$len2 ) { $ht=substr($seq2,0,$len); }
        else { $ht=substr($seq1,0,$len); }
        return ($len2-$len1,$ht); 
    }

    ($len,$ht) = sw_align($seq1,$seq2,{same_pos=>1,match=>2,mismatch=>-100,gap=>0,is_indel=>1});
    return ($len,$ht); 
}

sub sw_align
{
    my ($seq1,$seq2,$opts) = @_;

    if ( !$opts ) { $opts={}; }

    my $MATCH    = exists($$opts{match}) ? $$opts{match} : 2;
    my $MISMATCH = exists($$opts{mismatch}) ? $$opts{mismatch} : -1;
    my $GAP      = exists($$opts{gap}) ? $$opts{gap} : -1;

    my $score = 0;
    my $ptr   = 1;
    my $none  = 0;
    my $diag  = 1;
    my $up    = 2;
    my $left  = 3;

    my $len1 = length($seq1);
    my $len2 = length($seq2);

    # initialization
    my @matrix;
    $matrix[0][0][$score] = $$opts{same_pos} ? $MATCH : 0;
    $matrix[0][0][$ptr]   = $none;
    for(my $j=1; $j<=$len1; $j++) 
    {
        $matrix[0][$j][$score] = 0;
        $matrix[0][$j][$ptr]   = $none;
    }
    for (my $i=1; $i<=$len2; $i++) 
    {
        $matrix[$i][0][$score] = 0;
        $matrix[$i][0][$ptr]   = $none;
    }

    # fill
    my $max_i     = 0;
    my $max_j     = 0;
    my $max_score = 0;
    for(my $i=1; $i<=$len2; $i++) 
    {
        for(my $j=1; $j<=$len1; $j++) 
        {
            my ($diagonal_score, $left_score, $up_score);

            # calculate match score
            my $letter1 = substr($seq1, $j-1, 1);
            my $letter2 = substr($seq2, $i-1, 1);       
            if ($letter1 eq $letter2) 
            {
                $diagonal_score = $matrix[$i-1][$j-1][$score] + $MATCH;
            }
            else 
            {
                $diagonal_score = $matrix[$i-1][$j-1][$score] + $MISMATCH;
            }

            # calculate gap scores
            $up_score   = $matrix[$i-1][$j][$score] + $GAP;
            $left_score = $matrix[$i][$j-1][$score] + $GAP;

            if ( $diagonal_score<=0 and $up_score<=0 and $left_score<=0 ) 
            {
                $matrix[$i][$j][$score] = 0;
                $matrix[$i][$j][$ptr]   = $none;
                next; 
            }

            # choose best score
            if ( $diagonal_score>=$up_score ) 
            {
                if ($diagonal_score >= $left_score) 
                {
                    $matrix[$i][$j][$score] = $diagonal_score;
                    $matrix[$i][$j][$ptr]   = $diag;
                }
                else 
                {
                    $matrix[$i][$j][$score] = $left_score;
                    $matrix[$i][$j][$ptr]   = $left;
                }
            } 
            else 
            {
                if ( $up_score >= $left_score ) 
                {
                    $matrix[$i][$j][$score] = $up_score;
                    $matrix[$i][$j][$ptr]   = $up;
                }
                else 
                {
                    $matrix[$i][$j][$score] = $left_score;
                    $matrix[$i][$j][$ptr]   = $left;
                }
            }

            # set maximum score
            if ( $matrix[$i][$j][$score] > $max_score ) 
            {
                $max_i     = $i;
                $max_j     = $j;
                $max_score = $matrix[$i][$j][$score];
            }
        }
    }

    # trace-back
    if ( $$opts{is_indel} )
    {
        # The last position must be reached in case of an indel after the checks in is_indel
        if ( $max_j!=$len1 or $max_i!=$len2 ) { return (0,''); }

        my $ht;
        my $j = $max_j;
        my $i = $max_i;
        my ($ins_j,$ins_i);
        while (1) 
        {
            last if $matrix[$i][$j][$ptr] eq $none;

            if ( $matrix[$i][$j][$ptr] eq $diag ) 
            {
                $i--; $j--;
            }
            elsif ( $matrix[$i][$j][$ptr] eq $left ) 
            {
                # Is the insertion continuous?
                if ( defined $ins_j && $ins_j!=$j+1 ) { return (0,''); }
                # Is there a deletion as well?
                if ( defined $ins_i ) { return (0,''); }
                $ht .= substr($seq1, $j-1, 1);
                $ins_j = $j;
                $j--;
            }
            elsif ( $matrix[$i][$j][$ptr] eq $up ) 
            {
                # Is the insertion continuous?
                if ( defined $ins_i && $ins_i!=$i+1 ) { return (0,''); }
                # Is there a deletion as well?
                if ( defined $ins_j ) { return (0,''); }
                $ht .= substr($seq2, $i-1, 1);
                $ins_i = $i;
                $i--;
            }   
        }
        if ( $i!=0 or $j!=0 ) { return (0,''); }
        if ( defined $ins_i && defined $ins_j ) { return (0,''); }
        if ( !defined $ins_i && !defined $ins_j ) { return (0,''); }
        my $len = length($ht);
        $ht = reverse($ht);
        if ( defined $ins_j ) { $len = -$len; }
        return ($len,$ht);
    }

    my $align1 = '';
    my $align2 = '';
    my $j = $max_j;
    my $i = $max_i;
    while (1) 
    {
        last if $matrix[$i][$j][$ptr] eq $none;

        if ( $matrix[$i][$j][$ptr] eq $diag ) 
        {
            $align1 .= substr($seq1, $j-1, 1);
            $align2 .= substr($seq2, $i-1, 1);
            $i--; $j--;
        }
        elsif ( $matrix[$i][$j][$ptr] eq $left ) 
        {
            $align1 .= substr($seq1, $j-1, 1);
            $align2 .= '-';
            $j--;
        }
        elsif ( $matrix[$i][$j][$ptr] eq $up ) 
        {
            $align1 .= '-';
            $align2 .= substr($seq2, $i-1, 1);
            $i--;
        }   
    }
    $align1 = reverse $align1;
    $align2 = reverse $align2;
    return ($align1,$align2);
}

1;

