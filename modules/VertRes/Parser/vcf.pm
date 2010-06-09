=head1 NAME

VertRes::Parser::vcf - parse vcf files

=head1 SYNOPSIS

use VertRes::Parser::vcf;

# create object, supplying vcf file or filehandle
my $pars = VertRes::Parser::vcf->new(file => 'my.vcf');

# get the hash reference that will hold the most recently requested result
my $result_holder = $pars->result_holder();

# loop through the vcf, getting results
while ($pars->next_result()) {
    my $chrom = $result_holder->{CHROM};
    my $dp = $result_holder->{INFO}->{DP};
}

# only get results that match a given filter:
$pars->match_filter(0);
while ($pars->next_result()) {
    my %samples_hash = %{$result_holder->{SAMPLES}};
    while (my ($sample_name, $data_hash) = each %samples_hash) {
        my $unfiltered_sample_gt = $data_hash->{GT};
    }
}

=head1 DESCRIPTION

A parser for VCF files. Currently only supports VCF 3.3, which must be stated
in the first line of the file.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Parser::vcf;

use strict;
use warnings;

use base qw(VertRes::Parser::ParserI);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Parser::vcf->new(file => 'filename');
 Function: Build a new VertRes::Parser::vcf object.
 Returns : VertRes::Parser::vcf object
 Args    : file => filename -or- fh => filehandle

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    # unlike normal parsers, our result holder is a hash ref
    $self->{_result_holder} = {};
    
    return $self;
}

sub _get_header {
    my $self = shift;
    
    my $fh = $self->fh() || return;
    
    return 1 if $self->_header_parsed;
    
    my $ok = 0;
    while (<$fh>) {
        if (/^##/) {
            if (/^##fileformat=VCFv(\d+.\d+)/i) {
                $self->{vcf_version} = $1;
            }
        }
        elsif (/^#/) {
            $self->_set_header_parsed();
            $ok = 1;
            
            s/^#//;
            my @cols = split;
            my %col_to_name;
            my $filter_col = 6;
            foreach my $col (0..$#cols) {
                $col_to_name{$col} = $cols[$col];
                if ($cols[$col] eq 'FILTER') {
                    $filter_col = $col;
                }
            }
            $self->{col_to_name} = \%col_to_name;
            $self->{sorted_cols} = [sort { $a <=> $b } keys %col_to_name];
            $self->{filter_col} = $filter_col;
            
            last;
        }
        else {
            my $version = $self->{vcf_version} || '?';
            $self->warn("Not a properly formatted VCF v$self->{vcf_version} file");
            last;
        }
    }
    
    $self->{vcf_version} || $self->warn("There was no VCF version line in the header, can't cope!");
    unless ($self->{vcf_version} == 3.3) {
        $self->warn("Parser can only cope with VCF v3.3 atm...");
    }
    
    return $ok;
}

=head2 result_holder

 Title   : result_holder
 Usage   : my $result_holder = $obj->result_holder()
 Function: Get the data structure that will hold the last result requested by
           next_result()
 Returns : hash ref, where the keys are:
           CHROM -> string
           POS -> int
           ID -> string
           REF -> string
           ALT -> array ref of alternate alleles
           QUAL -> number
           FILTER -> string
           INFO -> hash ref with keys as INFO field names (eg. 'DP')
           FORMAT -> array ref of format field names

           There is also a special 'SAMPLES' key that has a hash ref for a
           value. This hash ref has sample ids as keys, and another hash ref
           as the value. These hash refs have format field names (eg. 'GT') as
           keys and the appropriate results as values.

           When values are '.' in the VCF file, meaning 'no value', then undef
           is returned instead of the string '.'
 Args    : n/a

=cut

=head2 match_filter

 Title   : match_filter
 Usage   : my $obj->match_filter(0);
 Function: Change the behaviour of next_result(), so that it only returns
           results that match the given filter.
 Returns : string (the filter set)
 Args    : string (the filter to set; not checked for validity)

=cut

sub match_filter {
    my $self = shift;
    
    if (@_) {
        $self->{match_filter} = shift;
    }
    
    return $self->{match_filter};
}

=head2 next_result

 Title   : next_result
 Usage   : while ($obj->next_result()) { # look in result_holder }
 Function: Parse the next entry from the vcf file.
 Returns : boolean (false at end of output; check the result_holder for the
           actual result information)
 Args    : n/a

=cut

sub next_result {
    my $self = shift;
    my $no_index = shift;
    
    # get the next line
    my $fh = $self->fh() || return;
    my $fh_id = $self->_fh_id;
    
    # make sure we've gotten our header first
    $self->_get_header() || $self->throw("Unable to parse header before first result - is this really a VCF file?");
    my $line = <$fh> || return;
    
    my @data = split(qr/\t/, $line);
    @data || return;
    chomp($data[-1]);
    
    # clear data first, incase we don't overwrite all fields
    my $result_holder = $self->{_result_holder};
    foreach my $key (keys %{$result_holder}) {
        delete $result_holder->{$key};
    }
    
    # should we skip this line?
    my $match_filter = $self->match_filter;
    if (defined $match_filter) {
        my $filter_col = $self->{filter_col};
        my $this_filter = $data[$filter_col];
        while ("$this_filter" ne "$match_filter") {
            $line = <$fh> || return;
            @data = split(qr/\t/, $line);
            @data || return;
            $this_filter = $data[$filter_col];
        }
        chomp($data[-1]);
    }
    
    # fill in the new data
    my %col_to_name = %{$self->{col_to_name}};
    my @format_field_names;
    foreach my $col (@{$self->{sorted_cols}}) {
        my $col_name = $col_to_name{$col};
        my $value = $data[$col];
        
        if ($col_name eq 'CHROM' || $col_name eq 'POS' || $col_name eq 'ID' || $col_name eq 'REF' || $col_name eq 'QUAL' || $col_name eq 'FILTER') {
            $result_holder->{$col_name} = $value eq '.' ? undef : $value;
        }
        elsif ($col_name eq 'ALT') {
            $result_holder->{$col_name} = [$value eq '.' ? () : (split(',', $value))];
        }
        elsif ($col_name eq 'INFO') {
            my %info;
            # NS=3;DP=14;AF=0.5;DB;H2
            foreach my $field (split(';', $value)) {
                my ($name, $sub_value) = split('=', $field);
                unless (defined $sub_value) {
                    $sub_value = 1;
                }
                $info{$name} = $sub_value;
            }
            $result_holder->{$col_name} = \%info;
        }
        elsif ($col_name eq 'FORMAT') {
            @format_field_names = split(':', $value);
            $result_holder->{$col_name} = \@format_field_names;
        }
        else {
            # a sample column; we should already have seen the FORMAT column
            my %data;
            my @data = split(':', $value);
            @data == @format_field_names || $self->throw("number of fields for sample $col_name did not match the format specification ([@format_field_names] vs [@data])");
            foreach my $i (0..$#format_field_names) {
                $data{$format_field_names[$i]} = $data[$i];
            }
            $result_holder->{SAMPLES}->{$col_name} = \%data;
        }
    }
    
    return 1;
}

1;
