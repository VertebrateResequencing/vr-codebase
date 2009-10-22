=head1 NAME

VertRes::Wrapper::fastqcheck - wrapper for fastqcheck

=head1 SYNOPSIS

use VertRes::Wrapper::fastqcheck;

my $wrapper = VertRes::Wrapper::fastqcheck->new();

$wrapper->run('my.fastq.gz', 'my.fastq.gz.fastqcheck');

=head1 DESCRIPTION

Runs fastqcheck in a nice way. Most importantly, lets you avoid the problem
with fastq files that have spaces in the sequence identifiers, and auto-handles
gzipped fastq input.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::fastqcheck;

use strict;
use warnings;
use VertRes::IO;

use base qw(VertRes::Wrapper::WrapperI);

=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::fastqcheck->new();
 Function: Create a VertRes::Wrapper::fastqcheck object.
 Returns : VertRes::Wrapper::fastqcheck object
 Args    : quiet => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(exe => 'fastqcheck', @args);
    
    return $self;
}

=head2 version

 Title   : version
 Usage   : my $version = $obj->version();
 Function: Returns the program version.
 Returns : undef - fastqcheck doesn't report its version
 Args    : n/a

=cut

sub version {
    my $self = shift;
    return;
}

=head2 run

 Title   : run
 Usage   : $obj->run($fastq_file, $output_file);
 Function: Fastqcheck a file.
 Returns : n/a
 Args    : paths to input and output files

=cut

sub _pre_run {
    my ($self, $fastq_file, $output_file) = @_;
    
    unless (defined $self->{_orig_exe}) {
        $self->{_orig_exe} = $self->exe;
    }
    
    # this is going to go to shell, so spaces in filename needs to be escaped
    $fastq_file =~ s/ /\\ /g;
    $self->register_output_file_to_check($output_file);
    $output_file =~ s/ /\\ /g;
    
    my $exe = "cat $fastq_file | ";
    
    # run via zcat pipe if fastq is gz
    if ($fastq_file =~ /\.gz$/) {
        $exe = 'z'.$exe;
    }
    
    # must be run via an awk pipe to avoid bug where fastqcheck can't cope with
    # spaces in sequence ids
    $exe .= "awk '{print \$1}' | ".$self->{_orig_exe};
    
    $self->exe($exe);
    
    $output_file = ' > '.$output_file;
    
    return ($output_file);
}

1;
