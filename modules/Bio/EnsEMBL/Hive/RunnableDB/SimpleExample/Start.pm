
=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Start

=head1 SYNOPSIS


=head1 DESCRIPTION

'SimpleExample::Start' is the first step of the SimpleExample pipeline that
splits, modifies, then cats input files.

It creates several SimpleExample::Split jobs for each file in the input fofn.

=cut

package Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Start;

use strict;

use base ('Bio::EnsEMBL::Hive::ProcessWithParams');

=head2 fetch_input

    Description : Implements fetch_input() interface method of
                  Bio::EnsEMBL::Hive::Process that is used to read in parameters
                  and load data.
                  Here the task of fetch_input() is to parse the fofn and create
                  a set of input_ids that will be used later.

    param('input_file_list'): a file of filenames
    param('lines_per_split'): the max number of lines per split file
    param('string_to_append'): the string to append to each line of the split
                               files

=cut

sub fetch_input {
    my $self = shift @_;
    
    my $fofn = $self->param('input_file_list')  || die "'input_file_list' is an obligatory parameter";
    my $lines = $self->param('lines_per_split')  || die "'lines_per_split' is an obligatory parameter";
    my $append = $self->param('string_to_append')  || die "'string_to_append' is an obligatory parameter";
    
    # output_ids of files to be split:
    open(my $fh, $fofn) || die "Could not open input_file_list '$fofn'";
    my @output_ids;
    while (<$fh>) {
        chomp;
        push(@output_ids, { orig_file => $_,
                            lines_per_split => $lines,
                            string_to_append => $append });
    }
    
    # store them for future use:
    $self->param('output_ids', \@output_ids);
}

=head2 run

    Description : Implements run() interface method of
                  Bio::EnsEMBL::Hive::Process that is used to perform the main
                  bulk of the job (minus input and output).
                  Here we don't have any real work to do, just input and output,
                  so run() remains empty.

=cut

sub run {}

=head2 write_output

    Description : Implements write_output() interface method of
                  Bio::EnsEMBL::Hive::Process that is used to deal with job's
                  output after the execution.
                  Here we dataflow all the split jobs whose input_ids were
                  generated in fetch_input() into the branch-1 ("fan out").

=cut

sub write_output {
    my $self = shift @_;
    
    # nothing to write out, but some dataflow to perform:
    my $output_ids = $self->param('output_ids');
    
    # "fan out" into branch-1
    $self->dataflow_output_id($output_ids, 1);
}

1;
