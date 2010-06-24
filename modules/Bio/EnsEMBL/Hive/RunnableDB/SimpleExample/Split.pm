
=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Split

=head1 SYNOPSIS


=head1 DESCRIPTION

'SimpleExample::Split' is the second step of the SimpleExample pipeline that
splits, modifies, then cats input files.

It takes an input file and splits into multiple chunks. For each chunk it also
creates an append job, and it creates one cat job that depends on those
completing.

=cut

package Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Split;

use strict;

use base ('Bio::EnsEMBL::Hive::ProcessWithParams');

=head2 fetch_input

    Description : Implements fetch_input() interface method of
                  Bio::EnsEMBL::Hive::Process that is used to read in parameters
                  and load data. Here we have nothing to fetch.

=cut

sub fetch_input {}

=head2 run

    Description : Implements run() interface method of
                  Bio::EnsEMBL::Hive::Process that is used to perform the main
                  bulk of the job (minus input and output).
                  Here we do the actual splitting of the input file.

    param('orig_file'): a file of filenames
    param('lines_per_split'): the max number of lines per split file
    param('string_to_append'): the string to append to each line of the split
                               files

=cut

sub run {
    my $self = shift @_;
    
    my $orig_file = $self->param('orig_file')  || die "'orig_file' is an obligatory parameter";
    my $max_lines = $self->param('lines_per_split')  || die "'lines_per_split' is an obligatory parameter";
    my $append = $self->param('string_to_append')  || die "'string_to_append' is an obligatory parameter";
    
    # split, and collect output_ids for the subsequent append jobs
    open(my $fh, $orig_file) || die "Could not open orig_file '$orig_file'\n";
    my @output_ids;
    my @split_files;
    my $lines = 0;
    my $chunk = 0;
    my $ofh;
    my $split_file;
    while (<$fh>) {
        if (! $ofh || $lines % $max_lines == 0) {
            if ($ofh) {
                close($ofh);
                push(@output_ids, { orig_file => $orig_file,
                                    split_file => $split_file,
                                    chunk => $chunk,
                                    lines_per_split => $max_lines,
                                    expected_lines => $lines,
                                    string_to_append => $append });
                push(@split_files, [$split_file, $lines]);
            }
            $chunk++;
            $split_file = "$orig_file.split_$chunk";
            open($ofh, '>', $split_file) || die "Could not write to $split_file\n";
            $lines = 0;
        }
        
        print $ofh $_;
        
        $lines++;
    }
    close($fh);
    close($ofh);
    push(@output_ids, { orig_file => $orig_file,
                        split_file => $split_file,
                        chunk => $chunk,
                        lines_per_split => $max_lines,
                        expected_lines => $lines,
                        string_to_append => $append });
    
    # check the split files aren't truncated
    foreach my $info (@split_files) {
        my ($file, $expected_lines) = @{$info};
        open(my $fh, $file) || die "Could not open split file '$file'\n";
        my $these_lines = 0;
        while (<$fh>) {
            $these_lines++;
        }
        close($fh);
        unless ($these_lines == $expected_lines) {
            unlink($file);
            die "split file '$file' was written but ended up truncated; deleted it\n";
        }
    }
    
    sleep(30);
    
    $self->param('output_ids', \@output_ids);
}

=head2 write_output

    Description : Implements write_output() interface method of
                  Bio::EnsEMBL::Hive::Process that is used to deal with job's
                  output after the execution.
                  Here we  dataflow all the append jobs into the branch-2
                  ("fan out"), and dataflow the split down branch-1 (create the
                  "funnel job").

=cut

sub write_output {
    my $self = shift @_;
    
    # some dataflow to perform:
    # "fan out" into branch-2
    $self->dataflow_output_id($self->param('output_ids'), 2);
    
    # then flow into the branch-1 funnel
    $self->dataflow_output_id($self->input_id, 1);
}

1;
