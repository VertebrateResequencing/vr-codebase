
=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Append

=head1 SYNOPSIS


=head1 DESCRIPTION

'SimpleExample::Split' is the third step of the SimpleExample pipeline that
splits, modifies, then cats input files.

It takes an split file and creates a new file with a string appended to each
line.

=cut

package Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Append;

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
                  Here we do the actual appending of a string to a split file.

    param('orig_file'): the path of the original file
    param('split_file'): the path of the split file
    param('expected_lines'): the number of lines that the file should have
    param('string_to_append'): the string to append to each line of the split
                               files
    param('chunk'): which split of the original this is

=cut

sub run {
    my $self = shift @_;
    
    my $split_file = $self->param('split_file')  || die "'split_file' is an obligatory parameter";
    my $expected_lines = $self->param('expected_lines')  || die "'expected_lines' is an obligatory parameter";
    my $append = $self->param('string_to_append')  || die "'string_to_append' is an obligatory parameter";
    
    # append and remember what we did
    open(my $fh, $split_file) || die "Could not open split file '$split_file'\n";
    my $append_file = $split_file.'.appended';
    open(my $ofh, '>', $append_file) || die "Could not write to append file '$append_file'\n";
    my ($in_lines, $out_lines) = 0;
    while (<$fh>) {
        $in_lines++;
        chomp;
        print $ofh $_, $append, "\n";
    }
    close($fh);
    close($ofh);
    
    if ($in_lines != $expected_lines) {
        unlink($append_file);
        die "split file '$split_file' was supposed to have $expected_lines lines, but was only able to read $in_lines lines\n";
    }
    
    # check for truncation
    open(my $fh, $append_file) || die "Could not open append file '$append_file'\n";
    while (<$fh>) {
        $out_lines++;
    }
    close($fh);
    if ($out_lines != $expected_lines) {
        unlink($append_file);
        die "append file '$append_file' was supposed to have $expected_lines lines, but was only able to write $out_lines lines; deleted it\n";
    }
    
    # cleanup as we go
    unlink($split_file);
    
    $self->param('result', [$append_file, $out_lines]);
}

=head2 write_output

    Description : Implements write_output() interface method of
                  Bio::EnsEMBL::Hive::Process that is used to deal with job's
                  output after the execution.
                  Here we store the names of the append files we made, and the
                  number of lines they contain, in the db.

=cut

sub write_output {
    my $self = shift @_;
    
    my $chunk = $self->param('chunk')  || die "'chunk' is an obligatory parameter";
    
    my $sql = "REPLACE INTO append_files (orig_file, max_lines, append, chunk, append_file, actual_lines) VALUES (?, ?, ?, ?, ?, ?) ";
    my $sth = $self->db->dbc->prepare($sql);
    $sth->execute($self->param('orig_file'), $self->param('lines_per_split'), $self->param('string_to_append'), $chunk, @{$self->param('result')});
}

1;
