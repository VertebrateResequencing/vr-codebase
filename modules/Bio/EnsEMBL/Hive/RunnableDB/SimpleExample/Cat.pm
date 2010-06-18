
=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Cat

=head1 SYNOPSIS


=head1 DESCRIPTION

'SimpleExample::Cat' is the final step of the SimpleExample pipeline that
splits, modifies, then cats input files.

It takes the list of split files that have been appended and cats them into
1 file.

=cut

package Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Cat;

use strict;

use base ('Bio::EnsEMBL::Hive::ProcessWithParams');

=head2 fetch_input

    Description : Implements fetch_input() interface method of
                  Bio::EnsEMBL::Hive::Process that is used to read in parameters
                  and load data. Here we get the list of append files we need
                  to cat

    param('orig_file'): the path of the original file
    param('lines_per_split'): the number of lines that the file should have
    param('string_to_append'): the string to append to each line of the split
                               files

=cut

sub fetch_input {
    my $self = shift;
    
    my $orig_file = $self->param('orig_file') || die "'orig_file' is an obligatory parameter";
    my $max_lines = $self->param('lines_per_split')  || die "'lines_per_split' is an obligatory parameter";
    my $append = $self->param('string_to_append')  || die "'string_to_append' is an obligatory parameter";
    
    my $sql = "SELECT chunk,append_file,actual_lines FROM append_files WHERE (orig_file, max_lines, append) = (?, ?, ?)";
    my $sth = $self->db->dbc()->prepare($sql);
    $sth->execute($orig_file, $max_lines, $append);
    my %append_files;
    my $total_lines = 0;
    while (my ($chunk, $append_file, $these_lines) = $sth->fetchrow_array()) {
        $append_files{$chunk} = $append_file;
        $total_lines += $these_lines;
    }
    $sth->finish();
    
    $self->param('append_files', \%append_files);
    $self->param('expected_lines', $total_lines);
}

=head2 run

    Description : Implements run() interface method of
                  Bio::EnsEMBL::Hive::Process that is used to perform the main
                  bulk of the job (minus input and output).
                  Here we do the actual concatenation of the append files.

=cut

sub run {
    my $self = shift @_;
    
    my %append_files = %{$self->param('append_files')};
    my @append_files = map { $append_files{$_} } sort { $a <=> $b } keys %append_files;
    my $orig_file = $self->param('orig_file');
    my $cat_file = $orig_file.'.result';
    
    # concatenate
    my ($in_lines, $out_lines) = (0, 0);
    open(my $ofh, '>', $cat_file) || die "Could not write to final result file '$cat_file'\n";
    foreach my $append_file (@append_files) {
        open(my $fh, $append_file) || die "Could not open append file '$append_file'\n";
        while (<$fh>) {
            $in_lines++;
            print $ofh $_;
        }
        close($fh);
    }
    close($ofh);
    
    my $expected_lines = $self->param('expected_lines');
    unless ($in_lines == $expected_lines) {
        unlink($cat_file);
        die "was only able to read in $in_lines lines from the split files, instead of the expected $expected_lines lines\n";
    }
    
    # check for truncation
    open(my $fh, $cat_file) || die "Could not read from result file '$cat_file'\n";
    while (<$fh>) {
        $out_lines++;
    }
    close($fh);
    unless ($out_lines == $expected_lines) {
        unlink($cat_file);
        die "was only able to read in $out_lines lines from the result file, instead of the expected $expected_lines lines; deleted it\n";
    }
    
    # cleanup
    foreach my $append_file (@append_files) {
        unlink($append_file);
    }
    
    $self->param('result', $cat_file);
}

=head2 write_output

    Description : Implements write_output() interface method of
                  Bio::EnsEMBL::Hive::Process that is used to deal with job's
                  output after the execution.
                  Just store path of final result file and its num of lines in
                  the db

=cut

sub write_output {
    my $self = shift @_;
    
    my $sql = "REPLACE INTO completed (orig_file, max_lines, append, result_file, actual_lines) VALUES (?, ?, ?, ?, ?) ";
    my $sth = $self->db->dbc->prepare($sql);
    $sth->execute($self->param('orig_file'), $self->param('lines_per_split'), $self->param('string_to_append'), $self->param('result'), $self->param('expected_lines'));
}

1;
