package Pathogens::Utils::FastqPooler;

use Moose;
use Utils;
use File::Spec;
use Pathogens::Exception;

has 'output_directory' =>  ( is => 'ro', isa => 'Str', required => 1 );
has 'lane_paths' => ( is => 'ro', isa => 'ArrayRef[Str]', required => 1);


sub _get_files_to_cat {
    my ($lane_paths) = @_;

    my @fwd_fastq;
    my @rev_fastq;

    for my $lane_path (@$lane_paths) {
        my ($base_directory, $base, $suff) = Utils::basename($lane_path);
        opendir(my $lane_dir_handle, $base_directory);
        my @fastq_files  = grep { /\.fastq\.gz$/ } readdir($lane_dir_handle);

        if(@fastq_files >=1 && -e $lane_path ."_1.fastq.gz") {
            push(@fwd_fastq, $lane_path."_1.fastq.gz");
        }
        if(@fastq_files >=2  && -e $lane_path."_2.fastq.gz") {
            push(@rev_fastq, $lane_path."_2.fastq.gz");
        }
    }

    return \@fwd_fastq, \@rev_fastq;
}


sub _cat_files {
    my ($infiles, $outfile) = @_;
    unless (scalar(@$infiles) == 0 or -e $outfile) {
        my $cmd = 'gzip -cd ' . join(' ', @{$infiles}) . " > $outfile";
        if (system($cmd)) {
            Pathogens::Exception::SystemCall->throw(error => "Error running command: $cmd");
        }
    }
}


sub run {
    my ($self) = @_;
    my ($fwd_fastqs, $rev_fastqs) = _get_files_to_cat($self->lane_paths);
    unless (-d $self->output_directory) {
        mkdir($self->output_directory);
    }
    my $fwd_out = File::Spec->catfile($self->output_directory, 'forward.fastq');
    my $rev_out = File::Spec->catfile($self->output_directory, 'reverse.fastq');
    _cat_files($fwd_fastqs, $fwd_out);
    _cat_files($rev_fastqs, $rev_out);
}


no Moose;
__PACKAGE__->meta->make_immutable;
