#!/usr/bin/env perl
#
# Take a list of bams and stream an interval merged sam out.
#
# bam_list_streamer.pl my.bam.list chr:start-end | something_expecting_sam
#

use strict;
use warnings;
use File::Spec;
use File::Basename;
use File::Copy;

# where are the picard jar files?
my $picard_dir = $ENV{PICARD} || die "I need the PICARD environment variable to point to the directory containing picard jar files\n";
my $jar = File::Spec->catfile($picard_dir, 'MergeSamFiles.jar');
-e $jar || die "Could not find MergeSamFiles.jar in $picard_dir\n";

# deal with user args
my $bam_list = shift;
my $interval = shift;
unless ($bam_list && -e $bam_list && $interval) {
    help();
}

my ($chr, $start, $end) = $interval =~ /^([^:\s]+):(\d+)\-(\d+)/;
unless ($chr && $start && $end) {
    help();
}
chomp($interval);

sub help {
    die "Usage: $0 my.bam.list chr:start-end | something_expecting_sam\n";
}

# get the list of bams
my @bams;
open(my $list_fh, $bam_list) || die "Could not open $bam_list\n";
while (<$list_fh>) {
    chomp;
    my $bam = $_;
    if (-s $bam && -s "$bam.bai") {
        push(@bams, $bam);
    }
    else {
        warn "The bam '$bam' in $bam_list didn't seem to exist or had no index; ignored\n";
    }
}

# merge the header once using picard; other concurrent jobs should not also
# attempt this. This section probably not very robust and should be improved...
my $merged_header = $bam_list.'.merged_header';
unless (-e $merged_header) {
    my $sentinal = $bam_list.'.sentinal';
    sleep(int(rand(20)) + 5);
    
    unless (-e $sentinal) {
        # use the sentinal file to stop other jobs doing this
        sub INT_handler {
            unlink($sentinal);
            exit(0);
        }
        $SIG{'INT'} = 'INT_handler';
        open(my $sfh, '>', $sentinal) || die "Could not write to $sentinal\n";
        
        # grab all the headers with samtools view and output to temp sam files,
        # whilst contructing the picard command line
        my $count = 0;
        my @bam_header_files;
        my $picard_command = 'java -Xmx1000m -jar '.$jar.' VALIDATION_STRINGENCY=SILENT O='.$merged_header.'.creating.sam ';
        foreach my $bam (@bams) {
            $count++;
            my $header_file = "$bam_list.bam$count.header.sam";
            my $samtools_command = "samtools view -H $bam > $header_file";
            system($samtools_command) && die "[$samtools_command] failed: $!\n";
            push(@bam_header_files, $header_file);
            $picard_command .= "I=$header_file ";
        }
        
        # do the merge
        system($picard_command) && die "header merge [$picard_command] failed: $!\n";
        
        # cleanup
        close($sfh);
        unlink($sentinal);
        foreach my $header (@bam_header_files) {
            unlink($header);
        }
        move("$merged_header.creating.sam", $merged_header) || die "Could not move $merged_header.creating.sam to $merged_header\n";
    }
    else {
        # wait up to 500 seconds for another job to finish merging the header
        my $attempts = 0;
        while ($attempts < 50) {
            last if -e $merged_header;
            warn "Waiting for another job to create $merged_header...\n";
            sleep(10);
            $attempts++;
        }
        unless (-e $merged_header) {
            die "Waited for another job to create $merged_header, but it never happened; giving up!\n";
        }
    }
}

# construct a big pipe command line to samtools view all the bams in the
# desired interval, going through unix sort and prepending the header before
# outputting
my $pipe_command = '(';
my @views;
foreach my $bam (@bams) {
    push(@views, "samtools view $bam $interval");
}
$pipe_command .= join('; ', @views);
$pipe_command .= qq[) | sort -k4,4n | perl -e 'open(\$fh, "$merged_header") || die "Could not open $merged_header\n"; while (<\$fh>) { print; } while (<>) { print }' | samtools view -Sbu - ];

# stream away
system($pipe_command);

exit;
