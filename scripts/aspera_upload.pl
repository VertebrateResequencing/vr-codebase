#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Carp;
use Net::FTP;
use VertRes::Wrapper::iRODS;
     
$|++;

my ($fofn, $help, $password, $server, $username, $ftp, $archive,$g1k);

GetOptions(
    'f|fofn=s'      => \$fofn,
    's|server=s'    => \$server,
    'p|password=s'  => \$password,
    'u|username=s'  => \$username,
    'h|help'        => \$help,
    'a|archive'     => \$archive,
    'ftp'           => \$ftp,
    'g1k'           => \$g1k,
    );

(-s $fofn && !$help) or die <<USAGE;
    Usage: $0   
                --fofn <file of filenames that you want to upload>
		--server <location of the aspera server [default: fasp.era.ebi.ac.uk] >
		--username <username to use for connecting [default: g1k-drop-si]>
		--password <password to use [6xIktskp]>
                --ftp <use ncftp instead of aspera>
                --archive <retrieve files from fuse or irods, depending on type>
                --g1k <default to g1k username & password>
                --help <this message>

Script to automate ERA upload of a set of files.  

Checks filesize of files on ERA server before upload to see if it is necessary,
and after to see if upload was successful.  Automatically retries upload 2
times if file sizes don't match.

USAGE


# Need to check file sizes via ftp, regardless of whether using aspera or ftp
# for upload
my $ftp_server = 'ebi-ftp.internal.sanger.ac.uk'; # 10GigE internal link 
#my $ftp_server = 'ftp.sra.ebi.ac.uk';   # default ftp server


# defaults
if ($ftp){
    if ($server){
        $ftp_server = $server;  # check using passed server as well
    }
    else {
        $server = $ftp_server;  # use default ftp server
    }
}
else {
    $server ||= 'fasp.sra.ebi.ac.uk';
}

if ($g1k){
    $username ||= 'era-drop-2';
    $password ||= 'nhjkR59p';
}
else {
    $username ||= 'g1k-drop-si';
    $password ||= 'Rn3znNB5';
}

open (my $fh, $fofn) or die "Can't open $fofn: $!\n";

while (<$fh>){
    chomp;
    my $file = $_;
    next unless $file;
    my $tmpfile;
    if ($archive){
        if ($file =~ /.srf$/){
            $file = "/fuse/mpsafs/all/$file";
        }
        elsif ($file =~ /.bam$/){
            # Have to do the match_filesize check first here, because getting
            # the file is expensive
            my $irodsfile = find_irods_file($file);
            die "Can't find $file\n" unless $irodsfile;
            if (match_irods_filesize($irodsfile, $file)){
                print "$file already uploaded with correct filesize.  Skipping\n";
                next;
            }
            $file = get_irods_file($irodsfile);
            unless ($file){
                carp "Problem retrieving irods file : $irodsfile, skipping\n";
                next;
            }
            $tmpfile = $file;
        }
        else {
            die "Filetype for $file not recognised\n";
        }
    }

    if( ! -f $file ){
        carp "Cant find file: $file";
        next;
    }

    my ($filename, $path) = fileparse($file);

    # check file isn't there already
    if (match_filesize($file, $filename)){
        print "$file already uploaded with correct filesize.  Skipping\n";
        next;
    }
    
    # Upload, and retry max twice.
    my $max_try = 3;

    for (my $try = 1; $try <= $max_try; ++$try){
        upload_file($file, $filename, $ftp);
        if (match_filesize($file, $filename)){
            print "$file uploaded successfully\n";
            last;
        }
        else {
            if ($try < $max_try){
                print "$file uploaded with wrong filesize.  Retrying ($try)\n";
            }
            else {
                print "$file upload failed\n";
            }
        }
    }
    if ($tmpfile){
        unlink $tmpfile;
    }
}


sub upload_file {
    my ($source, $dest, $ftp) = @_;
    if ($ftp){
        upload_ftp($source, $dest);
    }
    else {
        upload_aspera($source);
    }
}

sub match_filesize {
    my ($source, $dest) = @_;
    my $orig_size = -s $source;
    my $dest_size = get_dest_filesize($dest);
    my $match = $orig_size == $dest_size ? 1 : 0;
 
    return $match;
}

sub match_irods_filesize {
    my ($irodsfile, $dest) = @_;
    my $irods = VertRes::Wrapper::iRODS->new();
    my $orig_size = $irods->get_file_size($irodsfile);
    my $dest_size = get_dest_filesize($dest);
    my $match = $orig_size == $dest_size ? 1 : 0;
    print "$irodsfile $orig_size $dest_size\n";
 
    return $match;
}

sub get_dest_filesize {
    my ($dest) = @_;
    my $ftp = Net::FTP->new($ftp_server, Debug => 0);
    $ftp->login($username, $password) or die "Cannot login ", $ftp->message;

    my $dest_size = $ftp->size($dest) || 0;
    $ftp->quit();

    return $dest_size;
}

sub upload_ftp {
    my ($source, $dest) = @_;

    # -z auto-resume (don't copy already copied files)
    # -C take local path of file and put to remote file
    # ncftp resets the last modified time of the uploaded file to that of the
    # original file ERA don't want that, as they want to clear the dropboxes by
    # timestamp.  -X lets you send raw ftp after each file transfer, and we can
    # use the MDTM command (that this version of ncftp uses) to reset the lmt
    # to "now".
    # -S uploads to temporary file with suffix, moving file when upload complete
    # Note that if used with resume, this will re-copy already copied files, as
    # it checks for the suffix filename, not the complete filename


    my $uploadtimestamp = `date +%Y%m%d%H%M%S`;
    chomp $uploadtimestamp;
    my $command = qq(ncftpput -u "$username" -p "$password" -z -X "MDTM $uploadtimestamp $dest" -S ".tmp" -C $server $source $dest);

    unless (system($command) == 0){
        # carp "$command failed: ",($? & 127);
        #if ($? == -1) {
        #    carp "failed to execute: $!\n";
        #}
        #elsif ($? & 127) {
        #    carp sprintf "ncftp $file died with signal %d\n", ($? & 127);
        #}
        #else {
        #    carp printf "ncftp $file exited with value %d\n", $? >> 8;
        #}
    }

}


sub upload_aspera {
    my ($source) = @_;

    # add -L- for noisy logging version
    my $command = "|ascp -k 1 -QT -m 50M -l 400M $source $username\@$server:";
    open(P,$command) || confess("Failed to open ".$command);
    print P "$password\n";
    close P;
}


sub find_irods_file {
    my ($file) = @_;
    my $irods = VertRes::Wrapper::iRODS->new();
    $file =~ /^(\d+)_(\d+)/;
    my ($run, $pos) = ($1,$2);
    my @bam = @{$irods->find_files_by_run_lane($run,$pos)};
    my $irodsloc;
    if (scalar @bam > 1){
        die "More than one bam for $file\n";
    }
    elsif (scalar @bam  == 1){
        $irodsloc = $bam[0];
    }
    return $irodsloc;
}

sub get_irods_file {
    my ($irodsfile) = @_;
    my $irods = VertRes::Wrapper::iRODS->new();
    my $dest_file = $irodsfile;
    $dest_file =~ s|.*/||;
    $dest_file ="/tmp/$dest_file";

    my $retval = $irods->get_file($irodsfile,$dest_file);
    unless ($retval ==0){
        warn "iRODS get $irodsfile failed with error code: $retval\n";
        return undef;
    }

    return $dest_file;
}
