=head1 NAME

VertRes::IO - do basic IO on files/ filehandles

=head1 SYNOPSIS

use VertRes::IO;

# supply a file or filehandle to have it opened
my $io1 = VertRes::IO->new(file => $infile);
my $io2 = VertRes::IO->new(fh => \*STDIN);
my $io3 = VertRes::IO->new(file => ">$infile.out.gz");

# read in $infile line by line, print back out compressed
my $in_fh = $io1->fh();
my $out_fh = $io3->fh();
while (<$in_fh>) {
    print $out_fh $_;
}

# explicitly close the filehandles if you like
$io1->close();
$io3->close();

=head1 DESCRIPTION

Handles IO including reading and writing to compressed files.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::IO;

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Spec;
use IO::Uncompress::Gunzip;
use File::Fetch;
use Net::FTP::Robust;
use VertRes::Utils::FileSystem;

use base qw(VertRes::Base);

=head2 new

 Title   : new
 Usage   : my $self = $class->SUPER::new(@args);
 Function: Instantiates your object, opening user's file/ filehandle.
 Returns : $self hash-ref blessed into your class
 Args    : file => filename -or- fh => filehandle

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    # if file or fh args are supplied, those methods will get called and we'll
    # auto-open the file and set the fh
    $self->_set_from_args(\@args, methods => [qw(file fh)]);
    
    return $self;
}

=head2 file

 Title   : file
 Usage   : $obj->file('filename'); # open to read
           $obj->file('>filename'); # open to write
           $obj->file('>filename.gz'); # open to write compressed
 Function: Get/set filename; when setting also opens the file and sets fh().
           There is also read support for remote files like
           'ftp://ftp..../file.txt' and it will be downloaded to a temporary
           location and opened.
 Returns : absolute path of file
 Args    : filename

=cut

sub file {
    my ($self, $filename) = @_;
    
    if ($filename) {
        my $in_out = '<';
        if ($filename =~ /^(>+)/) {
            $in_out = $1;
            $filename =~ s/^>+//;
        }
        elsif ($filename =~ /^ftp:|^http:/) {
            $filename = $self->get_remote_file($filename) || $self->throw("Could not download remote file '$filename'");
        }
        
        # avoid potential problems with caller changing dir and things being
        # relative; also more informative and explicit to throw with full path
        $filename = abs_path($filename);
        
        # set up the open command, handling compressed files automatically
        my $open = $filename;
        if ($filename =~ /\.gz$/) {
            if ($in_out eq '<') {
                my $z = IO::Uncompress::Gunzip->new($filename);
                $self->{_filename} = $filename;
                $self->fh($z);
                return $filename;
            }
            else {
                $open = "| gzip -c > $filename";
            }
            $in_out = '';
        }
        
        # go ahead and open it (3 arg form not working when middle is optional)
        open(my $fh, $in_out.$open) || $self->throw("Couldn't open '$open': $!");
        
        $self->{_filename} = $filename;
        $self->fh($fh);
    }
    
    return $self->{_filename};
}

=head2 fh

 Title   : fh
 Usage   : $obj->fh($filehandle);
 Function: Get/set the file handle.
 Returns : filehandle
 Args    : filehandle

=cut

sub fh {
    my ($self, $fh) = @_;
    
    if ($fh) {
        $self->close();
        my $ref = ref($fh) || 'string';
        ($ref eq 'GLOB' || $ref eq 'IO::Uncompress::Gunzip') or $self->throw("fh() takes a filehandle GLOB or IO::Uncompress::Gunzip, not a '$ref'");
        $self->{_fh} = $fh;
        $self->{_fh_count}++; # acts as a unique id so you can know if your
                              # fh has been changed from out under you - _fh
                              # globs are not guaranteed to be unique per script
                              # run!
    }
    
    return $self->{_fh};
}

sub _fh_id {
    my $self = shift;
    return $self->{_fh_count};
}

=head2 seek

 Title   : seek
 Usage   : $obj->seek($pos, $whence);
 Function: Behaves exactly like Perl's standard seek(), except that if the
           filehandle was made by opening a .gz file for reading, you can
           effectively seek backwards.
 Returns : boolean (for success)
 Args    : position to seek to, position to seek from

=cut

sub seek {
    my ($self, $tell, $whence) = @_;
    my $fh = $self->{_fh} || return;
    
    if (ref($fh) eq 'IO::Uncompress::Gunzip') {
        # we can't go backwards, so close and re-open without changing our
        # fh id
        close($fh);
        my $z = IO::Uncompress::Gunzip->new($self->{_filename});
        $self->{_fh} = $z;
        $fh = $z;
    }
    
    CORE::seek($fh, $tell, $whence);
}

=head2 close

 Title   : close
 Usage   : $obj->close();
 Function: Closes the file handle associated with this IO instance.
 Returns : n/a
 Args    : n/a

=cut

sub close {
    my $self = shift;
    
    my $fh = $self->fh();
    undef $self->{_fh};
    
    if ($fh) {
        # don't close STD*
        return if ( \*STDOUT == $fh ||
                    \*STDERR == $fh ||
                    \*STDIN == $fh );
        
        CORE::close($fh);
    }
}

=head2  num_lines

 Title   : num_lines
 Usage   : $obj->num_lines(); 
 Function: Get the number of lines in the current file (that specified by
           file() or fh()).
           NB: bad things will probably happen if you try this on a piped
           filehandle or other unseekable case.
 Returns : int
 Args    : n/a

=cut

sub num_lines {
    my $self = shift;
    my $fh = $self->fh() || return;
    my $tell = tell($fh);
    
    my $lines = 0;
    $self->seek(0, 0);
    $fh = $self->fh(); # if seeking on a .gz, fh will change
    while (<$fh>) {
        $lines++;
    }
    $self->seek($tell, 0);
    
    return $lines;
}

=head2 parse_fod

 Title   : parse_fod
 Usage   : my @directories = $obj->parse_fod('file_of_directories'); 
 Function: Parse a file containing a list of directories.
 Returns : a list consisting of the absolute paths to the directories listed in
           the file
 Args    : filename. To get absolute paths (symlinks are followed), no other
           args. Supply a dir ('' for current dir) to get paths relative
           to it.'/' can be used to get an absolute path that does not follow
           symlinks.

=cut

sub parse_fod {
    my ($self, $fod, $relative) = @_;
    
    -s $fod || $self->throw("fod file '$fod' empty!");
    
    open(my $fodfh, $fod) || $self->throw("Couldn't open fod file '$fod'");
    my %dirs;
    while (<$fodfh>) {
        chomp;
        /\S/ || next;
        unless (-d $_) {
            $self->warn("fod file contained a line that wasn't a directory, ignoring: $_");
            next;
        }
        my $dir = defined $relative ? File::Spec->abs2rel($_, $relative) : abs_path($_);
        if ($relative && $relative eq '/') {
            $dir = '/'.$dir;
        }
        $dirs{$dir} = 1;
    }
    
    my @sorted = sort keys %dirs;
    return @sorted;
}

=head2 parse_fofn

 Title   : parse_fofn
 Usage   : my @files = $obj->parse_fofn('file_of_files'); 
 Function: Parse a file containing a list of files. Lines beginning with # are
           ignored.
 Returns : a list of paths to files in the fofn
 Args    : filename. To get absolute paths (symlinks are followed), no other
           args. Supply a dir ('' for current dir) to get paths relative
           to it. '/' can be used to get an absolute path that does not follow
           symlinks.

=cut

sub parse_fofn {
    my ($self, $fofn, $relative) = @_;
    
    -s $fofn || $self->throw("fofn file '$fofn' empty!");
    
    open(my $fofnfh, $fofn) || $self->throw("Couldn't open fofn file '$fofn'");
    my %files;
    while (<$fofnfh>) {
        chomp;
        /\S/ || next;
        /^#/ && next;
        -f $_ || -l $_ || $self->throw("fofn file contained a line that wasn't a file: $_");
        my $file = defined $relative ? File::Spec->abs2rel($_, $relative) : abs_path($_);
        if ($relative && $relative eq '/') {
            $file = '/'.$file;
        }
        $files{$file} = 1;
    }
    
    my @sorted = sort keys %files;
    return @sorted;
}

=head2  get_remote_file

 Title   : get_remote_file
 Usage   : $obj->get_remote_file('url', save => '/path/to/download/to'); 
 Function: Download a remote file from the internet. Tries to be robust: will
           attempt multiple times to get a file.
 Returns : path to downloaded file on success (otherwise the file won't exist)
 Args    : source url. Optionally:
           save => path to save to. With no path, saves to a temp file.
           md5 => the expected md5 of the file - if it doesn't match what was
                  downloaded, the download will be automatically reattempted up
                  to 3 times. Currently only applies to ftp downloads.

=cut

sub get_remote_file {
    my ($self, $url, %opts) = @_;
    return unless $url;
    $self->{fsu} = VertRes::Utils::FileSystem->new();
    my $local_dir = $self->{fsu}->tempdir();
    my $md5 = $opts{md5};
    
    my $ff = File::Fetch->new(uri => $url);
    my $scheme = $ff->scheme;
    my $host = $ff->host;
    my $path = $ff->path;
    my $basename = $ff->file;
    my $full_path = $path.$basename;
    
    my $out_file = $self->{fsu}->catfile($local_dir, $ff->output_file);
    
    if ($scheme eq 'ftp') {
        # use Net::FTP::Robust, since it's potentially better
        my $ftp;
        if (exists $self->{ftp_objs}->{$host}) {
            $ftp = $self->{ftp_objs}->{$host};
        }
        else {
            $ftp = Net::FTP::Robust->new(Host => $host);
            $self->{ftp_objs}->{$host} = $ftp;
        }
        
        $ftp->get($full_path, $local_dir);
        
        unless (-s $out_file) {
            $self->warn("After 10 automated attempts, failed to download $url at all");
            unlink($out_file);
            return;
        }
        
        if ($md5) {
            my $ok = $self->{fsu}->verify_md5($out_file, $md5);
            
            unless ($ok) {
                my $tries = 0;
                while (! $ok) {
                    unlink($out_file);
                    $ftp->get($full_path, $local_dir);
                    $ok = $self->{fsu}->verify_md5($out_file, $md5);
                    
                    $tries++;
                    last if $tries == 3;
                }
            }
            
            unless ($ok) {
                $self->warn("Tried downloading $url 3 times, but the md5 never matched '$md5'");
                unlink($out_file);
                return;
            }
        }
    }
    else {
        $out_file = $ff->fetch(to => $local_dir) or $self->throw($ff->error);
    }
    
    my $save_file = $opts{save};
    if ($save_file) {
        my $tmp_save_file = $save_file;
        if ($md5) {
            $tmp_save_file .= '.tmp';
        }
        
        my $success = File::Copy::move($out_file, $tmp_save_file);
        $success || $self->throw("Failed to move $out_file to $tmp_save_file");
        
        if ($md5) {
            # recheck md5 after the move, since we may have crossed filesystems
            my $ok = $self->{fsu}->verify_md5($tmp_save_file, $md5);
            if ($ok) {
                my $success = File::Copy::move($tmp_save_file, $save_file);
                $success || $self->throw("Failed to move $tmp_save_file to $save_file");
            }
            else {
                $self->warn("md5 match failed after moving good download file from $out_file to $tmp_save_file; will unlink it\n");
                unlink($tmp_save_file);
                return;
            }
        }
        
        return $save_file;
    }
    return $out_file;
}

1;
