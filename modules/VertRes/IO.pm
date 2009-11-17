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

# generally useful file-related functions
my $io = VertRes::IO->new();
my $base_dir = $io->catfile($G1K, 'DATA');
my @paths = $io->get_filepaths($base_dir, suffix => 'fastq.gz');
my ($tempfh, $tempfile) = $io->tempfile();
my $tempdir = $io->tempdir;
$io->rmtree($directory_structure_safe_to_delete); # !! CAREFULL !!

=head1 DESCRIPTION

Handles IO including reading and writing to compressed files. Also provides
aliases to commonly needed file-related functions: tempfile, tempdir, catfile,
rmtree.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::IO;

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Temp;
use File::Spec;
use File::Basename;
require File::Path;
require File::Copy;
use IO::Uncompress::Gunzip;
use File::Fetch;
use Net::FTP::Robust;
use Digest::MD5;

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

=head2 get_filepaths

 Title   : get_filepaths
 Usage   : my @paths = $obj->get_filepaths('base_dir'); 
 Function: Get the absolute paths to all files in a given directory and all its
           subdirectories.
 Returns : a list of filepaths
 Args    : path to base directory
           optionally, the following named args select out only certain paths
           according to if they match the supplied regex string(s)
           filename => regex (whole basename of file must match regex)
           prefix   => regex (basename up to the final '.' must match regex)
           suffix   => regex (everything after the final '.' must match regex;
                              the '.' in .gz is not treated as the final '.' for
                              this purpose)
           dir      => regex (return directory paths that match, instead of
                              files - disables above 3 options)
           subdir   => regex (at least one of a file/dir's parent directory must
                              match regex)

=cut

sub get_filepaths {
    my ($self, $dir, %args) = @_;
    
    $dir = abs_path($dir);
    my $wanted_dir = $args{dir};
    opendir(my $dir_handle, $dir) || $self->throw("Couldn't open dir '$dir': $!");
    
    my @filepaths;
    foreach my $thing (readdir($dir_handle)) {
        next if $thing =~ /^\.+$/;
        my $orig_thing = $thing;
        $thing = $self->catfile($dir, $thing);
        
        # recurse into subdirs
        if (-d $thing) {
            if ($wanted_dir && $orig_thing =~ /$wanted_dir/) {
                if (($args{subdir} && $thing =~ /$args{subdir}/) || ! $args{subdir}) {
                    push(@filepaths, $thing);
                }
            }
            
            push(@filepaths, $self->get_filepaths($thing, %args));
            next;
        }
        
        next if $wanted_dir;
        
        # check it matches user's regexs
        my $ok = 1;
        my ($basename, $directories) = fileparse($thing);
        my $gz = '';
        if ($basename =~ s/\.gz$//) {
            $gz = '.gz';
        }
        my ($prefix, $suffix) = $basename =~ /(.+)\.(.+)$/;
        unless ($prefix) {
            # we have a .filename file
            $suffix = $basename;
        }
        $suffix .= $gz;
        $basename .= $gz;
        while (my ($type, $regex) = each %args) {
            if ($type eq 'filename') {
                $basename =~ /$regex/ || ($ok = 0);
            }
            elsif ($type eq 'prefix') {
                unless ($prefix) {
                    $ok = 0;
                }
                else {
                    $prefix =~ /$regex/ || ($ok = 0);
                }
            }
            elsif ($type eq 'suffix') {
                unless ($suffix) {
                    $ok = 0;
                }
                else {
                    $suffix =~ /$regex/ || ($ok = 0);
                }
            }
            elsif ($type eq 'subdir') {
                $directories =~ /$regex/ || ($ok = 0);
            }
        }
        
        push(@filepaths, $thing) if $ok;
    }
    
    return @filepaths;
}

=head2 parse_fod

 Title   : parse_fod
 Usage   : my @directories = $obj->parse_fod('file_of_directories'); 
 Function: Parse a file containing a list of directories.
 Returns : a list consisting of the absolute paths to the directories listed in
           the file
 Args    : filename

=cut

sub parse_fod {
    my ($self, $fod) = @_;
    
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
        my $dir = abs_path($_);
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
           to it.

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
        $files{$file} = 1;
    }
    
    my @sorted = sort keys %files;
    return @sorted;
}

=head2 tempfile

 Title   : tempfile
 Usage   : my ($handle, $tempfile) = $obj->tempfile(); 
 Function: Get a temporary filename and a handle opened for writing and
           and reading. Just an alias to File::Temp::tempfile.
 Returns : a list consisting of temporary handle and temporary filename
 Args    : as per File::Temp::tempfile

=cut

sub tempfile {
    my $self = shift;
    
    my $ft = File::Temp->new(@_);
    push(@{$self->{_fts}}, $ft);
    
    return ($ft, $ft->filename);
}

=head2  tempdir

 Title   : tempdir
 Usage   : my $tempdir = $obj->tempdir(); 
 Function: Creates and returns the name of a new temporary directory. Just an
           alias to File::Temp::tempdir.
 Returns : The name of a new temporary directory.
 Args    : as per File::Temp::tempdir

=cut

sub tempdir {
    my $self = shift;
    
    my $ft = File::Temp->newdir(@_);
    push(@{$self->{_fts}}, $ft);
    
    return $ft->dirname;
}

=head2  catfile

 Title   : catfile
 Usage   : my ($path) = $obj->catfile('dir', 'subdir', 'filename'); 
 Function: Constructs a full pathname in a cross-platform safe way. Just an
           alias to File::Spec->catfile.
 Returns : the full path
 Args    : as per File::Spec->catfile

=cut

sub catfile {
    my $self = shift;
    return File::Spec->catfile(@_);
}

=head2  rmtree

 Title   : rmtree
 Usage   : $obj->rmtree('dir'); 
 Function: Remove a full directory tree - files and subdirs. Just an alias to
           File::Path::rmtree.
 Returns : n/a
 Args    : as per File::Path::rmtree

=cut

sub rmtree {
    my $self = shift;
    return File::Path::rmtree(@_);
}

=head2  copy

 Title   : copy
 Usage   : $obj->copy('source.file', 'dest.file'); 
 Function: Copy a file and check that the copy is identical to the source
           afterwards.
 Returns : boolean (true on success; on failure the destination path won't
           exist)
 Args    : source file path, output file path. Optionally, the number of times
           to retry the copy if the copy isn't identical, before giving up
           (default 3).

=cut

sub copy {
    my ($self, $source, $dest, $max_retries) = @_;
    unless (defined $max_retries) {
        $max_retries = 3;
    }
    
    for (1..$max_retries) {
        my $success = File::Copy::copy($source, $dest);
        if ($success) {
            my $diff = `diff $source $dest`;
            return 1 unless $diff;
        }
    }
    
    return 0;
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
    my $local_dir = $self->tempdir();
    my $md5 = $opts{md5};
    
    my $ff = File::Fetch->new(uri => $url);
    my $scheme = $ff->scheme;
    my $host = $ff->host;
    my $path = $ff->path;
    my $basename = $ff->file;
    my $full_path = $path.$basename;
    
    my $out_file = $self->catfile($local_dir, $ff->output_file);
    
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
            my $ok = $self->verify_md5($out_file, $md5);
            
            unless ($ok) {
                my $tries = 0;
                while (! $ok) {
                    unlink($out_file);
                    $ftp->get($full_path, $local_dir);
                    $ok = $self->verify_md5($out_file, $md5);
                    
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
            my $ok = $self->verify_md5($tmp_save_file, $md5);
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

=head2  verify_md5

 Title   : verify_md5
 Usage   : if ($obj->verify_md5($file, $md5)) { #... }
 Function: Verify that a given file has the given md5.
 Returns : boolean
 Args    : path to file, the expected md5 (hexdigest as produced by the md5sum
           program) as a string

=cut

sub verify_md5 {
    my ($self, $file, $md5) = @_;
    my $new_md5 = $self->calculate_md5($file);
    return $new_md5 eq $md5;
}

=head2  calculate_md5

 Title   : calculate_md5
 Usage   : my $md5 = $obj->calculate_md5($file)
 Function: Calculate the md5 of a file.
 Returns : hexdigest string
 Args    : path to file

=cut

sub calculate_md5 {
    my ($self, $file) = @_;
    
    open(my $fh, $file) || $self->throw("Could not open file $file");
    binmode($fh);
    my $dmd5 = Digest::MD5->new();
    $dmd5->addfile($fh);
    my $md5 = $dmd5->hexdigest;
    
    return $md5;
}

1;
