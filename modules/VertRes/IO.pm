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
        
        # avoid potential problems with caller changing dir and things being
        # relative; also more informative and explicit to throw with full path
        $filename = abs_path($filename);
        
        # set up the open command, handling compressed files automatically
        my $open = $filename;
        if ($filename =~ /\.gz$/) {
            if ($in_out eq '<') {
                $open = "zcat $filename |";
            }
            else {
                $open = "| gzip -c $filename";
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
        $ref eq 'GLOB' or $self->throw("fh() takes a filehandle GLOB, not a '$ref'");
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
    seek($fh, 0, 0);
    while (<$fh>) {
        $lines++;
    }
    seek($fh, $tell, 0);
    
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

1;
