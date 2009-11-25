=head1 NAME

VertRes::Utils::FileSystem - do filesystem manipulations

=head1 SYNOPSIS

use VertRes::Utils::FileSystem;

my $fsu = VertRes::Utils::FileSystem->new();

# ...


# generally useful file-related functions
my $base_dir = $fsu->catfile($G1K, 'META');
my @paths = $fsu->get_filepaths($base_dir, suffix => 'fastq.gz');
my ($tempfh, $tempfile) = $fsu->tempfile();
my $tempdir = $fsu->tempdir;
$fsu->rmtree($directory_structure_safe_to_delete); # !! CAREFULL !!

=head1 DESCRIPTION

Provides functions related to storing/getting things on/from the file-system.

Also provides aliases to commonly needed file-related functions: tempfile,
tempdir, catfile, rmtree.

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Utils::FileSystem;

use strict;
use warnings;
use Cwd qw(abs_path);
use File::Temp;
use File::Spec;
use File::Basename;
require File::Path;
require File::Copy;
use Digest::MD5;

use base qw(VertRes::Base);

=head2 new

 Title   : new
 Usage   : my $self = $class->SUPER::new(@args);
 Function: Instantiate a new VertRes::Utils::FileSystem object.
 Returns : $self hash-ref blessed into your class
 Args    : n/a

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);
    
    return $self;
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

=head2 tempdir

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

=head2 catfile

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

=head2 rmtree

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

=head2 copy

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

=head2 verify_md5

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

=head2 calculate_md5

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
