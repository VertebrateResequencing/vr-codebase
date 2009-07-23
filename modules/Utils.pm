package Utils;

use strict;
use warnings;
use Carp qw(confess);
use IPC::Open3 'open3';
local $SIG{CHLD} = 'IGNORE';

=pod

=head1 NAME

FastQ

=head1 SYNOPSIS

A module for fastq files manipulation.

=cut

=head1 AUTHORS

Petr Danecek, I<pd3@sanger.ac.uk>

=cut

=head1 DATA

Use C<$Utils::OK> and C<$Utils::Error> on exit, so that the caller can check the return status.

=cut

$Utils::OK    = 0;
$Utils::Error = 1;

$Utils::remove_on_exit = {};

=head1 METHODS

=head2 error

    Description     : Prints backtrace, error message (join('',@msg) is used on arrays) and exists.

=cut

sub error
{
    my (@msg) = @_;

    chomp(my $cwd=`pwd`);

    my $bt = backtrace();
    my $msg = join('',@$bt);
    $msg .= "Current directory: $cwd\n\n";

    if ( scalar @msg )
    {
        $msg .= join('',@msg) . "\n\n";
    }

    $? = $Utils::Error;
    die $msg;
}


sub cleanup
{
    my @removed = ();
    for my $file (keys %$Utils::remove_on_exit)
    {
        if ( -f $file ) 
        { 
            push @removed, $file;
            unlink $file; 
        }
        # todo: directories?
    }

    # Let the user know..
    if ( scalar @removed )
    {
        print STDERR "Cleaning files:\n";
        for my $file (@removed)
        {
            print STDERR "\t $file\n";
            delete $$Utils::remove_on_exit{$file};
        }
    }

    if ( scalar keys %$Utils::remove_on_exit)
    {
        print STDERR "Failed to remove:\n";
        for my $file (keys %$Utils::remove_on_exit)
        {
            print STDERR "\t $file\n";
        }
    }

    return;
}


sub exit_cleanup
{
    Utils::error(@_);
}


sub backtrace
{
    my @msg = ();

    push @msg, "Backtrace:\n";

    my ($package, $filename, $line, $subroutine);
    my $i=0;
    while ( (($package, $filename, $line, $subroutine)=caller($i)) )
    {
        push @msg, "\t $filename:$line  ->  $subroutine\n";
        $i++;
    }

    return \@msg;
}


=head2 CMD

    Arg[1]   : The command to execute
    Arg[2]   : Hash with options controlling its behaviour
    Returns  : The output returned by the command.
    Example  : my @out = Utils::CMD("ls",{'verbose'=>1,'exit_on_error'=>0});
    Options  :
                exit_on_error .. should the encountered errors be ignored?
                logfile       .. where should be the command and the output logged?
                verbose       .. print what's being done.

=cut

sub CMD
{
    my ($cmd, $options) = @_;

    $options = {} unless $options;
    $$options{'exit_on_error'} = 1 unless exists($$options{'exit_on_error'});
    $$options{'verbose'}       = 0 unless exists($$options{'verbose'});

    print "$cmd\n" unless !$$options{'verbose'};
    log_msg($$options{'logfile'},"$cmd\n") unless !exists($$options{'logfile'});
    my @out = `$cmd`;

    if ( $? )
    {
        my @msg = ();
        push @msg, "The command \"$cmd\" returned non-zero status $?";
        if ( $! ) 
        { 
            push @msg, ": $!\n"; 
        }
        else 
        { 
            push @msg, ".\n"; 
        }
        if ( scalar @out )
        {
            push @msg, @out;
        }
        my $bt = Utils::backtrace();
        my $msg = join('',@$bt) . "\n" . join('',@msg) . "\n";
        log_msg($$options{'logfile'},$msg) unless !exists($$options{'logfile'});

        Utils::error(@msg) unless !exists($$options{'exit_on_error'});
        print STDERR $msg;
    }

    return (@out);
}


=head2 create_dir

    Arg[1]          : The directory path.
    Description     : Creates a directory recursively, much the same like "mkdir -p"
    Returntype      : none, exits on error

=cut

sub create_dir
{
    my ($path) = @_;

    my @items = split m{/}, $path;
    my $tmp_path  = ($path =~ m{^/}) ? '/' : '';
    for my $item (@items)
    {
        $tmp_path .= $item . '/';
        if ( -d $tmp_path ) { next }
        mkdir($tmp_path) or error("create_dir \"$path\" ($tmp_path): $!\n");
    }
    return;
}


=head2 uncompressed_gz_size

    Arg[1]          : The gzipped file name
    Description     : Determines what is the size of the uncompressed file.
    Returntype      : The size of the uncompressed file in bytes.

=cut

sub uncompressed_gz_size
{
    my ($file) = @_;

    # From some reason, this does not work ("Use of uninitialized value in <HANDLE>")
    #
    #   my ($writer, $reader, $err);
    #   open3($writer, $reader, $err, "gzip -l $file");
    #   my @output = <$reader>;
    #   my @errors = <$err>;

    my ($writer, $reader);
    open3($writer, $reader, \*ERR, "gzip -l $file");
    my @output = <$reader>;
    my @errors = <ERR>;

    if ( scalar @errors ) { error("gzip -l $file:\n", @errors) }

    # Expected output:
    #   compressed        uncompressed  ratio uncompressed_name
    #   1772198718          4147914907  57.3% mouse-2470_1_2.fastq
    if ( scalar @output != 2 || !($output[1]=~/^\s+\d+\s+(\d+)/ ) ) { error("Uh, expected something else:\n", @output) }
    return $1;
}

=head2 shadow_hierarchy

    Arg[1]      : The root of the existing hierarchy
    Arg[2]      : The root of the new hierarchy
    Arg[3]      : The list of files to be symlinked from the new hierarchy to the old one.
    Arg[4]      : The list of the symlink names. If not present, the name of the original file will be used. [optional]
    Description : Creates a symbolic link using relative paths. The parameters must be absolute.
    Returntype  : None
    Example     : @files=('1/a','1/b','1/2/c'); shadow_hierarchy('1','2',\@files);

=cut

sub shadow_hierarchy
{
    my ($old_root,$new_root,$file_list,$new_names) = @_;

    my $nfiles = scalar @$file_list;
    for (my $ifile=0; $ifile<$nfiles; $ifile++)
    {
        my $file = $$file_list[$ifile];

        my $fname_relative = $file;
        if ( $fname_relative=~m{^/} ) { $fname_relative=~s{^$old_root/?}{}; }
        if ( $fname_relative=~m{^/} ) { error("FIXME: $file, $old_root, $fname_relative\n") }
        
        if ( !($fname_relative=~m{/?([^/]+)$}) ) { error("FIXME: could not parse $fname_relative.\n") }
        my $relative_dir = $`;
        my $base_name    = $1;
        my $new_name = defined($new_names) ? $$new_names[$ifile] : $base_name;

        create_dir("$new_root/$relative_dir");

        my $old = "$old_root/$fname_relative";
        my $new = "$new_root/$relative_dir/$new_name";

        if ( -e $new ) { next }

        symlink($old,$new) or error("symlink $old $new: $!");
    }
    return;
}


=head2 relative_symlink

    Arg[1]          : The existing file
    Arg[2]          : The symlink to be created
    Description     : Creates a symbolic link using relative paths. The parameters must be absolute.
    Examples        : 
                        relative_symlink('/a/b/c/d/xxx','/a/b/yyy') will create 'yyy -> c/d/xxx' in /a/b/
                        relative_symlink('/a/b/xxx','/a/b/c/d/yyy') will create 'yyy -> ../../xxx' in /a/b/c/d/
    Returntype      : None

=cut

sub relative_symlink
{
    my ($cur,$new) = @_;

    chomp(my ($cwd) = CMD('pwd'));

    # In fact, this is more complex to program (imagine arguments like "././a/b/../../c/d/"),
    #   but we do not need this to be so bullet-proof right now.

    if ( !($cur =~ m{^/} ) ) { $cur = "$cwd/$cur"; }
    if ( !($new =~ m{^/} ) ) { $new = "$cwd/$new"; }

    $cur =~ s{/$}{};
    $new =~ s{/$}{};

    # Replace multiple slashes //
    $cur =~ s{/+}{/}g;
    $new =~ s{/+}{/}g;

    my @cur_els = split m{/}, $cur;
    my @new_els = split m{/}, $new;

    my $common_prefix = '';
    while ( scalar @cur_els > 1 && $cur_els[0] eq $new_els[0] )
    {
        $common_prefix .= '/' . shift @cur_els;
        shift @new_els;
    }
    if ( !scalar @cur_els || !scalar @new_els ) { Utils::error("FIXME: $cur, $new\n") }

    # How many times should we go up ../../../..
    my @uphill_path = ();
    my @new_dir  = splice @new_els,0,-1;
    foreach my $el (@new_dir)
    {
        push @uphill_path, '..';
    }

    # The directory of the symlink
    my $dir = $common_prefix;
    if ( scalar @new_dir )
    {
        $dir .= '/' . join('/', @new_dir);
    }
    chdir($dir) or Utils::error("chdir $dir: $!");

    # The uphill + downhill path
    my $path = '';
    if ( scalar @uphill_path ) 
    {
        $path = join('/', @uphill_path) . '/';
    }
    $path .= join('/', @cur_els);

    #print STDERR "cd $dir; symlink '$path' '$new_els[0]'\n";
    symlink($path,$new_els[0]);

    chdir($cwd) or Utils::error("chdir $cwd: $!");
    return;
}


=head2 basename

    Arg[1]      : The file name.
    Returntype  : The path '/some/path/basename.suffix' splitted to ('/some/path','basename','.suffix')

=cut

sub basename
{
    my ($path) = @_;
    if ( !($path =~ m{/?([^/]+?)(\.?[^.]*)$}) ) { Utils::error("FIXME: could not parse \"$path\".\n") }
    my $dir  = $` ? $` : '';
    my $base = $1 ? $1 : '';
    my $suff = $2 ? $2 : '';
    return ($dir,$base,$suff);
}

sub file_newer
{
    my ($afile,$bfile) = @_;
    my (@astat) = stat($afile) or Utils::error("stat $afile: $!");
    my (@bstat) = stat($bfile) or Utils::error("stat $bfile: $!");

    if ( $astat[9]>$bstat[9] ) { return 1 }
    return 0;
}

=head2 fai_chromosome_lengths

    Arg[1]      : The .fai file.
    Description : Reads the .fai file and returns a hash with chromosome lengths.

=cut

sub fai_chromosome_lengths
{
    my ($fai_file) = @_;

    my $out = {};

    open(my $fh, '<', $fai_file) or Utils::error("$fai_file: $!");
    while (my $line=<$fh>)
    {
        if ( !($line=~/^(\S+)\s+(\d+)\s+/) ) { Utils::error("Unexpected format of $fai_file:\n$line\n"); }
        $$out{$1} = $2;
    }
    close $fh;

    return $out;
}


=head2 cmp_mixed

    Description : Compares both numbers and strings, numbers go first.
    Example     : @a = sort cmp_mixed ('1','2','X','MT');

=cut

sub cmp_mixed($$)
{
    my ($a,$b) = @_;

    if ( $a=~/^\d+$/ )
    {
        if ( $b=~/^\d+$/ ) { return $a<=>$b }
        return -1;
    }
    elsif ( $b=~/^\d+$/ )
    {
        return 1;
    }
    return $a cmp $b;
}


=head2 log_msg

    Description : 
    Arg[1]      : 

=cut

sub log_msg
{
    my ($log_file, $str) = @_;
    open(my $fh, '>>', $log_file) or error("$log_file: $!");
    print $fh scalar gmtime, "\n";
    print $fh $str, "\n";
    close $fh;
}


1;

