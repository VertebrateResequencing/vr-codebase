=head1 NAME

VertRes::Wrapper::iRODS 

=head1 SYNOPSIS

use VertRes::Wrapper::iRODS;

=head1 DESCRIPTION

Wrapper for Wrapper::iRODS (Integrated Rule Oriented Data Systems) icommands.
New sequencing data is written into the NPG iRODS system, and this module
provides a perl-friendly way to query, view, and retrieve those data.

NB: you need to have logged into irods first, with iinit or kinit, before using
the module.

=head1 AUTHOR

Jim Stalker jws@sanger.ac.uk

=cut

package VertRes::Wrapper::iRODS;

use strict;
use warnings;
use base qw(VertRes::Base);


our $defaults = { 'icommands'         => '/software/irods/icommands/bin',
                };


=head2 new

    Description: create new irods object
    Arg [1]    : optional iRODS settings, currently just location of irods icommand binaries
    Example    : my $irods = VertRes::iRODS->new('icommands'=>'/usr/bin/irods');
    Returntype : VertRes::iRODS object

=cut

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(%$defaults, @args);
    bless($self,$class);

    # set up iRODS environment
    # NB - this currently doesn't work.  Need to run iinit first
    #$ENV{'irodsEnvFile'} = $self->{'irodsEnvFile'};
    #$ENV{'irodsAuthFileName'} = $self->{'irodsAuthFileName'};

    return $self;
}


=head2 find_file_by_name

    Description : lists irods location for a specific filename
    Arg [1]     : file name, e.g. '1234_5.bam'
    Example     : my $ifile = $irods->find_file_by_name('1234_5.bam');
    Returntype  : string of irods file location, e.g. '/seq/1234/1234_5.bam'

=cut

use Memoize;
memoize('find_file_by_name');

sub find_file_by_name {
    my ($self, $name) = @_;
    my $cmd = join "/",($self->{icommands},qq(iquest -z seq "SELECT COLL_NAME, DATA_NAME WHERE DATA_NAME = '$name'"));
    open(my $irods, "$cmd |");
    my ($path, $filename);
    while (<$irods>) {
        chomp;
        # output looks like:
        # Zone is seq
        # COLL_NAME = /seq/5150
        # DATA_NAME = 5150_1#1.bam
        # -------------------------
        # 
        # or an error is thrown, and only:
        # Zone is seq
        # is output

        if (/^COLL_NAME = (.+)$/){
            $path = $1;
        }
        if (/^DATA_NAME = (.+)$/){
            $filename = $1;
        }
    }
    close $irods;
    my $file = undef;
    if ($path && $filename){
        unless ($filename eq $name){
            die "Error: $filename should be the same as $name\n";
        }
        $file = join "/",($path, $filename);
    }
    return $file;
}



=head2 find_files_by_run_lane

    Description : lists irods locations for a specific run & lane.  Can return multiple files, due to multiplexing
    Arg [1]     : run id, e.g. 1234
    Arg [2]     : lane id, e.g. 5
    Example     : my @files = $irods->find_files_by_run_lane('1234','5');
    Returntype  : arrayref of irods file locations

=cut

sub find_files_by_run_lane {
    my ($self, $run, $lane) = @_;
    unless ($run && $lane){
         $self->throw("Missing parameters: run and lane.\n");
    }
    my $cmd = join "/",($self->{icommands},"imeta -z seq qu -d id_run = $run and lane = $lane");

    open(my $irods, "$cmd |");

    my @out;
    my $path = '';
    while (<$irods>) {
        # output looks like:
        #collection: /seq/5253
        #dataObj: 5253_1.bam
        #----
        #collection: /seq/5253
        #dataObj: 5253_2.bam
        # or
        # No rows found

        if (/^collection: (.+)$/){
            $path = $1;
        }
        if (/^dataObj: (.+)$/){
            push @out, "$path/$1";
        }
    }
    close $irods;
    return \@out;
}


=head2 find_files_by_name

    Description : lists irods locations for a specific file.  Can return multiple files, due to multiplexing
    Arg [1]     : file name, e.g. 1234_5.bam
    Example     : my @files = $irods->find_files_by_run_lane('1234_5.bam');
    Returntype  : arrayref of irods file locations

=cut

sub find_files_by_name {
    my ($self, $file) = @_;
    unless ($file){
         $self->throw("Missing parameters: run and lane.\n");
    }

    if ( !($file=~/(\S+)_(\S+)\.bam$/) ) { $self->throw("TODO: $file"); }
    return $self->find_files_by_run_lane($1,$2);
}


=head2 get_file_md5
    Description : return the md5 for a file in irods
    Arg [1]     : irods file location
    Example     : my $md5 = $irods->get_file_md5('/seq/1234/1234_5.bam');
    Returntype  : md5 string

=cut

sub get_file_md5 {
    my ($self, $file) = @_;
    my $cmd = join "/",($self->{icommands},"ichksum $file");

    my $md5 = `$cmd`;
    chomp $md5;
    $md5 =~s/.*\s//;
    return $md5;
}


=head2 get_file_size
    Description : return the byte size of a file in irods
    Arg [1]     : irods file location
    Example     : my $size = $irods->get_file_size('/seq/1234/1234_5.bam');
    Returntype  : integer number of bytes

=cut

sub get_file_size {
    my ($self, $file) = @_;
    my $cmd = join "/",($self->{icommands},"ils -l $file");

    open(my $irods, "$cmd |");
    my $line = <$irods>;
    #   srpipe            0 res-g2                17084186273 2010-10-20.19:13 & 5330_1.bam
    my @fields = split ' ', $line;
    return $fields[3];
}


=head2 get_file

    Description : retrieves a file from irods.  Does checksum check after copy.  Will overwrite if file already exists at destination.
    Arg [1]     : irods file location
    Arg [2]     : optional local filesystem destination
    Example     : $irods->get_file('/seq/1234/1234_5.bam', "/datadir/allfiles/1234_5.bam');
    Returntype  : systemcall return value

=cut

sub get_file {
    my ($self, $file, $dest) = @_;
    my $cmd = join "/",($self->{icommands},"iget");
    # -K: checksum
    # -Q: use UDP rather than TCP
    # -f: force overwrite
    #my @args = ($cmd, "-K", "-Q", "-f", $file, $dest);
    my @args = ($cmd, "-K", "-f", $file, $dest);
    return system(@args);
}


=head2 get_total_reads

    Description : returns the total nubmer of reads in the bam file
    Arg [1]     : irods file location
    Example     : my $total_reads = $irods->get_total_reads('/seq/1234/1234_5.bam');
    Returntype  : integer number of reads

=cut

sub get_total_reads {
  my ($self, $file) = @_;
  my $cmd = join "/",($self->{icommands},"imeta ls -d $file total_reads | grep value | awk '{ print $2 }'");

  my $total_reads = `$cmd`;
  chomp $total_reads;
  return $total_reads;
}

