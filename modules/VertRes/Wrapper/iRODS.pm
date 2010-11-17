=head1 NAME

VertRes::Wrapper::iRODS 

=head1 SYNOPSIS

use VertRes::Wrapper::iRODS;

=head1 DESCRIPTION

Wrapper for Wrapper::iRODS (Integrated Rule Oriented Data Systems) icommands.
New sequencing data is written into the NPG iRODS system, and this module
provides a perl-friendly way to query, view, and retrieve those data.

=head1 AUTHOR

Jim Stalker jws@sanger.ac.uk

=cut

package VertRes::Wrapper::iRODS;

use strict;
use warnings;
use base qw(VertRes::Base);

# you need to install your .irodsEnv + .irodsA file on the server side in a
# given path <path1>.  then in your code, you can set up the env variables
# "irodsEnvFile" and "irodsAuthFileName" with the function putenv where you
# give the full path name of the 2 connexion files.  once you have that, please
# note that the "iinit" step is not required anymore. 

our $defaults = {   'irodsEnvFile'      => '/lustre/scratch102/conf/irodsEnv',
                    'irodsAuthFileName' => '/lustre/scratch102/conf/irodsA',
                    'icommands'         => '/software/irods/icommands/bin',
                };


=head2 new

    Arg [1]    : optional iRODS connection settings pointing to an irodsEnv file and irodsAuth file.
    Example    : my $irods = VertRes::iRODS->new('irodsEnvFile'=>'~/.irods/.irodsEnv','irodsAuthFileName'=>'~/.irods/.irodsA');
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


sub ils {
    my ($self) = @_;
    my $cmd = join "/",($self->{icommands},'ils');
    my $out = `$cmd`;
    return $out;
}

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


sub get_file_md5 {
    my ($self, $file) = @_;
    my $cmd = join "/",($self->{icommands},"ichksum $file");

    my $md5 = `$cmd`;
    chomp $md5;
    $md5 =~s/.*\s//;
    return $md5;
}


sub get_file_size {
    my ($self, $file) = @_;
    my $cmd = join "/",($self->{icommands},"ils -l $file");

    open(my $irods, "$cmd |");
    my $line = <$irods>;
    #   srpipe            0 res-g2                17084186273 2010-10-20.19:13 & 5330_1.bam
    my @fields = split ' ', $line;
    return $fields[3];
}


sub get_file {
    my ($self, $file, $dest) = @_;
    my $cmd = join "/",($self->{icommands},"iget");
    my @args = ($cmd, "-K", "-Q", "-f", $file, $dest);
    return system(@args);
}

