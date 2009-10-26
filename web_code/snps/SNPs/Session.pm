#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
#
#--------------- Session --------------------------------------------
#
# For storing and retrieving stored data
#

package SNPs::Session;

use strict;
use warnings;

use Compress::Zlib;
our $unique_id;

sub new
{
    my ($class,$args) = @_;

    my $self = $args ? $args : {};
    if ( !exists($$self{'dbstore'}) )
    {   
        if ( !exists($$self{'sw'}) ) { die "Missing the 'dbstore' or 'sw' option.\n"; }
        $$self{'dbstore'} = $$self{'sw'}->dbstore();
    }
    if ( !$unique_id ) { $unique_id = $$ . time; } 
    if ( !exists($$self{'prefix'}) ) { $$self{'prefix'}='modelorgs/mousegenomes/'; }
    if ( $$self{'hours'} ) { $$self{'hours'}=1; }
    bless $self, ref($class) || $class;
    return $self;
}

sub append
{
    my ($self,$data) = @_;
    $$self{'data'} .= $data;
    return;
}

sub store
{
    my ($self) = @_;
    if ( !$$self{'id'} ) { $self->id(); }
    if ( !$$self{'data'} ) { return 0; }
    my $data = compress($$self{'data'});
    my $status = $$self{'dbstore'}->set($data,$$self{'prefix'}.$$self{'id'},$$self{hours});
    if ( !$status ) { die "Unable to cache the data.\n"; }
    return $$self{'id'};
}

sub retrieve
{
    my ($self,$id) = @_;
    my $key = $$self{'prefix'}.$id;
    $$self{'dbstore'}->refresh($key,$$self{hours});
    my $data = $$self{'dbstore'}->get($key);
    if ( $data ) { $data = uncompress($data); }
    return $data;
}

sub id
{
    my ($self) = @_;
    if ( !$$self{'id'} ) { $$self{'id'} = $unique_id++; }
    return $$self{'id'};
}

1;

