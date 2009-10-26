#
# Author:    Petr Danecek (pd3@sanger.ac.uk)    Team 145
#
#--------------- DataWriter ---------------------------------------------
#
# The result files download, no header or footer.
#

package SNPs::DataWriter;

use strict;
use warnings;
use base qw(SNPs::Writer);

sub new
{
    my ($class, $args) = @_;
    my $self = $class->SUPER::new($args);
    if ( !exists($$self{'fname'}) ) { $$self{'fname'}=''; }
    return $self;
}

sub fname
{
    my ($self,$fname) = @_;
    if ( $fname ) { $$self{'fname'} = $fname; }
    return $$self{'fname'};
}

sub header
{
    my ($self,@args) = @_;

    return qq[Content-Type: application/download
Content-Disposition: attachment; filename=$$self{'fname'}
Content-Description: File Transfer
Content-Transfer-Encoding: binary

];
}

sub print_footer
{
    return;
}

1;

