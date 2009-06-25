#!/usr/bin/env perl
use lib qw(/nfs/team71/phd/dmc/myperl/lib/site_perl);
use lib qw(/nfs/team71/phd/dmc/myperl/share/perl);
use XML::Validator::Schema;
use XML::SAX::ParserFactory;

use strict;

# create a new validator object, using foo.xsd
my $xsd = shift @ARGV;
my $validator = XML::Validator::Schema->new(file => $xsd);

#
# create a SAX parser and assign the validator as a Handler
#
my $parser = XML::SAX::ParserFactory->parser(Handler => $validator);

#
# validate foo.xml against foo.xsd
#

foreach my $f (@ARGV) {
    eval { $parser->parse_uri($f) };
    die "File failed validation: $@" if $@;
}
