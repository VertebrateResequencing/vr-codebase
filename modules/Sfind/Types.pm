package Sfind::Types;
=head1 NAME

Sfind::Types - data types for Sfind.

=head1 SYNOPSIS

package Foo;
use Moose;
use Sfind::Types qw(MysqlDateTime);

# use the exported constants as type names
has 'bar' => (
    isa => MysqlDateTime,
    is => 'rw',
    coerce => 1,    # allow type to take a mysql format string
);

=head1 DESCRIPTION

Defines data types for Sfind.

=head1 CONTACT

jws@sanger.ac.uk

=cut

use MooseX::Types -declare => [qw( MysqlDateTime YesNoBool)];
use MooseX::Types::Moose qw(Str Bool);
use DateTime;
use DateTime::Format::MySQL;
use strict;
use warnings;

class_type MysqlDateTime, { class => 'DateTime' };

coerce MysqlDateTime,
    from Str,
    via { DateTime::Format::MySQL->parse_datetime($_) };

1;
