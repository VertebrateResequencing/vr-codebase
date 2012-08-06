=head1 NAME

Pathogens::Variant::Exception  - Basic exception class. It subclasses the CPAN's Error.pm (http://search.cpan.org/~shlomif/Error-0.17016/lib/Error.pm).

=head1 SYNOPSIS

use Pathogens::Variant::Exception;
throw Pathogens::Variant::Exception( {text => "I have failed for some reason"} );;


#or something like this:
use Pathogens::Variant::Exception qw(:try);

try {
    do_some_stuff();
    die "error!" if $condition;
    throw Error::Simple "Oops!" if $other_condition;
}
catch Error::IO with {
    my $E = shift;
    print STDERR "File ", $E->{'-file'}, " had a problem\n";
}
except {
    my $E = shift;
    my $general_handler=sub {send_message $E->{-description}};
    return {
        UserException1 => $general_handler,
        UserException2 => $general_handler
    };
}
otherwise {
    print STDERR "Well I don't know what to say\n";
}
finally {
    close_the_garage_door_already(); # Should be reliable
}; # Don't forget the trailing ; or you might be surprised

=cut

package Pathogens::Variant::Exception;
use base qw(Error);
use overload ('""' => 'stringify', fallback => 1); #this is what we do if called in string context.

sub new {

    my ($proto, $params) = @_;

    my $class = ref($proto) || $proto;
    my $text  = defined $$params{text} ? $$params{text} : ' ';
    my $value = defined $$params{value} ? $$params{value} : ' ';
    local $Error::Depth = $Error::Depth + 1;
    local $Error::Debug = 1;                   # Enables storing of stacktrace
    my $exception = $class->SUPER::new( -text => $text, -value => $value );
    return $exception;

}

sub stringify {

    my $self = shift;
    my $class = ref($self) || $self;
    return $self->stacktrace;
}
1;

=head1 BUGS

=head1 SUPPORT

=head1 AUTHOR

    Feyruz Yalcin
    CPAN ID: FYALCIN

=head1 SEE ALSO

=cut

package Pathogens::Variant::Exception::Argument;
use base 'Pathogens::Variant::Exception';
1;

package Pathogens::Variant::Exception::CommandExecution;
use base 'Pathogens::Variant::Exception';
1;

package Pathogens::Variant::Exception::Format;
use base 'Pathogens::Variant::Exception';
1;

package Pathogens::Variant::Exception::File;
use base 'Pathogens::Variant::Exception';
1;

package Pathogens::Variant::Exception::Default;
use base 'Pathogens::Variant::Exception';
1;

package Pathogens::Variant::Exception::ObjectUsage;
use base 'Pathogens::Variant::Exception';
1;

