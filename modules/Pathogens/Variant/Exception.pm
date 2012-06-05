package Pathogens::Variant::Exception;
use base qw(Error);
use overload ('""' => 'stringify', fallback => 1); #this is what we do if called in string context.

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

The CPAN Error.pm (http://search.cpan.org/~shlomif/Error-0.17016/lib/Error.pm) is the base class for this object.

=head1 USAGE

=cut

sub new {

    my ($proto, $params) = @_;

    my $class = ref($proto) || $proto;
    my $text  = defined $$params{text} ? $$params{text} :   ' ';
    my $value = defined $$params{value} ? $$params{value} : ' ';
    local $Error::Depth = $Error::Depth + 1;
    local $Error::Debug = 1;                   # Enables storing of stacktrace
    my $exception = $class->SUPER::new( -text => $text, -value => $value );
    return $exception;

}

sub stringify {

    my $self = shift;
    my $class = ref($self) || $self;
    my $description = $class . "=". $self->text . "\n" . $self->value;
    return $description;
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

package Pathogens::Variant::Exception::NotImplemented;
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

