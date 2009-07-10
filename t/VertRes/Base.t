#!/usr/bin/perl -w
use strict;
use warnings;

BEGIN {
    use Test::Most tests => 66;
    
    use_ok('VertRes::Base');
    use_ok('File::Spec');
}

my $base = VertRes::Base->new(foo => 'bar');
isa_ok $base, 'VertRes::Base';
is $base->{foo}, 'bar', 'auto setting of instance data';

# verbose and warn testing
is $base->verbose, 0, 'default verbose value';
warning_is { $base->warn('simple msg') } { carped => 'simple msg' }, 'verbose 0 warning is simple carp';
is $base->verbose(-1), -1, 'verbose -1 could be set';
warning_is { $base->warn('simple msg') } '', 'no warning when verbose -1';
is $base->verbose(1), 1, 'verbose 1 could be set';
warning_is { $base->warn('simple msg') } { carped => 'simple msg' }, 'verbose 1 warning is a cluck';
is $base->verbose(2), 2, 'verbose 2 could be set';
throws_ok { $base->warn('thrown msg') } qr/thrown msg/, 'verbose 2 warning is a throw';
is $base->verbose(0), 0, 'verbose 0 could be set';
is VertRes::Base::verbose(1), 1, 'verbose 1 could be set globally';
is $base->verbose, 1, 'global verbose applied to instance';
is $base->verbose(2), 1, 'instance set does not override global';
is VertRes::Base::verbose(undef), 0, 'global verbose can be unset';
is $base->verbose(2), 2, 'instance set works again after unsetting global';

# throw and debug
throws_ok { $base->throw('thrown msg') } qr/thrown msg/, 'throw works';
warning_like { $base->debug('debug msg') } qr/debug msg/, 'debug message when verbose > 0';
$base->verbose(0);
warning_is { $base->debug('debug msg') } '', 'no debug message when verbose <= 0';

# register_for_cleanup difficult to test for properly...
can_ok $base, qw(register_for_cleanup);
$base->register_for_cleanup('foo');
$base->register_for_cleanup('verbose');
is_deeply $base->{'_cleanup_methods'}, {verbose => 1}, 'existing methods can be registered for cleanup';
$base->unregister_for_cleanup('verbose');
is_deeply $base->{'_cleanup_methods'}, {}, 'method can be unregistered for cleanup';

# unlinking can be tested, alongside logging
my $tfile1 = File::Spec->catfile('t', 'data', 'vertres_base_unlink_test_1');
my $tfile2 = File::Spec->catfile('t', 'data', 'vertres_base_unlink_test_2');
system("touch $tfile1; touch $tfile2");
$base->register_for_unlinking($tfile1, $tfile2);
is_deeply $base->{'_unlink_files'}, {$tfile1 => 1, $tfile2 => 1}, 'files can be registered for unlinking';
$base->unregister_for_unlinking($tfile2);
is_deeply $base->{'_unlink_files'}, {$tfile1 => 1}, 'files can be unregistered for unlinking';
$base->register_for_unlinking($tfile2);
ok -e $tfile1, 'file1 exists prior to unlink attept';
ok -e $tfile2, 'file2 exists prior to unlink attept';
my $default_log_file = $base->log_file();
like $default_log_file, qr/vertres.log$/, 'default log file has correct name';
is $base->log_file($tfile1), $tfile1, 'log file location could be changed';
is $base->write_logs, 0, 'logging off by default';
$base->log('message');
ok ! -s $tfile1, 'logging while logging is off does nothing';
is $base->write_logs(1), 1, 'logging could be turned on';
$base->log('message');
my $log_size = -s $tfile1;
cmp_ok $log_size, '>', 0, 'logging while logging is on does something';
$base->log('message2');
my $log_size2 = -s $tfile1;
cmp_ok $log_size2, '>', $log_size, 'adding a log message appends';
undef $base;
ok ! -e $tfile1, 'file1 unlinked after object destruction';
ok ! -e $tfile2, 'file1 unlinked after object destruction';

# test inheritance and auto-setting of methods
package VertRes::Foo;
use base 'VertRes::Base';
sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self->_set_from_args(\@args, methods => ['bar']);
    return $self;
}
sub bar {
    my ($self, $val) = @_;
    if ($val) {
        $self->{_bar} = $val;
    }
    return $self->{_bar};
}
sub foo {
    my ($self, $val) = @_;
    if ($val) {
        $self->{_foo} = $val;
    }
    return $self->{_foo};
}
package main;
my $foo = VertRes::Foo->new(verbose => 1, bar => 'monkey', foo => 'coconut');
is $foo->verbose, 1, 'verbose set via _set_from_args()';
is $foo->bar, 'monkey', 'method set automatically via _set_from_args()';
is $foo->foo, undef, 'foo not set automatically since not asked for';
is $foo->{foo}, 'coconut', 'foo variable was set automatically via _set_from_args()';

# more _set_from_args tests (from BioPerl Root::RootI tests, modified)
package VertRes::Foo1;
use base qw(VertRes::Base);
sub new {
    my $class = shift;
    my $self = {};
    bless $self, ref($class) || $class;
    $self->_set_from_args(\@_);
    return $self;
};
package main;
my $obj = VertRes::Foo1->new(-verbose => 1, t1 => 1, '--Test-2' => 2);
ok ! $obj->can('t1'), 'arg not callable';

package VertRes::Foo2;
use base qw(VertRes::Base);
sub new {
    my $class = shift;
    my $self = {};
    bless $self, ref($class) || $class;
    $self->_set_from_args(\@_, create => 1);
    return $self;
};
package main;
$obj = VertRes::Foo2->new(-verbose => 1, t3 => 1, '--Test-4' => 2);
can_ok $obj, ('t3'); # 'arg callable since method was created'
is $obj->t3, 1, 'created method was set';
can_ok $obj, ('test_4'); # 'mal-formed arg callable since method was created with good name'
for my $m (qw(t3 test_4)) {
    can_ok 'VertRes::Foo2', ($m);
    ok ! UNIVERSAL::can('VertRes::Base', $m), "Methods don't pollute original VertRes::Base namespace";
}

package VertRes::Foo3;
use base qw(VertRes::Base);
sub new {
    my $class = shift;
    my $self = {};
    bless $self, ref($class) || $class;
    $self->_set_from_args(\@_, methods => ['verbose', 't5'], create => 1);
    return $self;
};
package main;
$obj = VertRes::Foo3->new(-verbose => 1, t5 => 1, '--Test-6' => 2);
can_ok $obj, ('t5');
ok ! $obj->can('test_6'), 'arg not in method list not created';
can_ok 'VertRes::Foo3', ('t5');
ok (!UNIVERSAL::can('VertRes::Base','t5'), "Methods don't pollute original VertRes::Base namespace");

package VertRes::Foo4;
use base qw(VertRes::Base);
sub new {
    my $class = shift;
    my $self = {};
    bless $self, ref($class) || $class;
    my %args = @_;
    $self->_set_from_args(\%args, methods => {(verbose => 'v',
                                               test7 => 't7',
                                               test_8 => 't8')},
                                               create => 1);
    return $self;
};
# with synonyms
package main;
$obj = VertRes::Foo4->new(-verbose => 1, t7 => 1, '--Test-8' => 2);
is $obj->verbose, 1, 'verbose was set correctly';
is $obj->t7, 1, 'synonym was set correctly';
is $obj->test7, 1, 'real method of synonym was set correctly';
is $obj->test_8, 2, 'mal-formed arg correctly resolved to created method';
is $obj->t8, 2, 'synonym of set method was set correctly';
for my $m (qw(t7 test7 test_8 t8)) {
    can_ok 'VertRes::Foo4', $m;
    ok ! UNIVERSAL::can('VertRes::Base', 't7'), "Methods don't pollute original VertRes::Base namespace";
}

exit;
