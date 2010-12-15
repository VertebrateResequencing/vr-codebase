package MyBuild;
# copied from Perl Testing: A developer's notebook (O'Reilly)
# http://proquestcombo.safaribooksonline.com/0596100922/130

BEGIN { unshift(@INC, './modules') }
use base 'Module::Build';
use File::Path;
use Data::Dumper;
use File::Spec::Functions;
use VRTrack::VRTrack;

sub create_vrconfig_file {
    my $self     = shift;
    
    my $configure_vrdb = $self->prompt('Do you wish to setup a test db for VRTrack? [y/n]: ', 'y');
    
    if ($configure_vrdb eq 'y') {
        my $config   =
        {   
            host     => $self->prompt( 'VRTrack test DB host: ', $ENV{VRTRACK_HOST} || 'localhost' ),
            port     => $self->prompt( 'VRTrack test DB port: ', $ENV{VRTRACK_PORT} || 3306 ),
            user     => $self->prompt( 'VRTrack test DB user: ', $ENV{VRTRACK_RW_USER} || '' ),
            password => $self->prompt( 'VRTrack test DB password: ', $ENV{VRTRACK_PASSWORD} || '' ),
            test_db  => $self->prompt( 'VRTrack test DB name: ', 'vrtrack_test' ),
        };
        $self->notes( db_config    => $config );
        mkpath( catdir( qw( modules VRTrack  ) ) );
        my $dd       = Data::Dumper->new( [ $config ], [ 'db_config' ] );
        my $path     = catfile(qw( modules VRTrack Testconfig.pm ));
        open( my $file, '>', $path ) or die "Cannot write to '$path': $!\n";
        printf $file <<'END_HERE', $dd->Dump();
package VRTrack::Testconfig;
my $db_config;
%s;

sub config
{
    my ($self, $key) = @_;
    return $db_config->{$key} if exists $db_config->{$key};
}
1;
END_HERE
    }
}

sub create_database {
    my ($self, $dbname) = @_;
    my $config          = $self->notes( 'db_config' );
    if ($config) {
        my @sql = VRTrack::VRTrack->schema();
        open(my $mysqlfh, "| mysql -h$config->{host} -u$config->{user} -p$config->{password} -P$config->{port}") || die "could not connect to database for testing\n";
        print $mysqlfh "drop database if exists $config->{test_db};\n";
        print $mysqlfh "create database $config->{test_db};\n";
        print $mysqlfh "use $config->{test_db};\n";
        foreach my $sql (@sql) {
            print $mysqlfh $sql;
        }
        close($mysqlfh);
    }
}

sub ACTION_test {
    my $self   = shift;
    my $config = $self->notes( 'db_config' );
    if ($config) {
        eval {
            $self->create_database( $config->{test_db} );
        };
        if ($@){
            print $@;
        }
        else {
            print "Created database $config->{test_db}\n";
        }
    }
    $self->SUPER::ACTION_test( @_ );
}
1;

