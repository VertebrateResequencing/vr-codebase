package Sfind::Sfind;
=head1 NAME

Sfind::Sfind - API for the denormalised SeqScape & NPG tracking database

=head1 SYNOPSIS
    my $track = Sfind::Sfind->new();

    #get arrayref of projects being tracked for traversing hierarchy
    my $projects = $track->projects();

    #also provides accessors for arbitrary objects in hierarchy
    my $lane = $track->get_lane_by_id

=head1 DESCRIPTION

Retrieves data from the sequencing tracking database.

=head1 CONTACT

jws@sanger.ac.uk

=head1 METHODS

=cut


use DBI;
use Sfind::Project;


=head2 new

  Arg [1]    : None
  Example    : my $track = Sfind::Sfind->new()
  Description: Returns Sfind object if can connect to database
  Returntype : Sfind::Sfind object

=cut

sub new {
    my ($class) = @_;

    my $self = {};
    bless ($self, $class);
    my $dbh = DBI->connect("DBI:mysql:host=mcs2a:port=3313;database=warehouse_production", "warehouse_ro",undef,
			     {'RaiseError' => 1, 'PrintError'=>0});
    $self->{_dbh} = $dbh;

    return $self;
}


=head2 get_project_by_id

  Arg [1]    : project id from sequencescape
  Example    : my $project = $track->get_project_by_id(140);
  Description: retrieve project object by sequencescape id
  Returntype : Sfind::Project object

=cut

sub get_project_by_id {
    my ($self, $id) = @_;
    my $obj = Sfind::Project->new($self->{_dbh},$id);
    return $obj;
}


=head2 get_project_by_name

  Arg [1]    : project name in sequencescape
  Example    : my $project = $track->get_project_by_name('1000Genomes-A1-CEU');
  Description: retrieve project object by sequencescape name
  Returntype : Sfind::Project object

=cut

sub get_project_by_name {
    my ($self, $name) = @_;
    my $sql = qq[select project_id from project_information where project_name=?];
    my $id_ref = $self->{_dbh}->selectrow_hashref($sql, undef, ($name));
    unless ($id_ref){
	warn "No project with name $name\n";
	return undef;
    }

    my $id = $id_ref->{project_id};
    return $self->get_project_by_id($id);
}


=head2 project_names

  Arg [1]    : None
  Example    : my $project_names = $track->project_names();
  Description: Returns a ref to an array of the project names that are being tracked
  Returntype : ref to array of project name strings

=cut

sub project_names {
    my ($self) = @_;

    unless ($self->{'project_names'}){
	my $sql = qq[select distinct project_name from project_information];
	my @projects;
	my $sth = $self->{_dbh}->prepare($sql);

	$sth->execute();
	foreach(@{$sth->fetchall_arrayref()}){
	    push @projects, $_->[0];
	}
	$self->{'project_names'} = \@projects;
    }

    return $self->{'project_names'};
}


=head2 get_lane_by_filename

  Arg [1]    : file name
  Example    : my $lane = $track->get_lane_by_filename('123_s_1.fastq');
  Description: retrieve lane object by filename
  Returntype : Sfind::Lane object

=cut

sub get_lane_by_filename {
    my ($self, $name) = @_;
    my $file = $self->get_file_by_name($name);
    my $obj;
    if ($file){
	my $lane_id = $file->lane_id;
	$obj = Sfind::Lane->new($self->{_dbh},$id);
    }
    return $obj;
}

1;
