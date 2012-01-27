=head1 NAME

BAS.pm   - Functionality for BAS files for a BAM

=head1 SYNOPSIS

use Pathogens::RNASeq::BAM;
my $file_meta_data_container = Pathogens::RNASeq::BAM->new(
  filename => '/abc/my_file.bam'
  );

=cut
package Pathogens::RNASeq::BAS;
use Moose;
use VertRes::Parser::bas;
use VertRes::Utils::Sam;

has 'total_mapped_reads' => ( is => 'rw', isa => 'Str', lazy_builder   => 1 );
has '_bas_parser'  => ( is => 'rw', isa => 'VertRes::Parser::bas', lazy_builder   => 1 );

sub _build__bas_parser
{
  my ($self) = @_;
  unless(-e $self->filename.'.bas')
  {
    $self->_create_stats_files();
  }
  
  VertRes::Parser::bas->new(file => $self->filename);
}


sub _build_total_mapped_reads
{
  my ($self) = @_;
  $self->_bas_parser->mapped_reads;
}


sub _create_stats_files
{
   my ($self) = @_;
   my $sam =  VertRes::Utils::Sam->new();
   $self->bas($self->filename, "$time{'yyyymmdd'}", $self->filename.'.tmp');
}