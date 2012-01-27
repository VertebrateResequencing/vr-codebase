=head1 NAME

BAMStats.pm   - Functionality for Stats files for a BAM

=head1 SYNOPSIS

use Pathogens::RNASeq::BAM;
my $bam_container = Pathogens::RNASeq::BAM->new(
  filename => '/abc/my_file.bam'
  );

=cut
package Pathogens::RNASeq::BAMStats;
use Moose;
use VertRes::Parser::flagstat;
use VertRes::Utils::Sam;
use Time::Format;

has 'total_mapped_reads' => ( is => 'rw', isa => 'Str', lazy_build   => 1 );
has '_parser'  => ( is => 'rw', isa => 'VertRes::Parser::flagstat', lazy_build   => 1 );

sub _build__parser
{
  my ($self) = @_;
  unless(-e $self->filename.'.flagstat')
  {
    $self->_create_stats_files();
  }
  
  VertRes::Parser::flagstat->new(file => $self->filename.'.flagstat');
}


sub _build_total_mapped_reads
{
  my ($self) = @_;
  $self->_parser->mapped_reads;
}


sub _create_stats_files
{
   my ($self) = @_;
   my $sam =  VertRes::Utils::Sam->new();
   $sam->stats("$time{'yyyymmdd'}", $self->filename);
}
1;