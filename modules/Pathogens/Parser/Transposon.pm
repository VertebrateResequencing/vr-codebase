=head1 NAME

Transposon.pm - Parse a fastq.gz file to find the percentage of tags that are present in reads

=head1 SYNOPSIS

use Pathogens::Parser::Transposon;
my $tradis_transposon = Pathogens::Parser::Transposon->new(
  filename   => 'myfile.fastq.gz',
  tag_length => 7,
);
$tradis_transposon->percentage_reads_with_tag ;

=head1 DESCRIPTION

Calculate the percentage of tags that have mapped to reads in a fastq file.  
If the tag is unknown, work it out by sampling the file and taking the most frequently occuring tag. Only tags which are 100% present are counted.

=head1 METHODS

=head2 new

  Arg [1]    : filename   => fastq.gz file
               tag_length => integer with the tag length
               tag => user provided tag (optional)
               num_reads_to_sample => if a tag isnt provided, subsample this number of reads to find the most common tag, defaults to 1000, 0 samples whole file
  Example    :   my $tradis_transposon = Pathogens::Parser::Transposon->new(
                 filename   => 'myfile.fastq.gz',
                 tag_length => 7,
                 tag => 'ACGT');
  Description: returns Transposon object
  Returntype : Pathogens::Parser::Transposon object


=head2 percentage_reads_with_tag

  Arg [1]    : None
  Example    : my percentage_reads = $tradis_transposon->percentage_reads_with_tag ;
  Description: Find out the percentage of reads in a file with a common tag
  Returntype : Float

=cut

package Pathogens::Parser::Transposon;

use Moose;

has 'filename'                  => ( is => 'rw', isa => 'Str', required => 1 );
has 'tag_length'                => ( is => 'rw', isa => 'Int', required => 1 );
has 'tag'                       => ( is => 'rw', isa => 'Maybe[Str]', predicate => 'has_tag',);

has 'inferred_tag'              => ( is => 'rw', isa => 'Str', lazy_build => 1 );
has 'num_reads_to_sample'       => ( is => 'rw', isa => 'Int', default => 1000);
has 'percentage_reads_with_tag' => ( is => 'rw', isa => 'Str', lazy_build => 1 );


sub _build_percentage_reads_with_tag
{
  my $self = shift;
  my $tag = $self->has_tag ? $self->tag : $self->inferred_tag;
  my $count = 0;
  my $tag_count = 0;
  
  open(INPUT_FILE, '-|', 'gzcat', $self->filename) or die 'Couldnt open input file';
  while(<INPUT_FILE>)
  {
    my $line = $_;
    if(($count % 4 == 1) && $line =~ m/^$tag/ )
    {
      $tag_count++;
    }
    $count++;
  }
  close(INPUT_FILE);

  return ($tag_count/($count/4))*100;
}

sub _build_inferred_tag
{
   my $self = shift;
   my %tag_frequency;
   my $count = 0;
   
   open(INPUT_FILE, '-|', 'gzcat', $self->filename) or die 'Couldnt open input file';
   while(<INPUT_FILE>)
   {
     my $line = $_;
     if($count % 4 == 1)
     {
       my $tag = substr($line,0, $self->tag_length);

       if(defined $tag_frequency{$tag})
       {
         $tag_frequency{$tag}++;
       }
       else
       {
         $tag_frequency{$tag} = 1;
       }
     }
     last if(($self->num_reads_to_sample > 0) && $count > ($self->num_reads_to_sample *4) );
     
     $count++;
   }
   close(INPUT_FILE);
   my @tag_frequency_keys_sorted_by_value = sort { $tag_frequency{$a} <=> $tag_frequency{$b} } keys %tag_frequency; 
   return pop(@tag_frequency_keys_sorted_by_value );
}


1;
