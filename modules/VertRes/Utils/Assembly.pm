=head1 NAME

VertRes::Utils::Assembly - wrapper for calling assembler and assembly optimiser

=head1 SYNOPSIS

=head1 DESCRIPTION

Provides a uniform interface for assembly so you dont need to worry about the differences between assemblers

=head1 AUTHOR

path-help@sanger.ac.uk

=cut

package VertRes::Utils::Assembly;

use strict;
use warnings;
use VertRes::IO;
use File::Basename;
use VertRes::Utils::Hierarchy;
use VertRes::Utils::FileSystem;
use VRTrack::Lane;

use base qw(VertRes::Base);

=head2 new

 Title   : new
 Usage   : my $obj = VertRes::Utils::Assembly->new();
 Function: Create a new VertRes::Utils::Assembly object.
 Returns : VertRes::Utils::Assembly object
 Args    : assember => 'velvet'|'abyss' 

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(@args);

    return $self;
}

=head2 find_module

 Title   : find_module
 Usage   : my $module = $obj->find_module();
 Function: Find out what assembly utility module to use
 Returns : class string (call new on it)

=cut

sub find_module {
    my ($self) = @_;
    return "VertRes::Utils::Assemblers::".$self->{assembler};
}


=head2 generate_files_str

 Title   : generate_files_str
 Usage   : my $module = $obj->generate_files_str();
 Function: create the input string for the files to go into the assembler
 Returns : string which can be passed into assembler

=cut
sub generate_files_str
{
  my $self = shift;
  $self->throw("This is supposed to be overriden");
}



=head2 estimate_memory_required

 Title   : estimate_memory_required
 Usage   : my $memory_required_in_kb = $obj->estimate_memory_required();
 Function: estimate the momory required for the assembler in KB
 Returns : integer in kb of memory requirement

=cut
sub estimate_memory_required
{
  my $self = shift;
  $self->throw("This is supposed to be overriden");
}


=head2 assembly_directories

 Title   : assembly_directories
 Usage   : my $module = $obj->assembly_directories();
 Function: Find out where the assemlbies are located
 Returns : array of paths

=cut
sub assembly_directories
{
  my $self = shift;
  $self->throw("This is supposed to be overriden");
}


=head2 generate_stats

 Title   : generate_stats
 Usage   : my $module = $obj->generate_stats($directory);
 Function: Generate stats for each fa file in the directory

=cut
sub generate_stats
{
  my ($self, $directory) = @_;
  my @output_files;
  opendir(DIR, $directory);
  my @files = grep {/\.fa$/} readdir(DIR);
  for my $file(@files)
  {
    next unless(-e $directory.'/'.$file);
    system("stats $directory/$file > $directory/$file.stats");
    system("/software/pathogen/external/bin/seqstat $directory/$file > $directory/$file.seqstat");
    # Use the GC plots from the QC pipeline instead
    system("~mh12/git/python/fastn2gc.py $directory/$file $directory/$file.png");
  }
}




1;
