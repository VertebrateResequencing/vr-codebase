=head1 NAME

VertRes::Wrapper::tophat - wrapper for tophat

=head1 SYNOPSIS

[stub]

=head1 DESCRIPTION

bowtie2-build ref.fa out_base

tophat ref read_1.fastq read_2.fastq

=head1 AUTHOR

Sendu Bala: bix@sendu.me.uk

=cut

package VertRes::Wrapper::tophat;

use strict;
use warnings;
use File::Copy;
use VertRes::IO;
use VertRes::Parser::fastqcheck;
use VertRes::Wrapper::fastqcheck;
use VertRes::Wrapper::bowtie2;
use File::Temp qw/ tempdir /;
use File::Basename;
use IO::File;

use base qw(VertRes::Wrapper::MapperI);

my %e_parameter = (
    37 	=> 70,
    63	=> 120,
    76  => 140,
    92 	=> 160,
    108 => 160,
    123	=> 200,
    156	=> 240
);

=head2 new

 Title   : new
 Usage   : my $wrapper = VertRes::Wrapper::tophat->new();
 Function: Create a VertRes::Wrapper::tophat object.
 Returns : VertRes::Wrapper::tophat object
 Args    : quiet   => boolean

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = $class->SUPER::new(exe => 'tophat', @args);
    $self->{orig_exe} = $self->exe;
    
    return $self;
}

=head2 version

 Title   : version
 Usage   : my $version = $obj->version();
 Function: Returns the program version.
 Returns : string representing version of the program 
 Args    : n/a

=cut

sub version {
  my $self = shift;

  my $exe = $self->{orig_exe};
  open(my $fh, "$exe -v 2>&1 |") || $self->throw("Could not start $exe");
  my $version = 0;
  while (<$fh>) {
      if (/TopHat v(\S+)/) {
          $version = $1;
          last;
      }
  }
  close($fh);

  return $version;
}

=head2 setup_reference

 Title   : setup_reference
 Usage   : $obj->setup_reference($ref_fasta);
 Function: Do whatever needs to be done with the reference to allow mapping.
 Returns : boolean
 Args    : n/a

=cut

sub setup_reference {
    my ($self, $ref) = @_;
    
    my $bowtie2 = VertRes::Wrapper::bowtie2->new();
    return $bowtie2->setup_reference($ref);
}

=head2 setup_fastqs

 Title   : setup_fastqs
 Usage   : $obj->setup_fastqs($ref_fasta, @fastqs);
 Function: Do whatever needs to be done with the fastqs to allow mapping.
 Returns : boolean
 Args    : n/a

=cut

sub setup_fastqs {
    my ($self, $ref, @fqs) = @_;
    
    foreach my $fq (@fqs) {
        if ($fq =~ /\.gz$/) {
            my $fq_new = $fq;
            $fq_new =~ s/\.gz$//;
            
            unless (-s $fq_new) {
                my $i = VertRes::IO->new(file => $fq);
                my $o = VertRes::IO->new(file => ">$fq_new");
                my $ifh = $i->fh;
                my $ofh = $o->fh;
                while (<$ifh>) {
                    print $ofh $_;
                }
                $i->close;
                $o->close;
            }
        }
    }
    
    return 1;
}

=head2 generate_sam

 Title   : generate_sam
 Usage   : $obj->generate_sam($out_sam, $ref_fasta, @fastqs);
 Function: Do whatever needs to be done with the reference and fastqs to
           complete mapping and generate a sam/bam file.
 Returns : boolean
 Args    : n/a

=cut

sub generate_sam {
    my ($self, $out, $ref, $fq1, $fq2, %other_args) = @_;
    my @fqs = ($fq1, $fq2);

    my($filename, $base_directory, $suffix) = fileparse($out);
    my $directories = tempdir( DIR => $base_directory );

    unless (-s $out) {
        my $e = 70;

        my $longest_read = 0;
        foreach my $fq (@fqs) {

            # create a fastqcheck file if it doesnt already exist
            unless(-e "$fq.fastqcheck")
            {
              my $fqc = VertRes::Wrapper::fastqcheck->new();
              $fqc->run($fq, "$fq.fastqcheck");
            }

            my $pars = VertRes::Parser::fastqcheck->new(file => "$fq.fastqcheck");
            my $length = $pars->max_length();
            if ($length > $longest_read) {
                $longest_read = $length;
            }
        }

        foreach (sort { $a <=> $b } keys %e_parameter) {
            if ($longest_read <= $_) {
                $e = $e_parameter{$_};
            }
        }

        foreach my $fq (@fqs) {
            $fq =~ s/\.gz$//;
        }

        my $X = 1000;
        my $max_intron_length_str = $self->_max_intron_length_option(%other_args);
        my $inner_mate_str = $self->_insert_size_option($longest_read, %other_args);
        my $min_intron_length_str = $self->_min_intron_length_option(%other_args);
        my $max_multihits_str = $self->_max_multihits_option(%other_args);
	my $library_type = $self->_library_type_option(%other_args);

        if((defined $other_args{is_paired}) && $other_args{is_paired} == 0)
        {
          # single ended
          $self->simple_run("$library_type $max_intron_length_str $min_intron_length_str -o $directories $max_multihits_str --no-convert-bam $ref $fqs[0]");
        }
        else
        {
          #paired_ended
          $self->simple_run("$library_type $max_intron_length_str $min_intron_length_str $inner_mate_str $max_multihits_str -o $directories --no-convert-bam $ref $fqs[0] $fqs[1]");
        }

	if ( -e "$directories/unmapped.sam" || -e "$directories/unmapped.bam" ) {

	  $self->add_unmapped($ref, $directories);

	  Utils::CMD( qq[ mv $directories/merged.sorted.fixed.sam $out ] );


	}
	else {

	  Utils::CMD( qq[ mv $directories/accepted_hits.sam $out ] );

	}
    }

    return -s $out ? 1 : 0;
}

=head2 add_unmapped

 Title   : add_unmapped
 Usage   : $obj->add_unmapped($ref, $directories);
 Function: Fix quality flags in the unmapped sam file
	   Merge accepted_hits and unmapped files
	   Sort the merged file by name
	   Apply samtools fixmate to the sorted merge file
 Returns : boolean
 Args    : $ref, $directories

=cut

sub add_unmapped {

    my ( $self, $ref, $directories ) = @_;

    #Operation file names list
    #We need to go through several intermediate steps to complete the merge operation
    #I am assuming only sam files here and no bams
    my %op_files = (
		    accepted_sam => "$directories/accepted_hits.sam",
		    accepted_bam => "$directories/accepted_hits.bam",
		    unmapped_sam => "$directories/unmapped.sam",
		    unmapped_bam => "$directories/unmapped.bam",
		    qfixed_unmapped_sam => "$directories/qfixed_unmapped.sam",
		    qfixed_unmapped_bam => "$directories/qfixed_unmapped.bam",
		    merged_bam => "$directories/merged.bam",
		    merge_sorted_bam => "$directories/merged.sorted.bam",
		    merge_sorted_fixed_bam => "$directories/merged.sorted.fixed.bam",
		    merge_sorted_fixed_sam => "$directories/merged.sorted.fixed.sam"
		   );

    my $samtools_sort_name = qq[ $directories/merged.sorted ];

    my $fix_flag_command = qq[samtools view -h $op_files{unmapped_bam} | ] . q[awk 'BEGIN{FS=OFS="\t"} $5 ~ /255/ { $5 = "0"; }1' > ] . qq[ $op_files{qfixed_unmapped_sam} ];
    my $ff_exit = system($fix_flag_command);

    #We will need to convert the sam files into bams, otherwise the merge won't happen

    my $acc_h_convert_command = qq[ samtools view -bS $op_files{accepted_sam} > $op_files{accepted_bam} ];
    my $acchc_exit = system($acc_h_convert_command);

    my $unmap_convert_command = qq[ samtools view -bS $op_files{qfixed_unmapped_sam} > $op_files{qfixed_unmapped_bam} ];
    my $uc_exit = system($unmap_convert_command);

    #This might be needed if files are huge. I'm not sure.
    #Let's wait for the files to be converted, in case perl decides to move on before file creation
    #while (!-e $op_files{accepted_bam} ) {
    #  sleep(1);
    #}
    #while (!-e $op_files{qfixed_unmapped_bam} ) {
    #  sleep(1);
    #}

    #We can now merge the accepted_hits and the unampped bam files into the merged bam file
    my $merge_command = qq[ samtools merge -f $op_files{merged_bam} $op_files{accepted_bam} $op_files{qfixed_unmapped_bam} ];
    my $m_exit = system($merge_command);


    #Once merged, we will need to sort the merged file by read name. samtools fixmate needs a sam/bam file sorted by read name

    Utils::CMD( qq[ samtools sort -n $op_files{merged_bam} $samtools_sort_name ] );

    #Once sorted, we can now run samtools fixmate on it. We will use the most recent version of fixmate (hence the full path to it's location)
    #This should be updated once the new verison of samtools is released and available in /software
    my $most_recent_samtools_fixmate_exe = qq[ /software/vertres/bin-external/samtools-exp-rc fixmate ];

    Utils::CMD( qq[ $most_recent_samtools_fixmate_exe $op_files{merge_sorted_bam} $op_files{merge_sorted_fixed_bam} ] );

    #To keep consistency across the vr-codebase, we convert the final bam file to sam


    my $convert_bam_to_sam = qq[ samtools view -h $op_files{merge_sorted_fixed_bam} > $op_files{merge_sorted_fixed_sam} ];
    my $cbts_exit = system($convert_bam_to_sam);

    return 1;
}



sub _library_type_option
{
  my ($self, %other_args) = @_;
  
  if(defined $other_args{additional_mapper_params})
  {
     $other_args{library_type} = $1 if $other_args{additional_mapper_params} =~ m/--library-type\s+([\w-]+)/;
  }
  if(defined $other_args{library_type} ) {

    return "--library-type ".$other_args{library_type};
  }
  else {

    return '';
  }
}

sub _max_intron_length_option
{
  my ($self, %other_args) = @_;
  
  $other_args{max_intron_length}  ||=10000;
  
  if(defined $other_args{additional_mapper_params})
  {
     $other_args{max_intron_length} = $1 if $other_args{additional_mapper_params} =~ m/-I\s*(\d+)/;
  }
  
  return "-I ".$other_args{max_intron_length};
}

sub _insert_size_option
{
  my ($self, $longest_read, %other_args) = @_;
  
  $other_args{insert_size}  ||=500;
  my $inner_mate_distance = $other_args{insert_size} - (2*$longest_read);
  if($inner_mate_distance <0)
  {
    $inner_mate_distance = 50;
  }
  return "-r $inner_mate_distance";
}

sub _min_intron_length_option
{
  my ($self, %other_args) = @_;
  my $min_intron_length;
  
  if(defined $other_args{additional_mapper_params})
  {
     $min_intron_length = $1 if $other_args{additional_mapper_params} =~ m/-i\s*(\d+)/;
  }
  
  return defined($min_intron_length) ? "-i $min_intron_length" : '';
}

sub _max_multihits_option
{
  my ($self, %other_args) = @_;
  my $max_multihits = 1; # default to 1
  
  if(defined $other_args{additional_mapper_params})
  {
    $max_multihits = $1 if $other_args{additional_mapper_params} =~ m/-g\s*(\d+)/;
  }
  
  return "-g $max_multihits";
}


=head2 do_mapping

 Title   : do_mapping
 Usage   : $wrapper->do_mapping(ref => 'ref.fa',
                                read1 => 'reads_1.fastq',
                                read2 => 'reads_2.fastq',
                                output => 'output.sam');
 Function: Run mapper on the supplied files, generating a sam file of the
           mapping. Checks the sam file isn't truncated.
 Returns : n/a
 Args    : required options:
           ref => 'ref.fa'
           output => 'output.sam'

           read1 => 'reads_1.fastq', read2 => 'reads_2.fastq'
           -or-
           read0 => 'reads.fastq'

=cut


sub do_mapping {
    my ($self, %args) = @_;
    
    my $orig_run_method = $self->run_method;
    $self->run_method('system');
    
    my $ref_fa = delete $args{ref} || $self->throw("ref is required");
    my @fqs;
    if (defined $args{read1} && defined $args{read2} && ! defined $args{read0}) {
        push(@fqs, delete $args{read1});
        push(@fqs, delete $args{read2});
    }
    elsif (defined $args{read0} && ! defined $args{read1} && ! defined $args{read2}) {
        push(@fqs, delete $args{read0});
    }
    else {
        $self->throw("Bad read args: need read1 and read2, or read0");
    }
    my $out_sam = delete $args{output};
    
    # setup reference-related files
    $self->setup_reference($ref_fa) || $self->throw("failed during the reference step");
    
    # setup fastq-related files (may involve alignment)
    $self->setup_fastqs($ref_fa, @fqs) || $self->throw("failed during the fastq step");
    
    # run the alignment/ combine alignments into final sam file
    unless (-s $out_sam) {
        my $tmp_sam = $out_sam.'_tmp';
        
        if (@fqs == 2) {
            unless (-s $tmp_sam) {
               $self->generate_sam($tmp_sam, $ref_fa, @fqs, %args) || $self->throw("failed during the alignment step");
            }
        }
        else {
            unless (-s $tmp_sam) {
               $self->generate_sam($tmp_sam, $ref_fa, $fqs[0], undef, %args) || $self->throw("failed during the alignment step");
            }
        }
        
        # add in unmapped reads if necessary
        #$self->add_unmapped($tmp_sam, @fqs) || $self->throw("failed during the add unmapped step");

        # tophat throws away some reads
        move($tmp_sam, $out_sam) || $self->throw("Failed to move $tmp_sam to $out_sam: $!");
        $self->_set_run_status(2);
    }
    else {
        $self->_set_run_status(1);
    }
    
    $self->run_method($orig_run_method);
    return;
}



=head2 run

 Title   : run
 Usage   : Do not call directly: use one of the other methods instead.
 Function: n/a
 Returns : n/a
 Args    : paths to input/output files

=cut

1;
