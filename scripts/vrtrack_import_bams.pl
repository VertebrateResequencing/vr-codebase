#!/usr/bin/env perl
# 
# Usage: vrtrack_import_bams.pl --index <bam.index> --database <database.name> --do_updates \\
#        --local_bam_dir <path.to.local.bam.dir> --verbose
# 
# Description: Update vrtrack database with BAMs on local disk from a non-tracked source. 
#     Meta data will be read from the BAM headers unless explicitly overridden 
#     in the index file or with a command line option.
#     
#     Index file is a tab delimited file which must contain the first two columns 
#     defining the md5 and BAM file. Other columns are optional and will override 
#     any values in the BAM headers. Run the script without the do_updates option to 
#     get a template index file which may be edited.
#     
#         [0]  md5 (required)
#         [1]  file (required)
#         [2]  sample (optional) Will override value in \@RG SM tag
#         [3]  population (optional) Will override value in \@RG DS tag
#         [4]  project (optional) Will override value in \@RG DS tag
#         [5]  library (optional) Will override value in \@RG LB tag
#         [6]  seq_centre (optional) Will override value in \@RG CN tag
#         [7]  seq_tech (optional) Will override value in \@RG PL tag
#         [8]  lane (optional) Will override value in \@RG ID or \@RG PU tag
#         [9]  insert_size (optional) Will override value in \@RG PI tag
#         [10] alias (optional) Alternate name for your sample that may already be in the db
# 
# Options:
#     --index         <file> (required) BAM index file as described above.
#     --database      <db name> (required) vrtrack database to update.
#     --local_bam_dir <path> (optional) root directory from which 
#     --do_updates    (optional) if this option is not set, the script will 
#                     simply print out the BAM index file that will be used.
#     -v --verbose    (optional) be verbose.
#     --strict        (optional) be strict about updates.
#     -h --help       This message.
# 
#     ** NB: If --local_bam_dir option is set, add this option to your Import 
#     pipeline conf file to ensure the import can locate your BAMs on disk.
# 
# Author: Shane McCarthy, sm15\@sanger.ac.uk
#
use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use VertRes::Parser::bam;
use VertRes::Utils::FileSystem;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Library;
use Time::Format;

my ($index, $database, $local_bam_dir, $do_updates, $verbose, $strict, $help);
GetOptions('index=s'          => \$index,
           'database=s'       => \$database,
           'local_bam_dir=s'  => \$local_bam_dir,
           'do_updates'       => \$do_updates,
           'v|verbose'        => \$verbose,
           'strict'           => \$strict,
           'h|help'           => \$help);

my $missing_opts = 0;
unless ($index && $database) {
    $missing_opts = 1;
}

($help || $missing_opts) and die <<USAGE;

Usage: $0 --index <bam.index> --database <database.name> --do_updates \\
       --local_bam_dir <path.to.local.bam.dir> --verbose

Description: Update vrtrack database with BAMs on local disk from a non-tracked source. 
    Meta data will be read from the BAM headers unless explicitly overridden 
    in the index file or with a command line option.
    
    Index file is a tab delimited file which must contain the first two columns 
    defining the md5 and BAM file. Other columns are optional and will override 
    any values in the BAM headers. Run the script without the do_updates option to 
    get a template index file which may be edited.
    
        [0]  md5 (required)
        [1]  file (required)
        [2]  sample (optional) Will override value in \@RG SM tag
        [3]  population (optional) Will override value in \@RG DS tag
        [4]  project (optional) Will override value in \@RG DS tag
        [5]  library (optional) Will override value in \@RG LB tag
        [6]  seq_centre (optional) Will override value in \@RG CN tag
        [7]  seq_tech (optional) Will override value in \@RG PL tag
        [8]  lane (optional) Will override value in \@RG ID or \@RG PU tag
        [9]  insert_size (optional) Will override value in \@RG PI tag
        [10] alias (optional) Alternate name for your sample that may already be in the db

Options:
    --index         <file> (required) BAM index file as described above.
    --database      <db name> (required) vrtrack database to update.
    --local_bam_dir <path> (optional) root directory from which 
    --do_updates    (optional) if this option is not set, the script will 
                    simply print out the BAM index file that will be used.
    -v --verbose    (optional) be verbose.
    --strict        (optional) be strict about updates.
    -h --help       This message.

    ** NB: If --local_bam_dir option is set, add this option to your Import 
    pipeline conf file to ensure the import can locate your BAMs on disk.

Author: Shane McCarthy, sm15\@sanger.ac.uk

USAGE

my %platform_to_tech = (
    ILLUMINA => 'SLX',
    LS454 => '454',
    ABI_SOLID => 'SOLID',
    SLX => 'SLX',
    '454' => '454',
    SOLID => 'SOLID'
);

my $fsu = VertRes::Utils::FileSystem->new();

my %data;
open my $ifh, "<$index" || die "Could not open index file $index\n";
while (<$ifh>) {
    chomp;
    /\S/ || next;
    /^#/ && next;
    my ($md5, $bam, $sample, $population, $project, $library, $seq_centre, $seq_tech, $lane, $insert_size, $alias, $date) = split /\t/;
    die "Index file must contain at least two columns containing md5 and the bam file" unless ($md5 && $bam);
    die "$md5 does not appear to be a md5sum\n" unless ($md5 =~ /^[0-9a-z]{32}$/);
    die "$bam does not appear to be a bam file\n" unless ($bam =~ /.bam$/);
    my $bam_path = $bam;
    $bam_path = $fsu->catfile($local_bam_dir, $bam) if ($local_bam_dir);
    die "$bam does not exist! Try using the local_bam_dir option if $bam is a relative path\n" unless (-s $bam_path);
    
    my $bp = VertRes::Parser::bam->new(file => $bam_path);
    my %rg_info = $bp->readgroup_info();
    $bp->close;
    my @rgs = keys %rg_info;
    die "'$bam' does not contain a single readgroup\n" unless (scalar @rgs == 1);
    my $rg = shift @rgs;
    
    # Set lane name to @RG->ID unless this is '1', in which case try to set it to @RG->PU
    my $lane_name = ($rg eq '1' && exists($rg_info{$rg}{PU})) ? $rg_info{$rg}{PU} : $rg;
    
    my $run_date = $date || $data{$bam}{DT}  || $time{'yyyy-mm-dd hh:mm:ss'};
    
    $data{$bam}{path}        = $bam_path;
    $data{$bam}{md5}         = $md5;
    $data{$bam}{sample}      = $sample      || $rg_info{$rg}{SM};
    $data{$bam}{population}  = $population  || $rg_info{$rg}{DS};
    $data{$bam}{project}     = $project     || $rg_info{$rg}{DS};
    $data{$bam}{lane}        = $lane        || $lane_name;
    $data{$bam}{library}     = $library     || $rg_info{$rg}{LB};
    $data{$bam}{seq_centre}  = $seq_centre  || $rg_info{$rg}{CN};
    $data{$bam}{seq_tech}    = $seq_tech    || $platform_to_tech{ $rg_info{$rg}{PL} };
    $data{$bam}{insert_size} = $insert_size || $rg_info{$rg}{PI};
    $data{$bam}{alias}       = $alias       || $data{$bam}{sample};
    $data{$bam}{run_date}    = Time::Format::time_format('yyyy-mm-dd hh:mm:ss', $run_date);
}
close $ifh;

unless ($do_updates) {
    foreach my $bam (keys %data) {
        print join "\t", ($data{$bam}{md5}, $bam,
            $data{$bam}{sample}, $data{$bam}{population}, $data{$bam}{project}, 
            $data{$bam}{library}, $data{$bam}{seq_centre}, $data{$bam}{seq_tech}, 
            $data{$bam}{lane}, $data{$bam}{insert_size}, $data{$bam}{alias}, 
            $data{$bam}{run_date});
            print "\n";
    }
    exit;
}

my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $database, mode => 'rw');
unless ($vrtrack) {
    die "DB connection failed: ".$DBI::errstr."\n";
}

my %object_cache;
my $lanes = 0; 
foreach my $bam (keys %data) 
{
    $vrtrack->transaction_start();
    
    my $project_name = $data{$bam}{project};
    my $population_name = $data{$bam}{population};
    my $sample_name = $data{$bam}{sample};
    my $library_name = $data{$bam}{library};
    my $seq_tech_name = $data{$bam}{seq_tech};
    my $seq_centre_name = $data{$bam}{seq_centre};
    my $lane_name = $data{$bam}{lane};
    my $alias = $data{$bam}{alias};
    
    my $md5 = $data{$bam}{md5};
    my $insert_size = $data{$bam}{insert_size};
    my $run_date = $data{$bam}{run_date};
    
    # Create project if necessary
    my $project = $object_cache{$project_name} || $vrtrack->get_project_by_name($project_name);
    unless ($project) {
        $project = $vrtrack->add_project($project_name);
        print "adding project $project_name\n" if $verbose;
        $project || die "Could not initialise project object '$project_name'";
        $project->hierarchy_name($project_name);
        $project->update || die "Could not update project";
    }
    $object_cache{$project_name} = $project;
    
    # Add sample if necessary
    my $sample = $project->get_sample_by_name($sample_name) || $project->get_sample_by_name($alias);
    unless ($sample) {
        $sample = $project->add_sample($sample_name);
        print "adding sample $sample_name\n" if $verbose;
        $sample || die "Could not initialise sample object '$sample_name'";
        $sample->hierarchy_name($sample_name);
        $sample->update;
        $project->update;
    }
    die "sample $sample_name is a Sanger sample and should be imported by other means" if (is_sanger_sample($sample));
    
    # Create individual and population if necessary
    my $individual = $sample->individual();
    unless ($individual) {
        $individual = $sample->individual($sample_name);
        $individual = $sample->add_individual($sample_name) unless $individual;
        print "adding individual $sample_name\n" if $verbose;
        $individual->update;
        $sample->update;
        
        my $population = $individual->population();
        unless ($population)  {
            $population = $individual->population($population_name);
            $individual->add_population($population_name) unless $population;
            print "adding population $population_name\n" if $verbose;
        }
        $population->update;
        $individual->update;
    }
    
    # Add library
    my $library = $sample->get_library_by_name($library_name);
    unless ($library) {
        $library = VRTrack::Library->new_by_name($vrtrack, $library_name);
        if ($library) 
        {
            die "Library name '$library_name' already exists\n";
        }
        else 
        {
            $library = $sample->add_library($library_name);
            print "adding library $library_name\n" if $verbose;
            my $hname = $library_name;
            $hname =~ s/[^\w\d]/_/g;
            $library->hierarchy_name($hname);
            $library->fragment_size_from($insert_size);
            $library->fragment_size_to($insert_size);
            
            # Add seq_tech
            my $seq_tech = $library->seq_tech();
            unless ($seq_tech) {
                $seq_tech = $library->seq_tech($seq_tech_name);
                $seq_tech = $library->add_seq_tech($seq_tech_name) unless $seq_tech;
                print "adding seq_tech $seq_tech_name\n" if $verbose;
                $seq_tech->update;
            }

            # Add seq_centre
            my $seq_centre = $library->seq_centre();
            unless ($seq_centre) {
                $seq_centre = $library->seq_centre($seq_centre_name);
                $seq_centre = $library->add_seq_centre($seq_centre_name) unless $seq_centre;
                print "adding seq_centre $seq_centre_name\n" if $verbose;
                $seq_centre->update;
            }
            
            $library->update;
            $sample->update;
        }
    }
    die "Did not create library object, $bam\n" unless $library;
    if ( $library->insert_size ne $insert_size ) {
        if ($strict) {
            print STDERR "Warning: insert size does not match for library $library_name, ".$library->insert_size." vs $insert_size value, skipping...\n";
            $vrtrack->transaction_rollback();
            next;
        } else {
            print "Warning: insert size does not match for library $library_name, ".$library->insert_size." vs $insert_size value in db will be overridden\n";
        }
    }
    
    # Add lane and file
    my $lane = $library->get_lane_by_name($lane_name);
    if ($lane) {
        my $file = $lane->get_file_by_name($bam);
        unless ($file) {
            $file = $lane->add_file($bam);
            $file->md5($md5);
            $file->hierarchy_name(File::Basename::basename($bam));
            $file->update;
            
            $lane->update;
        }
    } else {
        $lane = $library->add_lane($lane_name);
        print "adding lane $lane_name\n" if $verbose;
        $lane->hierarchy_name($lane_name);
        $lane->run_date($run_date);
        $lane->update;
        ++$lanes;
        
        my $file = $lane->add_file($bam);
        print "adding file $bam\n" if $verbose;
        $file->md5($md5);
        $file->hierarchy_name(File::Basename::basename($bam));
        $file->update;
        
        $lane->update;
    }
    
    $vrtrack->transaction_commit();
}
print "$lanes lanes added\n" if $verbose;

exit;

sub is_sanger_sample {
    my $sample = shift;
    my $libraries = $sample->libraries();
    my $sc = 0;
    foreach my $lib (@$libraries) {
        my $seq_centre = $lib->seq_centre;
        if ($seq_centre && ($seq_centre->name eq 'SC')) {
            $sc = 1;
            last;
        }
    }
    return $sc;
}
