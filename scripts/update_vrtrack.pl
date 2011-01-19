#!/usr/bin/env perl

=head1 NAME

update_vrtrack.pl 

=head1 SYNOPSIS

=head1 DESCRIPTION

This script updates the specified tracking database using details from the warehouse.
Please note that references to VRTrack 'projects' are actually studies in the warehouse.

=head1 CONTACT

jws@sanger.ac.uk

=cut

use strict;
use warnings;
no warnings 'uninitialized';
#use lib ".";

use Getopt::Long;
use Sfind::Sfind;
use VertRes::Utils::VRTrackFactory;
use VRTrack::Lane;
use VRTrack::Library;
use VRTrack::Individual;
use VRTrack::Multiplex_pool;
use VRTrack::Library_Multiplex_pool;

my ($projfile, $spp, $update_files, $samplemap, $help, $database,$create_individuals);

GetOptions(
    'studies|p|projects=s'  =>  \$projfile,
    's|spp=s'       =>  \$spp,
    'f|files'       =>  \$update_files,
    'd|database=s'  =>  \$database,
    'c|create_individuals'          =>  \$create_individuals,
    'm|sample_map=s'=>  \$samplemap,
    'h|help'        =>  \$help,
    );

my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track',
                   'g1k'    => 'g1k_track',
                 );

# For testing
#my %db_for_spp = ( 'mouse'  => 'mouse_reseq_track_test',
#                    'g1k'   => 'g1k_track_test',
#          );

my $db;
if( $database )
{
    $db = $database;
}
else
{
    $db = $db_for_spp{$spp};
}

($projfile && $db && !$help) or die <<USAGE;
    Usage: $0   
                --studies  <study name or file of SequenceScape study names>
                --spp       <species, i.e. g1k or mouse>
                [--database  <vrtrack database name override>]
                [--files     <force update of files (usually skipped as doesn't change, and is slow)>]
                [--create_individuals  <if set, generates an individual for each new sample name>]
                [--sample_map  <a file of individual -> samplename mappings. Cannot be used with --create_individuals>]
                --help      <this message>

Updates vrtrack databases from sequencescape.

One of spp or database must be supplied.

project file is
project name
acc

The information for new samples (sex, accession, etc) comes from the individual
table - this must be populated for new samples before update will load the
samples.  The exception is if --create_individuals is passed, in which case
each new sample name will trigger creation of a corresponding individual (sex
will be set to 'unknown').

By default, samples will be mapped to existing individuals by the prefix of the
sample and individual name up to the first non-word or underscore character
(e.g. NOD_sample_5 and NOD_brain_tissue will both be mapped to individual
NOD_mouse).  To override this behaviour, --sample_map takes a filename of
tab-separated individual, sample names which explicitly sets the mapping.  This
cannot be used in conjunction with --create_individuals.

USAGE

print "Database: $db\n";
my $vrtrack = VertRes::Utils::VRTrackFactory->instantiate(database => $db,
                                                          mode => 'rw');

unless ($vrtrack){
    die "Can't connect to tracking database\n";
}

my %projects;
if (-s $projfile){
    open (my $PROJ, "$projfile") or die "Can't open $projfile: $!\n";
    while (<$PROJ>){
        chomp;
        my ($project, $acc) = split "\t", $_;
        $projects{$project} = $acc;
    }
    close $PROJ;
}
else {
    $projects{$projfile} = undef;
}

# Set up mapping from samples to individuals by name

my %individual_names; # lookup of individual names by sample key

if ($samplemap){
    die "Can't use --sample_map and --create_individuals together" if $create_individuals;
    open (my $MAPFILE, "<$samplemap") or die "Can't open sample_map $samplemap: $!\n";
    while (<$MAPFILE>){
        chomp;
        my ($ind, $samp) = split "\t", $_;
        $individual_names{$samp} = $ind;
    }
    close $MAPFILE;
}
else {
    foreach my $ind_name(@{$vrtrack->individual_names}){
        if ($create_individuals){
            # just hash on the individual name - we only want to match full names
            if (exists $individual_names{$ind_name}){ 
                die "Individual name collision - $ind_name in DB more than once";
            }
            $individual_names{$ind_name} = $ind_name;
        }
        else {
            # trim off everything after a non-alphanum this is to match mouse
            # samples - g1k are more easily related to individuals
            my $lookup_key = uc($ind_name);
            $lookup_key =~ s/[\W_].*//; 

            if (exists $individual_names{$lookup_key}){ 
                die "Individual lookup collision - $lookup_key from $ind_name matches $lookup_key from ",$individual_names{$lookup_key},"\n";
            }
            $individual_names{$lookup_key} = $ind_name;
        }
    }
}

# important that these are checked in descending length - i.e. more specific
# matches first, down to less specific.  Hopefully this will catch the mouse
# samples like A/J which will be keyed as 'A'
my @individual_lookups = sort{length($b) <=> length($a)} keys %individual_names;

my $strack = Sfind::Sfind->new();

foreach my $pname (keys %projects){
    my $study = $strack->get_study_by_name($pname);
    unless ($study){
        print "Can't get warehouse study $pname, skipping\n";
        next;
    }
    my $vproject = get_core_object('project',$vrtrack,$study);
    unless ($vproject){
        print STDERR "Can't get tracking project for ",$study->name,"; skipping\n";
        next;
    }

    #Get and set the EBI accession number
    my $study_acc = $study->get_accession();

    if ($study_acc){
        my $vstudy = $vproject->study($study_acc);
        unless ($vstudy){
            $vstudy = $vproject->add_study($study_acc);
        }
    }
    print "Project $pname updating\n" if $vproject->dirty;
    $vproject->update;

    my $samples;
    eval {
        $samples = $study->samples;
    };
    if ($@){
        print STDERR "Error getting samples off $pname : ",$@,". Skipping\n";
        next;
    }

    foreach my $sample(@$samples){
        # Check we can map this sample name to an individual
        my $vind_name;  # full individual name
        if ($create_individuals){

            # match only on full sample name - we want one individual/sample
            my $sample_name = $sample->name;
            $vind_name = $individual_names{$sample_name};

            if (! $vind_name){
                # no matching individual.
                # make new individual for this sample
                my $vind = VRTrack::Individual->new_by_name($vrtrack,$sample_name);
               
                next if $vind; # already there, so skip out. Shouldn't happen!

                print "Adding individual: $sample_name.\n";
                $vind = VRTrack::Individual->create($vrtrack, $sample_name);
                # default sex to M - QC pipeline needs sex set
                $vind->sex('M');
                $vind->alias($sample_name);

                $vind->update;
                $vind_name = $vind->name;

                if( ! $vind->population($pname) )
                {
                    #create a population with the project name
                    $vind->add_population($pname);
                    $vind->population($pname);
                    print qq[Created a new population called $pname\n];
                }

                $vind->update;
                $vind_name = $vind->name;

                if( ! $vind->population($pname) ){
                    print STDERR qq[Failed to get population for new individual!];
                    exit;
                }
            }
        }
        elsif ($samplemap){
            $vind_name = $individual_names{$sample->name};
            unless ($vind_name) {
                print "Can't identify individual from ",$sample->name,". Skipping\n";
                next;
            }
        }
        else {
            my $hname = uc($sample->name);  # hierarchy name - i.e. for filesystem
            $hname =~ s/\W+/_/g;
            foreach (@individual_lookups){
                if ($hname =~ /^$_/){
                    $vind_name = $individual_names{$_};
                    last;
                }
            }
            unless ($vind_name) {
                print "Can't identify individual from ",$sample->name,". Skipping\n";
                next;
            }
        }
        
        my $vsample = get_core_object('sample',$vproject,$sample);
        unless ($vsample){
            print STDERR "Can't get tracking sample for ",$sample->name,"; skipping\n";
            next;
        }

        my $vind = $vsample->individual($vind_name); 
        unless ($vind){
            print STDERR "Can't get individual $vind_name for sample ",$sample->name,". Skipping\n";
            next;
        }

        $vsample->hierarchy_name($vind->hierarchy_name);    # this is so we pool samples under individual
        print "Sample ",$vsample->name," updating\n" if $vsample->dirty;
        $vsample->update;

        # Library Requests
        my $slibrequests;
        $slibrequests = $sample->library_requests();
        foreach my $slibrequest(@$slibrequests){
            my $vlibrequest = get_core_object('library_request',$vsample,$slibrequest);
            unless ($vlibrequest){
                print STDERR "Can't get tracking librequest for ",$slibrequest->id(),"; skipping\n";
                next;
            }
            $vlibrequest->ssid($slibrequest->id());
            $vlibrequest->prep_status($slibrequest->status());
            print "Library Request ",$vlibrequest->id()," updating\n" if $vlibrequest->dirty;
            $vlibrequest->update;

            #Libraries
            my $libs;
            eval {
                $libs = $slibrequest->libraries;
            };
            if ($@){
                print STDERR "Error getting libraries off ",$sample->name," : ",$@,". Skipping\n";
                next;
            }
                    
            foreach my $lib(@$libs){
                print "Library ",$lib->id,"\n";
                my $vlib = get_core_object('library',$vlibrequest,$lib);
                unless ($vlib){
                    print "Can't get tracking lib for ",$lib->name,"; skipping\n";
                    next;
                }
             

                $vlib->ssid($lib->id);
                $vlib->name($lib->name);
                $vlib->prep_status($lib->prep_status);
                my $is_tagged_library = $lib->is_tagged;
                if ($is_tagged_library){
                    $vlib->library_tag_sequence($lib->tag_sequence);
                    $vlib->library_tag_group($lib->tag_group_id);
                    $vlib->library_tag($lib->tag_id);
                }

                my $vmultiplex_pool;
                my $multiplex_pool_asset_ids=$lib->multiplex_pool_asset_ids();
                if($is_tagged_library && (@$multiplex_pool_asset_ids>0)){
                    # there can be more than 1 multiplex pool created from a single library tube
                    foreach my $multiplex_pool_asset_id (@$multiplex_pool_asset_ids){
                        # find/create the multiplex pool, and link it to the library with a library_multiplex_pool 
                        #print STDOUT "\t\t\tmultiplex_library_tube_asset_id=".$multiplex_library_tube_asset_id."\ttag_id=".$ref->{$multiplex_library_tube_asset_id}."\n";
                        $vmultiplex_pool=VRTrack::Multiplex_pool->new_by_ssid($vrtrack,$multiplex_pool_asset_id);
                        unless($vmultiplex_pool){
                            print "New multiplex pool ",$multiplex_pool_asset_id,"\n";
                            $vmultiplex_pool=VRTrack::Multiplex_pool->create($vrtrack,$multiplex_pool_asset_id);
                        }
                        # TODO: set multiplex_pool->name.  Not available from Sfind right now.
                        my $vlib_multiplex_pool=VRTrack::Library_Multiplex_pool->new_by_library_id_multiplex_pool_id($vrtrack,$vlib->id,$vmultiplex_pool->id);
                        unless($vlib_multiplex_pool){
                            print "New library multiplex pool ",$lib->id,",",$multiplex_pool_asset_id,"\n";     
                            $vlib_multiplex_pool=VRTrack::Library_Multiplex_pool->create($vrtrack,$vlib->id,$vmultiplex_pool->id);
                        }
                        print "Library<->Multiplex pool ",$vlib_multiplex_pool->library_id,"<->",$vlib_multiplex_pool->multiplex_pool_id," updating\n" if $vlib_multiplex_pool->dirty;
                        $vlib_multiplex_pool->update;
                    }                           
                    print "Multiplex Pool ",$vmultiplex_pool->id," updating\n" if $vmultiplex_pool->dirty;
                    $vmultiplex_pool->update;
                }

                my $type=$lib->type;
                my $libtype = $vlib->library_type($type);
                unless ($libtype){
                print "New library type $type\n";
                    $vlib->add_library_type($type);
                }
                my $seq_centre = $vlib->seq_centre('SC');
                unless ($seq_centre){
                    $vlib->add_seq_centre('SC');
                }
                my $seq_tech = $vlib->seq_tech('SLX');
                unless ($seq_tech){
                    $vlib->add_seq_tech('SLX');
                }
                my ($f_from,$f_to) = @{$lib->fragment_size};
                if($f_from && $f_to){
                    $vlib->fragment_size_from($f_from);
                    $vlib->fragment_size_to($f_to);
                }
                
                        
                print "Library ",$vlib->name," updating\n" if $vlib->dirty;
                $vlib->update;
                
                
                # Sequencing requests
                my $seq_requests;
                eval {
                    $seq_requests=$lib->seq_requests;
                };
                if ($@){
                    print STDERR "Error getting sequencing requests  ",$lib->id," : ",$@,". Skipping\n";
                    next;
                }
                foreach my $seq_request(@$seq_requests){
                    my $vseq_request;
                    my $is_multiplexed_seq_request=($seq_request->library_id==$lib->id) ? 0 : 1;
                    # multiplexd sequencing
                    if($is_multiplexed_seq_request){
                        $vmultiplex_pool=VRTrack::Multiplex_pool->new_by_ssid($vrtrack,$seq_request->library_id);
                        $vseq_request = get_core_object('seq_request',$vmultiplex_pool,$seq_request); 
                    }
                    # non multiplexed sequencing
                    else{
                        $vseq_request = get_core_object('seq_request',$vlib,$seq_request);      
                    }
                    unless ($vseq_request){
                        print STDERR "Can't get tracking seq_request for ",$seq_request->id,"; skipping\n";
                        next;
                    }
                    
                    $vseq_request->ssid($seq_request->id);
                    $vseq_request->seq_type($seq_request->type);
                    $vseq_request->seq_status($seq_request->status);

                    print "seq_request ",$vseq_request->ssid," updating\n" if $vseq_request->dirty;
                    $vseq_request->update;


                    #Lanes
                    my $lanes;
                    eval {
                        $lanes = $seq_request->lanes;
                    };
                   
                    if ($@){
                        print STDERR "Error getting lanes off ",$seq_request->id," : ",$@,". Skipping\n";
                        next;
                    }
        
                    foreach my $lane (@$lanes){
                        ########################################################
                        # Don't add lanes without resulting files.  These are
                        # either cancelled (and will never have files), or have
                        # not been archived yet, so will get added later
                        ########################################################
        
                        my $tag_id=$lib->tag_id();

                        # Check for BAM first, then fastq
                        my $files;
                        eval {
                            $files = $lane->bam;
                        };
                        if ($@){
                            print "Error getting bam from ",$lane->name," : ",$@,".  Skipping\n";
                            next;
                        }
                        else {
                            # have the bam, but if this is multiplexed, we'll have them all, so filter on the tag
                            # should maybe do this in Wrapper::iRODS?
                            if ($is_multiplexed_seq_request){
                                #print $lane->name, " $tag_id files:\n\t";
                                #print $_->name."\n\t" foreach @$files;

                                @$files = grep($_->name =~ /#$tag_id\.bam$/, @$files);
                            }
                        }


                        unless (@$files){
                            eval {
                                $files = $lane->fastq;
                            };
                            if ($@){
                                print "Error getting fastq from ",$lane->name," : ",$@,".  Skipping\n";
                                next;
                            }
                        }

                        unless(@$files){
                            print $lane->name, " has no files\n";
                            next;
                        }

                        my $vlane = $vlib->get_lane_by_name($lane->name);
                        my $deplexed_lane;
                        
                        if($is_multiplexed_seq_request){
                            $deplexed_lane=$lane->name."#".$tag_id;
                            $vlane = $vlib->get_lane_by_name($deplexed_lane);
                        }
                        if ($vlane){
                            # Don't do file stuff - it should all have been done when
                            # the lane was added, and shouldn't have changed.  Also,
                            # (a) i/o is slow and (b) if the file has changed (split
                            # old _s_ files, or gzipped or something, then updating
                            # will just bugger things up.
                            # Use files flag to force file update anyway.
                            unless ($update_files){
        
                                # Need to update qc status of pending lanes, but old
                                # pending lanes will never update, so need to skip
                                # these.  A date-based skip would be better.
                                if ($lane->run_name > 2500 && $vlane->npg_qc_status eq 'pending'){
                                    $vlane->npg_qc_status($lane->npg_qc);
                                    print "Lane ",$vlane->name," npg qc updating\n" if $vlane->dirty;
                                    $vlane->update();
                                }
                                next;
                            }
                        }
                        else {
        
                            # this lane doesn't exist on this library.  Does it exist
                            # on any?  If so, barf & skip
                            my $checklane = VRTrack::Lane->new_by_name($vrtrack, $lane->name);
                            if ($checklane && !$is_multiplexed_seq_request){
                                my $checklib = VRTrack::Library->new($vrtrack,$checklane->library_id);
                                print $lane->name," should be on ",$vlib->name," but is on ",$checklib->name," Skipping\n";
                                next;
                            }
                                                                    
                            # jws 2011-01-18 removed check for basepairs.  This is 0 for bam as 
                            # there currently is no bamcheck in place in irods
                            #if ($lane->basepairs && !$is_multiplexed_seq_request){
                            if ($is_multiplexed_seq_request){
                                print "New lane ",$deplexed_lane,"\n";
                                #$vlane = $vlib->add_lane($deplexed_lane);
                                $vlane = $vseq_request->add_lane($deplexed_lane);
				$vlane->hierarchy_name($deplexed_lane);
                                $vlane->library_id($vlib->id);
                                $vlane->update;
                            }
                            else {
                                print "New lane ",$lane->name,"\n";
                                #$vlane = $vlib->add_lane($lane->name);
                                $vlane = $vseq_request->add_lane($lane->name);
                            }
                        }
                        $vlane->npg_qc_status($lane->npg_qc);
                        $vlane->read_len($lane->read_len);
                        $vlane->run_date($lane->created);
                        if ($lane->is_paired){
                            $vlane->is_paired(1);
                        }
                        else {
                            $vlane->is_paired(0);
                        }
                        if ($lane->basepairs < 0){
                            warn "problem with lane $lane fastqcheck - probably overflow\n";
                        }
                        $vlane->raw_reads($lane->reads);
                        $vlane->raw_bases($lane->basepairs);
                        print "Lane ",$vlane->name," updating\n" if $vlane->dirty;
                        $vlane->update;
                        
                        foreach my $file(@$files){
                            # OK, so file names in tracking get changed as the file
                            # is imported, possibly split (for old _s_ files), and 
                            # compressed.  Need to check for the likely names this file
                            # may now have.
                            # check in reverse order of processing to speed things up -
                            # i.e. it's most likely that a file is present in the 
                            # compressed form, so check for that first
                            my $filename = $file->name;
                            my $vfile;
                            my @check_names;
        
                            # handle _s_ files.  These are where the fwd & rev reads
                            # are concatenated on a single line, and are processed by
                            # splitting into _1.fastq and _2.fastq files.
                            # We'll assume the file is in the db if we hit either fwd
                            # or reverse split files, or the _s_ file.

                            my $deplex_filename;    # need this later if we have to add the file
                            if($is_multiplexed_seq_request){
                                # handle s files
                                if ($filename =~ /(\d+)_s_(\d)\.fastq$/){
                                    die "Shouldn't have multiplexed _s_ files: $filename on ".$vlane->name."\n";
                                }

                                if ($filename =~ /\.bam$/){
                                    # should already be deplexed, so no problem
                                    $deplex_filename = $filename;
                                }
                                else {  # presumably fastq
                                    $deplex_filename = deplex_name_from_fastq_tag($filename, $tag_id);
                                    next unless $deplex_filename;
                                }
                                push @check_names, $deplex_filename;
                            }
                            # non multiplexed requests
                            else{
                                # handle s files
                                if ($filename =~ /(\d+)_s_(\d).fastq/){
                                    my $splitfwd;
                                    my $splitrev;
                                    $splitfwd = "$1_$2_1.fastq";
                                    $splitrev = "$1_$2_2.fastq";
                                    @check_names = ($splitfwd,$splitrev);
                                }
                                # handle other files
                                else{
                                   push @check_names, $filename;   
                                }   
                            }

                            # Right, now check if we already have either .gz or not .gz files
                            foreach my $check_name (@check_names){
                                last if $vfile;
                                $vfile = $vlane->get_file_by_name($check_name);
                                unless ($vfile){
                                    $vfile = $vlane->get_file_by_name("$check_name.gz");
                                }
                            }

                            if ($vfile){
                                # we already have the file in the vrtrack database.
                                # Don't update unless forced by --update_files.  It's slow, and shouldn't change.
                                unless ($update_files){
                                    next;
                                }
                            }
                            else {
                                # We have a new file that isn't in the vrtrack db, so add new one
                                if($is_multiplexed_seq_request){
                                    print "New file: ",$deplex_filename,"\n"; 
                                    $vfile = $vlane->add_file($deplex_filename);
                                    $vfile->hierarchy_name($deplex_filename);
                                }
                                else{
                                    print "New file: ",$filename,"\n"; 
                                    $vfile = $vlane->add_file($filename);
                                }
                            }
                            $vfile->raw_reads($file->reads);
                            $vfile->raw_bases($file->basepairs);
                            $vfile->read_len($file->read_len);
                            $vfile->mean_q($file->mean_quality);
                            $vfile->md5($file->md5);
                            print "File ",$vfile->name," updating\n" if $vfile->dirty;
                            # determine type
                            #   0 is single-end
                            #   1 fwd
                            #   2 rev
                            #   3 is for _s_ files
                            #   4 is for bam files
                            if ($filename =~ /^\d+_s_\d.fastq/){
                                $vfile->type(3);
                            }
                            elsif ($filename =~ /^\d+_\d_1.fastq/){
                                $vfile->type(1);
                            }
                            elsif ($filename =~ /^\d+_\d_2.fastq/){
                                $vfile->type(2);
                            }
                            elsif ($filename =~ /^\d+_\d.fastq/){
                                $vfile->type(0);
                            }
                            else {
                                print "Can't determine type of file $filename\n";
                            }
        
                            $vfile->update;
                        } # Fastqs
                    } #Lanes
                } # Seq Requests       
            } # Library
        } # Library Requests
    } # Sample 
} # Study
print "Done\n";

# sub to get a core object, or create a new one if required
# checks that the object has the right ssid and name and is not a duplicate
# e.g. my $proj = get_core_object('project',$vrtrack,$sproj);
sub get_core_object {
    my ($type, $parent, $ss_obj) = @_;
    #print "Called with $type ",$ss_obj->name,"\n";
    $type = lc($type);
    my $class = ucfirst($type);
    my $get_ssid_method="get_${type}_by_ssid";       
    my $add_method="add_${type}";       

    # Right, try and get an existing object with the same sequencescape id
    my $obj = $parent->$get_ssid_method($ss_obj->id);
    #print "$parent $get_ssid_method, $ss_obj:", $ss_obj->id,"\n";
    if(!($type eq 'library_request' || $type eq 'seq_request' || $type eq 'multiplex_seq_request')){
        if ($obj){
            if ($obj->name ne $ss_obj->name){
                print "$class ",$ss_obj->id," name has changed (",$obj->name," to ",$ss_obj->name,")\n";
                $obj->name($ss_obj->name);
            }
        }
        else {
            # No obj with this ssid belongs to this parent, but let's see if we
            # have an obj with this ssid that we've swapped due to genotype changes
            $obj = "VRTrack::$class"->new_by_ssid($parent->vrtrack,$ss_obj->id);
            if ($obj){
                # TODO: if lane, check genotype for swap.  What to do if library
                # or sample?
                print "$class ",$obj->name," exists (by ssid) but does not belong to sequencescape parent ",$parent->id,"\n";
            }
        }
    }
    unless ($obj) {
        # no obj by ssid, so make new one by name, but first check if it
        # already exists by name
        my $nameobj;
        my $name;
        if (!($class eq 'Library_request' || $class eq 'Seq_request' )) {
            $name=$ss_obj->name;
            if ($class eq 'Sample'){    
                $nameobj = VRTrack::Sample->new_by_name_project($parent->vrtrack,$name,$parent->id);
            }
            else {
                $nameobj = "VRTrack::$class"->new_by_name($parent->vrtrack,$name);
            }
        }
        if ($nameobj){
            print "$class name $name is already in the database, but has the wrong ssid\n";
            $obj=undef;
        }
        else {
            print "New $type $name\n";
            if($type eq 'library_request' || $type eq 'seq_request' ){
                $obj = $parent->$add_method($ss_obj->id);
            }
            else{
                $obj = $parent->$add_method($name);
            }
            $obj->ssid($ss_obj->id);
        }
    }

    return $obj;
}


# sub to work out what the deplexed filename should be for a given fastq file
# and library tag
sub deplex_name_from_fastq_tag {
    my ($filename, $tag_id) = @_;
    my $deplex_name;

    if($filename =~ /^(\d+_\d)_(1.fastq)$/){
        $deplex_name = $1.'#'.$tag_id.'_'.$2;
    }
    elsif ($filename =~ /^(\d+_\d)_(2.fastq)$/){
        $deplex_name = $1.'#'.$tag_id.'_'.$2;
    }
    elsif ($filename =~ /^(\d+)_(\d.fastq)$/){
        $deplex_name = $1.'#'.$tag_id.'_'.$2;
    }
    else {
        print "Can't parse multiplex fastq name to generate deplexed filename from $filename\n";
    }
    return $deplex_name;
}
