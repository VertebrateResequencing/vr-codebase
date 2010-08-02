=head1 NAME

=head1 VertRes::File

# create: at minimum we need to know the path and type
my $file_obj = VertRes::File->new(path => '/abs/path/filename.bam',
                                  type => 'bam');

# new() calls:
$file_obj->refresh_stats();
# above stores stats about file in db, including the time the stats were
# obtained. Even if the file doesn't exist yet we still store what we can.
# when stats are subsequently requested, by default they come from db without a
# disc-check.

# reference:
my $file_id = $file_obj->id;
# gives us a unique id that other systems and dbs can refer to or store. The
# same id is always returned for every file_obj created with the same 'path'
# option.

# stringify:
print "$file_obj\n"; # /abs/path/filename.bam

# before getting stats about a file, we can control if we'll use the existing
# db-stats or update from a new disc-check
my $updated = $file_obj->refresh_stats();
# always checks on disc and stores answer in db, return true if the file had
# been modified since the previous refresh (or if this is the first time the
# file was detected as existing). A refresh does not write to the db unless
# the file had been modified, and so things like md5 and num_records are not
# touched; otherwise they are cleared and will have to be calculated again
$file_obj->non_existant_behaviour('always_refresh');
# by default, if a file is non-existant in the db, a refresh_stats() is
# triggered prior to any get stat method to see if it exists now. This can be
# turned off above temporarily (while $file_obj is in scope) with
# 'never_refresh'.
$file_obj->expiry_time(86400); VertRes::File::expiry_time(86400);
# stats that were taken more than this many seconds ago trigger a
# refresh_stats(), can work across all file instances with the class method.

# get stats about the file from db
if ($file_obj->e) { # yay! the file exists! }
my $size_in_bytes = $file_obj->s;
my $modified = $file_obj->last_modify_time; # in seconds since the epoch
my $last_refresh = $file_obj->last_disc_check; # "
my $md5sum = $file_obj->md5; # only calculated when first called
my $num_records = $file_obj->num_records; # the number of records in the file,
# only calculated when first called. it will store record counts in db so only
# has to parse any given file once:
unless ($file_obj->num_records) {
    my $type = $file_obj->type;
    require "VertRes::FileTypeHandler::$type";
    my $handler = VertRes::FileTypeHandler::$type->new;
    $file_obj->refresh_stats(); # so that we'll have the correct last modified time for the answer we'll get
    my $num_records = $handler->num_records($file_obj); # parse file and get answer
    $file_obj->num_records($num_records); # answer stored in db
}

# truncation checking
if ($file_obj->complete(expected_records => 50)) { ... }
if ($file_obj->complete(expected_records_from => [$file_obj2, $file_obj3])) { ... }
# compares summed total num_records in the given files to the num in this file



=head1 VertRes::FileTypeHandler::*




=cut



=head2 split_requires

 Title   : split_requires
 Usage   : my $required_files = $obj->split_requires('/path/to/lane');
 Function: Find out what files the split action needs before it will run.
 Returns : array ref of file names
 Args    : lane path

=cut


1;
