=head1 NAME

=head1 VertRes::Persistent

 We want to be able to store/retrieve persistent values in a database, but
 don't want to be specific to a particular driver. We don't even want to have
 to use SQL: we could store our information in a file or memory if need be.
 We'll implement using the CPAN Persistent framework. We'll wrap around it and
 use it 'indirectly' so that we can avoid having to specifiy which Persistent
 module ('driver') we're using, leaving that up to an environment variable,
 class or instance method setting to decide at run-time.

 This is only used to implement certain VertRes:: modules, not for direct use by
 end-users.

 # setup:
 use base VertRes::Persistent;
 sub initialize {
  my $self = shift;
  my $driver = $self->driver; # eg. MySQL, chosen by an environment variable
  VertRes::Persistent::driver('File'); # override globablly
  $self->driver('File'); # override locally
  $self->datastore(...); # override driver-specific settings, like db name

  $self->SUPER::initialize(@_);
  # At this point, $self stores in itself a correctly setup Persistent::$driver.

  $self->add_attribute(...);
  # passes through to Persistent, doing the eval dance internally within
  # VertRes::Persistent.
 }

 # get/set:
 $self->$attribute(); # if not overriden, passes through to VertRes::Persistent
                      # AUTOLOAD, which passes to:
 $self->value($attribute, $value); # which passes to the Persistent object

 # Does an auto-save on every value set. Does an auto-restore during object
 # creation. Does not bother wrapping any of the other Persistent methods.
 # If you really need them, provides access to the the Persistent object with
 # $self->persistent(). auto-restore works if the necessary keys have been
 # supplied to new(), or if new(id => $id) was called. A unique (auto-increment)
 # id is associated with every sub-class key set.

=head1 VertRes::FileTypeHandler*

 These classes provide 'best' ways of doing common operations on files of a
 certain filetype. We have VertRes::Parser::*, but not all filetypes necessarily
 need a parser, or there is a fast way to do a particular operation which does
 not use normal ParserI methods. As with Parser modules, being compressed should
 be transparent; that is, .gz is not a filetype.

 # create a handler without knowing the filetype; it will guess, defaulting to
 # text type:
 my $handler = VertRes::FileTypeHandler->new(file => "myfile");

 # or if you know the filetype:
 $handler = VertRes::FileTypeHandler->new(file => "myfile", type => "bam");
 # == VertRes::FileTypeHandler::bam->new(file => "myfile");
 # VertRes::FileTypeHandler::<filetype> modules implement
 # VertRes::FileTypeHandlerI and we become them automatically.
 # can also accept a VertRes::File instead of string for file option, and get
 # type from that; VertRes::IO would get the path from it

 # do some common operations:
 my $records = $handler->num_records; # get the number of records in the file
 my $header_lines = $handler->num_header;
 unless ($handler->is_indexed) {
    $handler->index; # eg. tabix for tab-seperated files
 }

=head1 VertRes::File

 This class describes a file-on-disc, storing information about it in our
 database. The idea is to minimise accessing the filesystem as much as possible.
 It is implemented using VertRes::Persistent.

 # create: at minimum we need to know the path and type
 my $file_obj = VertRes::File->new(path => '/abs/path/filename.bam',
                                   type => 'bam');
 # or by id:
 $file_obj = VertRes::File->new(id => $a_valid_file_id);

 # new() calls:
 $file_obj->refresh_stats();
 # above stores stats about file in db, including the time the stats were
 # obtained. Even if the file doesn't exist yet we still store what we can.
 # when stats are subsequently requested, by default they come from db without a
 # disc-check. This method clears out any existing data in the db for this file
 # and replaces it with the results of a check-on-disc if the last modify time
 # on disc is greater than that in the db (or if the file does not exist).

 # get a reference to the file entry in the db:
 my $file_id = $file_obj->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # same id is always returned for every file_obj created with the same 'path'
 # option.

 # stringify:
 print "$file_obj\n"; # /abs/path/filename.bam

 # open:
 my $fh = $file_obj->open;
 # it attempts an open if the db says the file exists, and if it fails due to
 # the file not existing it will refresh_stats() before throwing an error. It
 # will also throw if called when the db says the file does not exist.

 # unlink:
 $file_obj->unlink;
 # deletes a file and then does refresh_stats().


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
     my $handler = VertRes::FileTypeHandler->new(file => $file_obj);
     # so that we'll have the correct last modified time for the answer:
     $file_obj->refresh_stats();
     # parse file and get answer:
     my $num_records = $handler->num_records();
     $file_obj->num_records($num_records); # answer stored in db
 }

 # truncation checking
 if ($file_obj->complete(expected_records => 50)) { ... }
 if ($file_obj->complete(expected_records_from => [$file_obj2, $file_obj3])) { }
 # compares summed total num_records in the given files to the num in this file

=head1 VertRes::JobManager::Job

 This class describes information about a command we want to run via the shell.
 The information is stored in a db (implemented using VertRes::Persistent).

 # create, with at least the command line specified, also requiring the working
 # directory, but this defaults to the current working directory:
 my $job = VertRes::JobManager::Job->new(cmd => "myshellcommand -a foo bar",
                                         dir => '/abs/path/');
 # or by id:
 $job = VertRes::JobManager::Job->new(id => $a_valid_job_id);

 # get a reference to the job entry in the db:
 my $job_id = $job->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # same id is always returned for every job object created with the same 'cmd'
 # and 'dir' pair.

 # find out about the job:
 my $cmd = $job->cmd; # $job also stringifies to the cmd
 my $dir = $job->dir;

 if ($job->pending) { # no attempt has been made to run the cmd yet }
 elsif ($job->running) {
    # the cmd is running right now
    my $start_time = $job->start_time; # in seconds since the epoch
    my $pid = $job->pid; #
    my $host = $job->host; # the machine we're running on
    my $beat = $job->last_heartbeat; # elapsed time since the object sent a
                                     # heartbeat to signify it was still OK
 }
 elsif ($job->finished) {
    # the cmd had been run in the past, and is now no longer running
    my $end_time = $job->end_time;
    my $elapsed_time = $job->wall_time;
    if ($job->ok) { # the cmd returned a good exit value }
    else {
        # the cmd failed
        my $exit_code = $job->exit_code;
        my $error_message = $job->err; # not STDERR of the cmd, but the error
                                       # message generated at the time of cmd
                                       # failure
        my $retries = $job->retries; # the number of times this cmd has reached
                                     # finished() state
    }
 }

 # run the command:
 $job->run;
 # This sets the running state in db, then chdir to ->dir, then forks to run
 # run the ->cmd (updating ->pid, ->host and ->start_time), then when finished
 # updates ->running, ->finished and success/error-related methods as
 # appropriate.
 # If ->run() is called when the job is already in running() state, it will by
 # default wait until the job has finished (polling db every minute) and then
 # run it again itself. Change that behaviour:
 $job->run(block_and_skip_if_ok => 1);
 # this waits until the job has finished running, runs it again if it failed,
 # otherwise does nothing so that $job->ok will be true. This is useful if you
 # start a bunch of tasks that all need to first index a reference file before
 # doing something on their own input files: the reference index job would be
 # run with block_and_skip_if_ok.
 # While running we update heartbeat() in the db every 30mins so that other
 # processes can query and see if we're still alive.

 # kill a job, perhaps one you suspect has got 'stuck' and won't die normally:
 if ($job->running && $job->last_heartbeat > 3600) {
    $job->kill;
    # this logs into the ->host if necessary and does a kill -9 on the ->pid.
    # if this doesn't result in the $job becoming ->finished after giving it
    # some time to react, then we forcibly set finished and custom err and
    # exit_code ourselves.
 }

=head1 VertRes::JobManager::Requirements

 This class describes the expected requirements of a job. The information is
 stored in a db (implemented using VertRes::Persistent).

 # create, specifying at least memory and time (others have defaults as shown,
 # all are used for the 'key'):
 my $req = VertRes::JobManager::Requirements->new(
                memory  => 3800, # maximum expected memory in MB
                time => 8, # maximum expected hrs
                cpus => 1, # the number of cpus we should reserved
                tmp_space => 0, # tmp space that should be reserved in MB
                local_space => 0, # as above, but for special local non-tmp
                custom => '' # custom stuff that might be understood by a
                             # certain scheduler
           );
 # or by id:
 $req = VertRes::JobManager::Requirements->new(id => $a_valid_requirements_id);

 # get a reference to the requirements entry in the db:
 my $req_id = $req->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # same id is always returned for every Requirements object created with the
 # same full specification of requirements.

 # get the requirements:
 my $memory = $req->memory;
 # and so on for the other things passed to new();

 # make a new requirement based on an existing one
 my $new_req = $req->clone(time => 9);

=head1 VertRes::JobManager::Submission

 This class describes information about a job submission, which is a request to
 run a particular job within a particular job scheduler with certain hints to
 that scheduler about what the job requries.
 The information is stored in a db (implemented using VertRes::Persistent).

 # create, with at least the job_id and requirements_id specified.
 my $sub = VertRes::JobManager::Submission->new(
                job_id => $vertres_jobmanager_job->id,
                requirements_id => $vertres_jobmanager_requirements->id
           );
 # or by id:
 $sub = VertRes::JobManager::Submission->new(id => $a_valid_submission_id);

 # get a reference to the submission entry in the db:
 my $sub_id = $sub->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # same id is always returned for every Submission object created with the
 # same job_id and requirements_id. Note that the scheduler attribute is stored
 # in the db, but does not form part of the key, since it can be dynamically
 # overridden (see below) or actually changed without affecting the essential
 # properties and outcome of a submission.

 # override scheduler:
 ...->new(..., scheduler => 'Local');
 # above sets the scheduler for this submission only
 VertRes::JobManager::Submission::scheduler('Local');
 # above changes the scheduler stored in the db for newly created Submission
 # objects when scheduler => isn't supplied to new(), and it also changes the
 # return value of $sub->scheduler on existing objects, without changing the
 #Êvalue in the db.
 # For possible schedulers, see VertRes::JobManager::Scheduler*

 # find out about and manipulate the submission:
 my $job_id = $sub->job_id;
 my $requirements_id = $sub->requirements_id;
 my $scheduler = $sub->scheduler; # eg. the string 'LSF';

 if ($sub->scheduled) {
    # submission has been scheduled by the scheduler
    my $sid = $sub->sid; # the id of the job given to us by the scheduler
    if ($sub->running) {
        # our associated job is running (according to the job, not the
        # scheduler)
    }
    elsif ($sub->finished) {
        # our associated job had been run in the past, and is now no longer
        # running (according to the job, not the scheduler)
    }
 }
 else { # no action has been taken with this submission yet }
 
 # *** above not quite right...

=head1 VertRes::JobManager::Scheduler*

 This class 

=cut

1;
