=head1 Overview

 This document describes all the modules and methods necessary to build a new
 database-driven pipeline system. We go from the fundamental building blocks
 at the start up to the modules and methods you're more likely to use yourself
 on a daily-basis by the end of this document. Later modules assume knowledge of
 earlier ones in the POD (and most likely use or inherit from them in the code),
 so it is recommended to read this document in the order presented.

 On overview diagram is available here: ...

=head1 VertRes::SiteConfig & VertRes::Utils::Config & Build.pl

 We need set-once, site-specific configuration. Probably implement something
 along the lines of http://www.perlmonks.org/?node_id=464358 with config stored
 as perl data structure in modules/VertRes/SiteConfig.pm, VertRes::Utils::Config
 providing get/set methods and the Build.pl script taking us through an
 interactive process to set appropriate values. Once all set
 modules/VertRes/SiteConfig.pm will be created and will be installable in the
 normal way.

 We need to know what the available job schedulers are, which one should be used
 by default, and which queues should be used in order of preference. Though it's
 possible to query LSF what the limits on a queue are, we may as well also just
 ask the user to supply max memory, time and discspace per queue. We'll also ask
 for the default root directory in which we'll store STDOUT and STDERR files.
 And we'll ask for the global maximum number of jobs we can have scheduled per-
 user across all queues (default 0 for infinity).

 We'll also store MySQL access details in the config file. One of the details
 will be a database name prefix that any VertRes:: module wanting to create a
 database will use. We'll do something special to store the password encrypted
 (with the key stored in a user-protectable file outside of the repository they
 passed in as an earlier config option) and we might also adjust our version
 control settings so that modules/VertRes/SiteConfig.pm is ignored for adding
 and commits.

 Details stored in the SiteConfig.pm file will be easily accessible to other
 VertRes::* modules via VertRes::Utils::Config when they need get a default
 value for something.

 # basic set/getting:
 my $config = VertRes::Utils::Config->new;
 $config->set('databases', 'mysql', 'username', $value);
 # above sets $VertRes::SiteConfig::CONFIG{databases}->{mysql}->{username}
 my $username = $config->get('databases', 'mysql', 'username');

 # encryption of passwords:
 $config->set_encrypted(['databases', 'mysql', 'password', $raw_value],
                        $config->get('key_file'));

 # set/get a list/hash
 $config->set('schedulers', [@values]); # or {%values}
 my @known_schedulers = $config->get('schedulers');

=head1 VertRes::Persistent

 We want to be able to store/retrieve persistent values in a database, but
 don't want to be specific to a particular driver. We don't even want to have
 to use SQL: we could store our information in a file or memory if need be.
 We'll implement using something like the CPAN Persistent framework. We'll wrap
 around it and use it 'indirectly' so that we can avoid having to specifiy which
 Persistent module ('driver') we're using, leaving that up to SiteConfig
 variable, class or instance method setting to decide at run-time.

 This is only used to implement certain VertRes:: modules, not for direct use by
 end-users.

 # setup:
 use base VertRes::Persistent;
 sub initialize {
  my $self = shift;
  my $driver = $self->driver; # eg. MySQL, chosen by a SiteConfig variable
  VertRes::Persistent::driver('File'); # override globablly
  $self->driver('File'); # override locally
  $self->datastore(...); # override driver-specific settings, like db name

  $self->SUPER::initialize(@_);
  # At this point, $self stores in itself a correctly setup Persistent::$driver.

  $self->add_attribute(...);
  # passes through to Persistent, doing the eval dance internally within
  # VertRes::Persistent. Regardless of the underlying system, there is some
  # mechanism for saying which attributes form the unique key
 }

 # get/set:
 $self->$attribute(); # if not overriden, passes through to VertRes::Persistent
                      # AUTOLOAD, which passes to:
 $self->value($attribute, $value); # which passes to the Persistent object
 # unique key attributes cannot be set

 # we can get/set arbitrarily complex data structures (mixtures of arrays and
 # hashes and scalars only... perhaps even code refs?).
 # most likely implemented by simple string concatenation with special
 # seperators, and automatically turned back into the right structure (if a list
 # goes in to a set, a list comes out during the get).

 # Does an auto-save on every value set. Does an auto-restore during object
 # creation. Does not bother wrapping any of the other Persistent methods.
 # If you really need them, provides access to the the Persistent object with
 # $self->persistent(). Auto-restore works if the necessary keys have been
 # supplied to new(), or if new(id => $id) was called. A unique (auto-increment)
 # id is associated with every sub-class key set.

 # For consistency and speed, objects are cached, so in
 # [$a = $class->new(id => 5); $b = $class->new(id => 5);], $a and $b point to
 # the same location in memory. To protect against multiple processes accessing
 # the same id, gets are never cached. You can store a list of all the
 # Persistent objects created like this:
 my @object_array;
 VertRes::Persistent::capture_start(\@object_array);
 # ... use classes that inherit from VertRes::Persistent ...
 VertRes::Persistent::capture_stop;
 # @object_array now contains every Persistent-based object that was created
 # between start and stop.

 # All VertRes::Persistent-based classes will also have a last_updated()
 # auto-set during a save() which will be used to delete old 'rows' from the db
 # when expunge_old($num_of_days) is called. expunge_old() calls the class-only-
 # method expunge($id), which children can overwrite. A touch() method updates
 # last_updated() without having to set anything else (so stuff you're reading
 # from but aren't writing to don't get expunged by expunge_old).

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
 database. The idea is to minimise accessing the filesystem as much as possible;
 VertRes::File would be used pervasively every single time any kind of file
 access is done in all VertRes::* code. It is implemented using
 VertRes::Persistent.

 # create: at minimum we need to know the path and type
 my $file_obj = VertRes::File->new(path => '/abs/path/filename.bam',
                                   type => 'bam');
 # or by id:
 $file_obj = VertRes::File->new(id => $a_valid_file_id);

 # new() calls:
 $file_obj->refresh_stats();
 # Above stores stats about file in db, including the time the stats were
 # obtained. Even if the file doesn't exist yet we still store what we can.
 # When stats are subsequently requested, by default they come from db without a
 # disc-check. This method clears out any existing data in the db for this file
 # and replaces it with the results of a check-on-disc if the last modify time
 # on disc is greater than that in the db (or if the file does not exist).

 # get a reference to the file entry in the db:
 my $file_id = $file_obj->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # same id is always returned for every file_obj created with the same 'path'
 # option.

 # open:
 my $fh = $file_obj->open;
 # it attempts an open if the db says the file exists, and if it fails due to
 # the file not existing it will refresh_stats() before throwing an error. It
 # will also throw if called when the db says the file does not exist.

 # unlink:
 $file_obj->unlink;
 # deletes a file and then does refresh_stats().

 # move:
 $file_obj->mv('/abs/path/new/location');
 # does a safe copy to the location, checks the md5s match, creates a new
 # file_obj for the new file, then does $self->unlink


 # before getting stats about a file, we can control if we'll use the existing
 # db-stats or update from a new disc-check

 my $updated = $file_obj->refresh_stats();
 # always checks on disc and stores answer in db, return true if the file had
 # been modified since the previous refresh (or if the file now exists when it
 # did not at the previous refresh). A refresh does not write to the db unless
 # the file had been modified, and so things like md5 and num_records are not
 # touched; otherwise they are cleared and will have to be calculated again

 $file_obj->non_existant_behaviour('always_refresh');
 # by default, if a file is non-existant in the db, a refresh_stats() is
 # triggered prior to any get stat method to see if it exists now. This can be
 # turned off above temporarily (while $file_obj is in scope) with
 # 'never_refresh', or with the class method to apply globally.

 $file_obj->expiry_time(86400); VertRes::File::expiry_time(86400);
 # stats that were taken more than this many seconds ago trigger a
 # refresh_stats(), can work across all file instances with the class method.

 # get stats about the file from db
 if ($file_obj->e) { # yay! the file exists! }
 my $size_in_bytes = $file_obj->s;
 my $modified = $file_obj->last_modify_time; # in seconds since the epoch
 my $last_refresh = $file_obj->last_disc_check; # "

 my $md5sum = $file_obj->md5; # only calculated when first called in 'get' mode

 my $num_records = $file_obj->num_records; # the number of records in the file,
 # only calculated when first called in 'get' mode. It will store record counts
 # in db so only has to parse any given file once:
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
 my $cmd = $job->cmd;
 my $dir = $job->dir;

 if ($job->pending) {
    # no attempt has been made to run the cmd yet, or we were reset
 }
 elsif ($job->running) {
    # the cmd is running right now
    my $start_time = $job->start_time; # in seconds since the epoch
    my $pid = $job->pid; # process id
    my $host = $job->host; # the machine we're running on
    my $user = $job->user; # the user the process belongs to
    my $beat = $job->last_heartbeat; # elapsed time since the object sent a
                                     # heartbeat to signify it was still OK
 }
 elsif ($job->finished) {
    # the cmd had been run in the past, and is now no longer running
    my $end_time = $job->end_time;
    my $elapsed_time = $job->wall_time;
    # if called while ->running, it's current time() - ->start_time, otherwise
    # it's ->end_time - ->start_time.

    if ($job->ok) { # the cmd returned a good exit value }
    else {
        # the cmd failed
        my $exit_code = $job->exit_code;
        my $error_message = $job->err; # not STDERR of the cmd, but the error
                                       # message generated at the time of cmd
                                       # failure
        my $retries = $job->retries; # the number of times this cmd has been
                                     # ->reset()
    }
 }

 # run the command:
 $job->run;
 # This sets the running state in db, then chdir to ->dir, then forks to run
 # run the ->cmd (updating ->pid, ->host, ->user and ->start_time), then when
 # finished updates ->running, ->finished and success/error-related methods as
 # appropriate. If ->run() is called when the job is not in ->pending() state,
 # it will by default return if ->finished && ->ok, otherwise will throw an
 # error. Change that behaviour:
 $job->run(block_and_repeat => 1);
 # this waits until the job has finished (polling db every minute) and then
 # calls ->reset and runs it again.
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
    # some time to react, then if we are the matching ->user we forcibly set
    # finished and custom err and exit_code ourselves and return true. If we're
    # not the user that started the job we still try our best to kill (we might
    # be root), and on failure return false without setting finished etc. We
    # hope that the owner will come along and do $job->kill instead.
 }
 
 # reset state:
 $job->reset;
 # this clears the various success/fail states and puts us back into pending.
 # The count returned by get-only method ->retries() is incremented. It will not
 # function if the job is ->running. If you want to reset a running job, you
 # have to ->kill it first.

=head1 VertRes::JobManager::Requirements

 This class describes the expected requirements of a job. The information is
 stored in a db (implemented using VertRes::Persistent).

 # create, specifying at least memory and time (others have defaults as shown,
 # all are used for the 'key'):
 my $req = VertRes::JobManager::Requirements->new(
                memory  => 3800, # maximum expected memory in MB
                time => 8, # maximum expected hrs
                cpus => 1, # the number of cpus we should reserve
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
 # and so on for the other things passed to new(). These cannot be set! They're
 # supposed to be immutable.

 # make a new Requirements based on another with one or more attributes changed:
 my $new_req = $req->clone(time => 10);

=head1 VertRes::JobManager::Submission

 This class describes information about a job submission, which is a request to
 run a particular job within a particular job scheduler with certain hints to
 that scheduler about what the job requries.
 The information is stored in a db (implemented using VertRes::Persistent).

 # create, with at least a Job and Requirements specified. A Scheduler is also
 # needed, but the default Scheduler (from SiteConfig) will be set if not
 # specified. Internally it just stores the object ids, and only the Job->id
 # forms the key.
 my $sub = VertRes::JobManager::Submission->new(
                job => $vertres_jobmanager_job,
                requirements => $vertres_jobmanager_requirements
           );
 # or by id:
 $sub = VertRes::JobManager::Submission->new(id => $a_valid_submission_id);

 # get a reference to the submission entry in the db:
 my $sub_id = $sub->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # same id is always returned for every Submission object created with the
 # same Job->id. Note that the scheduler and requirements attributes are
 # stored in the db, but do not form part of the key, since scheduler can be
 # dynamically overridden (see below) or actually changed without affecting the
 # essential properties and outcome of a submission, and if requirements change
 # we want to still get back the same submission object (which will also prevent
 # multiple processes trying to run the same job with different requiremnts at
 # once).

 # override scheduler:
 ...->new(..., scheduler => $vertres_jobmanager_scheduler);
 # above sets the scheduler for this submission only
 VertRes::JobManager::Submission::scheduler($vertres_jobmanager_scheduler);
 # above changes the scheduler stored in the db for newly created Submission
 # objects when scheduler => isn't supplied to new(), and it also changes the
 # return value of $sub->scheduler on existing objects, without changing the
 #Êvalue in the db.

 # find out about and manipulate the submission:
 my $job = $sub->job;
 my $req = $sub->requirements;
 my $sch = $sub->scheduler;
 my $retries = $sub->retries;

 # it provides get-only methods to all the requirements of our requirements
 # object, which just pass through:
 my $memory = $sub->memory;
 # however we also have extra_* methods for time, memory, cpus, tmp_space and
 # local_space which ->requirements->clone($time => $current_time + $extra) each
 # time they're called, then associate us with the new Requirements object:
 my $time = $sub->time; # eg. 4, == $sub->requirements->time
 $sub->extra_time(2);
 $sub->extra_time(3);
 $time = $sub->time; # now this returns 9
 # This is done so that a particular submission can have its time increased
 # on-the-fly (but we remain the same submission object with the same id),
 # whilst a scheduler might then see the new time should put the submission in
 # a different queue, and potentially the scheduler could then be instructed to
 # switch the queue for this submission. You might call extra_memory() for a
 # submission that failed due to running out of memory, so that a second attempt
 # running the submission might work.

 if ($sub->failed) {
    # we're currently scheduled, we had been run to ->job->finished status,
    # but the job died.

    # If we failed for some known reason like running out of
    # memory, we might change our requirements and retry:
    if ($sub->scheduler->ran_out_of_memory(fh => $sub->last_stdout) {
        # double the required memory
        $sub->extra_memory($sub->memory);
        $sub->reset;
    }

    # of if we failed for some unknown reason we might just try it again for
    # luck, 3rd time's the charm:
    if ($sub->retries < 3) { $sub->reset }
    
    # reset() increments the count reported by ->retries(). It only works when
    # ->failed or ->done. It clears scheduled() (and double-checks that sid() is
    # clear), and does ->job->reset. Once reset, we're now free to retry the
    # submission as normal.
 }
 elsif ($sub->done) {
    # nothing to do, we're done!
 }
 elsif ($sub->scheduled) {
    # when a user uses us to schedule a job, they set sid() to the id they got
    # back from the scheduler
    my $sid = $sub->sid; 
    # and doing so auto-sets scheduled() to the time that sid() was set.
    # when we are scheduled, ->scheduler() can no longer be set, and neither can
    # its return value be overriden with the class method. scheduled() cannot
    # be set except by special internal method for use by sid() and reset()
    # only.
    
    if ($sub->job->running) {
        # user's scheduler might kill the submission if it runs too long in the
        # queue it was initially submitted to; user could do something like this
        # every 15mins:
        if ($sub->close_to_time_limit(30)) {
            # where close_to_time_limit returns true if $sub->job->wall_time >
            # $job->requirements->time - 30.
            $sub->extra_time(2);
            # now make your scheduler recalculate the appropriate queue of a job
            # that takes 2 extra hours, and if the queue changes, switch the
            # queue:
            $sub->scheduler->switch_queues($sub);
        }
    }
    elsif ($sub->job->finished) {
        $sub->update_status();
        # this checks ->job->ok and sets the status to 'done' if ok, otherwise
        # 'failed'. The value of the status we just set is retrieved via the
        # get-only boolean methods ->done and ->failed. It only does something
        # when ->job->finished. It also clears ->sid() using a special
        # internal-only method. scheduled() is still set though, so $sub remains
        # unclaimable and so won't get resubmitted. Prior to clearing sid() it
        # also auto-calls the following 2 methods:
        
        $sub->sync_scheduler();
        # this will wait a reasonable amount of time for the scheduler to report
        # that the sid is finished. After that amount of time it will get the
        # scheduler to kill the sid, and failing that, force remove it from its
        # list of jobs. We don't care about the return status from the
        # scheduler (we trust the ->job), but we do this to make sure our
        # STDERR and STDOUT files have been finished writing to and are closed.
        
        $sub->archive_output();
        # this method asks the scheduler where our STDERR and STDOUT files were
        # written to, then appends them to a location within a hashed directory
        # structure based on $sub->id and deletes the original files. It adds
        # its own clear marker (which includes the seek points for every
        # previous marker) so you can easily grab just the last STD* for an
        # attempt to complete this submission, including being able to handle
        # appending empty output. If a file is over 10000 lines long, it only
        # appends the first and last 5000 lines and throws away the rest.
        
        # at this point you can now get at the STD* files:
        my $o_file = $sub->stdout(); # and ->stderr()
        my $fh = $sub->last_stdout(); # and ->last_stderr()
        # the above returns a filehandle to the $o_file opened and seeked to the
        # start of the last block of stdout
        my $fh = $sub->tail_stdout(); # and ->tail_sterr()
        # above returns a filehandle to a pipe of tail -50 (by default, the
        # method accepts an int to override this) on the $o_file.
    }
    else {
        # we must be pending in the scheduler
        my $pend_time = $sub->pend_time;
        # this also works when we're running/finished by taking the difference
        # between $job->start_time and $self->scheduled. Otherwise its the
        # difference between the time now and $self->scheduled.
    }
 }
 else {
    # no action has been taken with this submission yet (or we've been reset).
    # When you want to use a submission to submit to a scheduler, you first
    # claim it:
    if ($sub->claim) {
        # calling claim() sets a boolean in the db to true, which causes it
        # to return false on subsequent calls, but this call returns true.
        
        # based on the requirments, you submit the job to your scheduler and
        # get back some kind of id from it, which is then stored:
        my $received_sid = $sub->submit;
        # this submit method is just a convienence pass-through to
        # $sub->scheduler->submit(submission => $sub) which stores the id it
        # gets from the job scheduling system like this:
        $sub->sid($received_sid);
        # sid cannot be set to undef or '' or 0. When set, it calls ->release.
        # sid cannot be set unless we have the claim.
        
        $sub->release;
        # this sets the claim boolean back to false. claim() will still return
        # false though since it checks for sid() being set. ->release will also
        # be called automatically if $sub is destroyed without sid() getting
        # set.
        
        # at this point your job is scheduled in your scheduler, and the $sub
        # is effectivly 'locked', preventing other processes from submitting it
        # again
    }
    else {
        # Another process has claimed it; sid() cannot be set by us.
        # claim() will return false if claim has already been called and sid()
        # has not yet been set, and it will also return false if sid() or
        # scheduled() have been set. You'll normally just silently ignore this
        # $sub and move on dealing with other subs you may be dealing with,
        # trusting that some other process will take care of this $sub.
    }
 }

 # This class extends expunge() to delete the ->stdout and ->stderr files (and
 # parent dirs if empty) prior to removal from the db.

=head1 VertRes::PersistentArray

 This class represents an array of Persistent objects. It is implemented by
 storing the object ids and class in it's "table" in the Persistent db against
 an auto-increment self id. It can even store other PersistentArray objects.

 An example usage is a list of Submission objects that you want to submit
 all at once to your Scheduler for efficiency reasons. Though it is
 (deliberatly) not enforced, the idea is that you'd group together many
 Submission objects that all share the same Requirements object, and submit
 those as an array.

 # create one using a required array ref of Persistent objects:
 my $array = VertRes::PersistentArray->new(members => [$sub1, $sub2...]);
 # or by id:
 $array = VertRes::PersistentArray->new(id => $a_valid_array_id);

 # get a reference to the array entry in the db:
 my $array_id = $array->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # id auto-increments every time new(members => []) is called, so if called
 # twice in a row with the exact same list of persistents it will give back
 # 2 different PersistentArray objects with 2 different ids.

 # get a list of Persistent objects from the PersistentArray:
 my @persistents = $array->members;
 my @submissions = $array->members('VertRes::JobManager::Submission');
 # (members is read-only)

 # get an individual Persistent object from the PersistentArray:
 my $persistent = $array->member($index);
 
 # merge a list of Persistent objects with an old PersistentArray, creating
 # a new PersistentArray and expunging the old one:
 my $new_array = $old_array->merge(@obj_list);

=head1 VertRes::JobManager::Scheduler*

 A Scheduler is a class that presents a simple consistent interface to a
 particular job scheduler, such as LSF. There will also be a Scheduler called
 'Local', implementing the interface for a single-CPU system that has no job
 scheduling system - useful for running tests.
 (This class is implemented using VertRes::Persistent.)

 # create; no options necessary since it will take defaults from SiteConfig, but
 # the key args are:
 my $sch = VertRes::JobManager::Scheduler->new(type => 'LSF',
                                               output_root => '/abs/path');
 # or by id:
 $sch = VertRes::JobManager::Scheduler->new(id => $a_valid_scheduler_id);
 # these new() calls create VertRes::JobManager::Schedulers::LSF objects

 # get a reference to the scheduler entry in the db:
 my $sch_id = $sch->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # same id is always returned for every Scheduler object created with the
 # same type and output_root. Unlike most other persistent methods, these are
 # the only 2 variables we store persistently; most methods are more like class
 # methods and you supply all the data needed by them yourself from other
 # objects, and results are stored in those objects.

 # get info about a scheduler:
 my $type = $sch->type; # eg. the string 'LSF'
 my $output_root = $sch->output_root; # the absolute path to the root directory
                                      # that the user wants STDOUT & STDERR of
                                      # this scheduler stored

 # submit one or more jobs to the scheduler (queue them to be run):
 my $scheduled_id = $sch->submit(submission => $vertres_jobmanager_submission);
 # this gets requirements from the Submission object. It also auto-sets
 # $submission->sid() to $scheduled_id (the return value, which is the value
 # returned by the scheduler).
 $scheduled_id = $sch->submit(array => $vertres_persistentarray,
                              requirements => $vertres_jobmanager_requirements);
 # for running more than one job in an array, you pass a PersistentArray object
 # and a Requirments object (ie. where you have arranged that all the Submission
 # objects in the array share this same Requirements). It gets each Submission
 # object from the PersistentArray and auto-sets $submission->sid() to the
 # $scheduled_id with the PersistentArray id and index suffixed in a special
 # format understood by other Scheduler methods.
 # In both cases, it claims all Submission objects prior to interacting with
 # the scheduler, and if the scheduler has a problem and the submission fails,
 # we release all the Submission objects.

 # submit() uses the Requirements object to decide on a queue to submit to:
 my $queue = $sch->queue($vertres_jobmanager_requirements);
 # and it figures out exactly where STDOUT & STDERR should go based on hashing
 # an id (the ->id of the Submission or PersistentArray object):
 my $output_dir = $sch->output_dir($id); # a dir within ->output_root()
 # and it generates whatever command-line args correspond to the Requirements,
 # queue and output_dir:
 my $args_string = $sch->args($vertres_jobmanager_requirements);
 # When dealing with a PersistentArray, each Submission will get its own output
 # files in the output_dir() based on which index it was in the PersistentArray.
 # You access a particular Submission's output files with:
 my $o_file = $sch->stdout($vertres_jobmanager_submission);
 # and likewise with ->stderr. It works out the hashing-id and index from the
 # Submission->sid (for PersistentArray submits) or Submission->id (for single
 # Submission submits). Having dealt with the output files, you can delete
 # them:
 $sch->unlink_output($vertres_jobmanager_submission);
 # this also prunes empy dirs.

 # the command that submit() actually schedules in the scheduler is a little
 # perl -e that takes the referenced Submission id (working it out from the
 # PersistentArray object and index if this was an PersistentArray submit),
 # pulls out the Job and does ->run on that.

 # though some schedulers may have some kind of limit in place on the maximum
 # number of jobs a user or the system as a whole can keep scheduled at once,
 # we deliberatly don't model or support that limit here. If we do go over the
 # limit we treat it as any other error from the scheduler: ->submit would just
 # return false and ->last_submit_error would give the error as a string. The
 # user is then free to try the submit again later. To help the user decide if
 # we're under the limit or not before attempting a submit:
 my $max_jobs = $sch->max_jobs;
 # this is the get-only max jobs as set in SiteConfig
 my $current_jobs = $sch->queued_jobs;
 # this is the total number of jobs being tracked in the scheduling system for
 # the current user.
 # For testing purposes, VertRes::JobManager::Schedulers::Local will accept
 # a practically-infinite-sized list (max_jobs is v.high), and handle tracking
 # and serially running the jobs one-at-a-time itself.

 # query and manipulate a particular submission (these methods access the
 # scheduler itself, so should normally be avoided; use Submission and Job
 # methods to track state instead):
 if ($sub->scheduled) {
    # $sub->scheduled means that $sub has a ->sid, which is what the following
    # methods extract to work the answers out. It does intellegent caching of
    # data extracted from the scheduler for a given sid, so all these repeated
    # calls may only result in ~1 system call.
    
    if ($sch->pending($sub)) {
        my $current_queue = $sch->queue($sub);
    }
    elsif ($sch->running($sub)) {
        # perhaps you want to clear out a submission that you know from the
        # job itself has already finished, but the scheduler has gotten
        # 'stuck':
        $sch->kill($sub);
        # this tries a normal kill and waits for the ->sid to be cleared, and
        # failing that tries again with a more severe kill/clear if the
        # scheduler supports such a thing. Eventually it will return boolean
        # depending on if the kill was successful (->sid is no longer visible
        # in the list of sids presented by the scheduler).
        
        # or perhaps you've noticed the run time is approaching the limit you
        # set in Requirements, and you want to switch queues before the
        # scheduler kills your job:
        $sch->switch_queues($sub);
        # this method recalcluates the queue given the current Requirments
        # object associated with $sub, and if this is different from ->queue,
        # then it will attempt a queue switch and return 1 if successful, 0 if
        # no switch was necessary, and -1 if the switch failed.
    }
    elsif ($sch->finished($sub)) {
        # these methods fall back on reading the output files on disc if the
        # scheduler has forgotten about the sid in question:
        if ($sch->done($sub)) { # yay! }
        elsif ($sch->failed($sub)) {
            my $exit_code = $sch->exit_code($sub);
            if ($sch->ran_out_of_memory) { ... }
            elsif ($sch->ran_out_of_time) { ... }
            elsif ($sch->killed) {
                # Submission finished due to something like ->kill being called;
                # ie. the user killed the job deliberly, so a resubmission could
                # work
            }
            else {
                # Submission ended due to the Job crashing for some other
                # unknown reason. You should make sure $sub->update_status() got
                # called, and then use Submission methods to investigate the
                # STDOUT and STDERR files.
            }
        }
    }
 }

=head1 VertRes::PipelineManager*

 A "pipeline" is a series of "actions" (specified in a VertRes::Pipelines::*
 module) that are run on some part of a dataset, and PipelineManager* modules
 exist to make it easy to supply the correct data to the pipeline, run multiple
 instances of the pipeline at once on every part of the whole dataset at once
 (where actual running of commands makes use of VertRes::JobManager*), keep
 track of progress and make sure that everything completes successfully.

 A pipeline essentially boils down to an ordered series of
 VertRes::JobManager::Submission objects where each one will only run once
 certain other Submissions have completed and input files are available etc.

=head1 VertRes::PipelineManager::DataSource

 This class describes a source of data that a VertRes::Pipelines::* module
 could be run on. A datasource is ultimately comprised of a list of things
 (directories, files, database keys, whatever), and a particular instance of
 a VertRes::Pipelines::* module will run on a single element from that list.
 An instance of this class holds the information and methods necessary to
 generate these elements as required. (It is implemented using
 VertRes::Persistent.)
 For example, a datasource might be a list of 1000 bam files, and the
 information stored in an instance of this class might be the path to a file
 that contains the list of those 1000 bam files. An element might be 1 of those
 bam files, and a particular VertRes::Pipelines::* module might be capable of
 doing something with a single bam file. Another example would be chaining two
 pipelines together, where the second pipeline utilises a datasource that is the
 output of the first. The information stored in an instance of the DataSource
 would be the Setup of the first pipeline, and an element might be one of the
 files that the last action of the first pipeline provides.

 # create; there are 3 required args 'type', 'source', and 'method', and an
 # optional 'options' hash that can be set:
 my $dat = VertRes::PipelineManager::DataSource->new(type => 'VRTrack',
                source => 'a_vrtrack_database_name',
                method => 'unmapped_lanes',
                options => {platforms => ['SLX', '454']});
 # another example:
 $dat = VertRes::PipelineManager::DataSource->new(type => 'Fofn',
                source => '/abs/path/to/my.fofn',
                method => 'all_files');
 # and another example:
 $dat = VertRes::PipelineManager::DataSource->new(type => 'Chain',
                source => $vertres_pipelinemanager_setup_obj,
                # source for type => 'Chain' can also be a PersistentArray
                # instance that lists multiple Setup (see later) instances, to
                # make one Setup dependent upon multiple other Setups.
                method => 'last_action')
 # These create, respectively, a VertRes::PipelineManager::DataSources::VRTrack
 # object, a VertRes::PipelineManager::DataSources::Fofn object and a
 # VertRes::PipelineManager::DataSources::Chain object.
 # Or create using an id:
 $dat = VertRes::PipelineManager::DataSource->new(id => $a_valid_datasource_id);

 # get a reference to the datasource entry in the db:
 my $dat_id = $pip->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # same id is always returned for every DataSource object created with the
 # same type, source and method.

 # get the core info on this pipeline (these are get-only):
 my $type = $dat->type;
 my $source = $dat->source; # this might be a database name string or a file
                            # path etc.
 my $method = $dat->method; # this is a string corresponding to the name of a
                            # method defined in $dat's class, and is therefore
                            # $type specific

 # get/set the optional arbitrary settings that control how the desired method
 # will behave:
 $dat->options(platforms => ['SLX', '454'], ...); # you set a hash
 my %options = $dat->options; # and get back a hash

 # get the list of data elements:
 my @data_elements = $dat->data_elements;
 # this method causes $dat->$method($source, %options) to be run, and by
 # interface convention, these methods return a list of hash refs, where each
 # hash_ref has at least the key 'dataelement_key' (with a value unique to all
 # elements in the source, that will stay unique as the source grows) and can be
 # used as the instance-specific-input to at least 1 VertRes::Pipelines::*
 # module.

=head1 VertRes::PipelineManager::Setup

 This class describes a pipeline run on a certain dataset with certain settings.
 The settings are stored in a db (implemented using VertRes::Persistent). Before
 you can run a pipeline, you "set up a pipeline" by creating an instance of
 Setup.

 # create; there are 3 required key-forming options 'name', 'module' and
 # 'data_source', and new() will also accept an abitrary set of key => value
 # pairs that describe the settings needed by the pipeline module in question:
 my $set = VertRes::PipelineManager::Setup->new(name => 'chimp_mapping',
                module => 'VertRes::Pipelines::Mapping',
                data_source => $vertres_pipelinemanager_datasource,
                root => '/abs/working/dir',
                reference => 'chimp_ref.fa',
                other_mapping_setting => 'other_mapping_value',
                list_config => ['val1', 'val2'],
                hash_config => { key1 => 'val1', key2 => 'val2' }, ...);
 # or by id:
 $set = VertRes::PipelineManager::Setup->new(id => $a_valid_pipeline_id);

 # get a reference to the Setup entry in the db:
 my $set_id = $set->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # same id is always returned for every Setup object created with the same
 # name, module and data_source.

 # get the core info on this Setup (these are get-only):
 my $name = $set->name; # (a string)
 my $module = $set->module; # (a string)
 my $data_source = $set->data_source; # (a VertRes::PipelineManager::DataSource)

 # deal with settings:
 my $ref = $set->get('reference');
 $set->set('reference', 'newchimp_ref.fa'); # overwrite or create as appropriate
 $set->unset('reference'); # the reference key and any value is removed from the
                           # db entirely

 # set a list/hash of values:
 $set->set('samples', ['NA01', 'NA02']);
 $set->set('platforms', {'NA01' => 'SLX', 'NA02' => '454'});
 # and get back:
 my @samples = $set->get('samples');
 my %sample_to_platform = $set->get('platforms');

 # copy the arbitrary settings of one Setup into another:
 $set->copy($set2);
 # this does a "union" preferring $set2 values, leaving $set values for keys
 # not found in $set2 untouched.

 # pull out all the settings into a hash ref, for passing into a
 # VertRes::Pipelines::* module:
 my $settings = $set->settings;

 # a Setup can store on itself the concept of being 'active'; when inactive
 # other systems can choose to ignore this pipeline (this is not enforced in
 # any way):
 $set->active($boolean);
 if ($set->active) { # true by default }

 # it can also store on itself a PersistentArray, which can be used
 # to describe its "children" (objects created in order to run this Setup):
 my $vertres_persistentarray = $set->children($vertres_persistentarray);
 # the primary purpose of doing this is to be able to later clean out the db
 # of all objects created for this Setup, and easily access all associated
 # Submission objects.

=head1 VertRes::ResultStore

 This class makes it very easy to store (a small amount of) data from the output
 of some job in the database, instead of in a file on disc. Implemented using
 VertRes::Persistent.

 # create; 'name' being required and you'll have to arrange that your name is
 # going to be globally unique across all pipelines, setups and data sources
 # (*** how?!... this module needs more thought...):
 my $res = VertRes::ResultStore->new(name => '/abs/path/to/file');
 # or by id:
 $res = VertRes::ResultStore->new(id => $a_valid_resultstore_id);

 # get a reference to the resultstore entry in the db:
 my $res_id = $res->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # same id is always returned for every ResultStore object created with the
 # same name.

 # store something:
 $res->store($my_resultfile_contents_as_one_big_string);
 $res->store(%my_results_as_a_hash);
 $res->store(@my_results_as_a_list);
 
 # get it back:
 my $string = $res->retrieve;

=head1 VertRes::Pipeline*

 VertRes::Pipeline is the base class of the modules (VertRes::Pipelines::*) that
 actually define the work to be done to go from raw input to final output, as
 a series of 'actions'. These modules do not directly store any persistent data.
 The base class provides the following implemented methods:

 # find out what VertRes::PipelineManager::DataSource types and methods are
 # supported:
 my @valid = $obj->valid_datasources; # a list of array refs, each ref being a
                                      # type,method tuple.

 # what are the required and optional settings for this particular
 # VertRes::Pipelines::* module?
 my %configurable_settings = $obj->valid_settings;
 # configurable_settings is some hash with required and optional keys and then
 # some data structure to describe what is valid

 # get an ordered list of action names:
 my @actions = $obj->actions;

 # deal with an action:
 foreach my $action_name (@actions) {
    my @required_files = $obj->required_files($action_name);
    my @provided_files = $obj->provided_files($action_name);
    
    # if you want to run the action, you'll check your own records to see what
    # the state of this action is for this particular pipeline setup and
    # datasource element. You won't check if the required files exist, since you
    # want to be able to skip past this action if it has completed but you since
    # deleted its required (or provided) files because they were only temporary
    # to make a later action work. If your records say we haven't started
    # running this action, we do:
    $obj->run_action($action_name);
    
    # prior to derefrencing and running the action subroutine, run_action does:
    my @missing = $obj->missing_required_files($action_name);
    # if any files are missing, then run_action returns false and does nothing.
    # Otherwise it runs the desired subroutine. If the subroutine returns true,
    # that means the action completed and we can move on to the finish_action
    # step. If it returned false, it probably used dispatch() (see below) so we
    # see what it dispatched:
    my @dispatches = $obj->dispatched($action_name);
    # this returns a list of [$cmd_line_string, $requirements_object] refs that
    # you could use to create VertRes::JobManager::Submission objects and
    # eventually actually get the cmd_line to run. You'll keep records on how
    # those Submissions do, and when they're all successfull you'll do:
    $obj->finish_action($action_name);
    # this runs a post-processing method for the action (that is not allowed
    # to do any more dispatching) and returns true if the post-process was fine.
    # finish_action appends its own auto-check that all the provided files
    # exist. So if finish_action returns true you record in your system that
    # this action (on this pipeline setup and datasource element) finished, so
    # the next time you come to this loop you'll skip and try the following
    # action. (VertRes::PipelineManager::Manager does all this sort of stuff for
    # you.)
 }

 # child classes don't override actions(), required_files() etc., they just
 # define a class array ref called $actions which contains hash refs which
 # contain name, action, requires and provides keys, with the last 3 having
 # values of subroutine refs. 'finish' is another optional key, which also has a
 # subroutine ref as a value. valid_datasources() gets its info from a different
 # class array ref called $valid_datasources. valid_settings() gets its info
 # from a class hash ref called $valid_settings;

 # when you create a new instance of a VertRes::Pipeline subclass, you must
 # suppy the following:
 my $subclass = VertRes::Pipelines::subclass->new(
                    input => $hashref_from_a_datasource_element,
                    setup => $hashref_from_pipelinemanager_pipeline_setup);
 # given this information the subclass should then be able to work out what its
 # required files are, where it will write its provided files, and be able to
 # proceed with doing work when run_action($name) is called.


 # The action, requires, provides and finish subroutines of a subclass (which
 # receive no arguments directly) are implemented with the help of some
 # VertRes::Pipeline methods:

 # get the root working directory:
 my $root_dir = $self->root_dir; # == $self->{config}->{root} || cwd
 # if you work in a subdirectory of the root, perhaps elements from your
 # datasource are hashrefs with a 'rel_path' key pointing to the relative
 # directory; set that:
 $self->relative_dir($self->{input}->{rel_path});
 my $working_dir = $self->working_dir; # root_dir and relative_dir combined
 my $out_file = $self->working_file('my.out.file');
 # $out_file is a VertRes::File with the absolute path to my.out.file in
 # working_dir(). Whenever dealing with files generally you'll use
 # VertRes::File.

 # inside an action subroutine, when you want to do something that takes more
 # than a few seconds, you make a VertRes::JobManager::Requirements object
 # that describes the needs of what you want to do, and come up with a command
 # line that would actually do what you want to do. These get passed to:
 $self->dispatch(cmd => $cmd_line_string, requirments => $requirements_object);
 # nothing actually happens as a direct result of doing that. The caller must
 # subsequently get all the dispatched things using ->dispatched and do
 # something with the cmd and requirements.
 # if you dispatched things you return false, otherwise if you successfully
 # completed everything your action needed to do, you return true.

 # if an action needs to carry out some post-processing after the things it
 # dispached have actually been run, it can setup another subroutine and include
 # a ref to it under the 'finish' key of the $actions class ref.

 # you could make use of VertRes::ResultStore to store (small) final output
 # instead of writing to disc.


 # NB: compared to the original VertRes::Pipeline* system, this is deliberatly
 # very lightweight: no checking of file existance is done here to determine
 # if an action is complete, there is no accessing of LSF, there is no locking
 # on-disc, and there is no actual running of shell commands. Higher-level
 # VertRes::PipelineManager* stuff is supposed to handle all that. These modules
 # just define what should be done and try and do as little work themselves as
 # possible.
 # Also, the idea is to change each VertRes::Pipelines::* module so that it does
 # a single well defined thing (so VertRes::Pipelines::Mapping does mapping
 # only, not recalibration), and to chain different Setups together. So you'd
 # create a Setup where module => 'VertRes::Pipelines::Mapping' and data_source
 # => $a_datasource_with_type_vrtrack, and then you'd create another Setup where
 # module => 'VertRes::Pipelines::Recalibration' and it's data_source would be
 # of type 'chain' with a source of the first Setup. 

=head1 VertRes::PipelineManager::Action

 This class provides state tracking for running a particular
 VertRes::Pipelines::* module action with a certain pipeline setup against a
 specific VertRes::PipelineManager::DataSource->data_elements element. State is
 stored persistently thanks to implementing with VertRes::Persistent.

 # create; there are 4 required key-forming options:
 my $act = VertRes::PipelineManager::Action->new(action => 'name_of_action',
                module => 'VertRes::Pipelines::Mapping',
                dataelement_key => $dataelement_key,
                setup => $vertres_pipelinemanager_setup_obj);
 # (where $dataelement_key is the value of the 'dataelement_key' key in one of
 # the hash refs returned by the data_elements() method on an instance of a
 # VertRes::PipelineManager::DataSource object)
 # or by id:
 $act = VertRes::PipelineManager::Action->new(id => $a_valid_action_id);

 # get a reference to the action entry in the db:
 my $act_id = $act->id;
 # gives us a unique id that other systems and dbs can refer to or store. The
 # same id is always returned for every Action object created with the
 # same action, module, dataelement_key and setup.

 # get the core info on this Action (these are get-only):
 my $module = $act->module;
 my $action = $act->action; # a string name of an action in $module
 my $dataelement_key = $act->dataelement_key;
 my $setup = $act->setup; # a VertRes::PipelineManager::Setup instance

 # do stuff depending on the overall state of the Action:
 if ($act->complete) {
    # at some point we decided that this action completed fully successfully.
    # maybe we changed our minds and want to force it to restart?
    $act->reset;
    # this attempts reset() on every associated Submission (killing etc. first
    # if necessary), and then sets the complete boolean back to false. This
    # reset() can actually be called at any point, letting us force-stop an
    # action that we realised was doing something wrong, for example.
 }
 else {
    # get the Submission objects for this Action:
    my @submissions = $act->submissions;
    if (@submissions) {
        if ($act->all_ok) {
            # all_ok returns true when every Submission object says that it
            # completed ok. Now you can do whatever post-processing the action
            # might need, and finally:
            $act->complete(1);
        }
        else {
            # take a look at your Submission objects yourself and do stuff...
            # most likely there's still some running, so you just wait. For
            # convience, you can grab just the ones that are having problems:
            my @fails = $act->failed_submissions;
        }
    }
    else {
        # you didn't yet run your action, get the ->dispatched things, turn them
        # into VertRes::JobManager::Submission objects and associate them with
        # $act yet. Do that now:
        # ...
        $act->submissions(@submission_objects);
        # this turns the submission objects into a VertRes::PersistentArray
        # and stores that on itself, ready to turn them back into a list of
        # objects when submissions() is called.
    }
 }

=head1 VertRes::PipelineManager::Manager

 The high-level manager interface for script and end-user use. It does need to
 store one bit of persistent data, but it isn't instance-specific, so while this
 module makes use of VertRes::Persistent, it isn't like all the others that
 inherit from it. (Infact, it probably implements its persistent storage by
 abusing VertRes::ResultStore ?)

 # create:
 my $man = VertRes::PipelineManager::Manager->new;

 # get a list of all VertRes::PipelineManager::Setup objects:
 my @setups = $man->setups;
 # limited to those for a certain pipeline module:
 @setups = $man->setups(module => 'VertRes::Pipelines::Mapping');

 # trigger/find out about the pipelines:
 foreach my $set (@setups) {
    if ($set->active) {
        # run the pipeline (if it is complete, nothing happens; since new data
        # could have turned up, there is no general concept of a pipeline ever
        # automatically being considered 'complete' so there is no "first check
        # we haven't already completed" method at this level):
        my $num_of_jobs_still_to_go = $man->run($set);
        # This first does VertRes::Persistent::capture_start and when done will
        # call capture_stop, and will turn the array into a PersistentArray
        # which it will store in $set->children, merging the array with any
        # existing PersistentArray. This capture_stop&merge code will also run
        # during destruction of $man to protect against being killed before
        # run() completes normally. Under normal circumstances it will also
        # touch() every child of $pip, to prevent expunge_old destroying our
        # lives.
        # It also calls $self->running(1) which sets a boolean to true in the
        # persistent storage, along with the time we called this. Before
        # returning from run(), or on destruction, we set running(0) if the
        # timestamp was ours. If another process comes along and calls ->run
        # while we are running(), then the other process will block for some
        # reasonable amount of time until running is 0. If it takes too long,
        # then it will stop blocking and update the timestamp on running() and
        # continue with its own run().
        #
        # Next, it uses $set->data_source (having checked it is one of the
        # $set->module's valid ones) to loop through each data input. For each
        # data_source element it creates an instance of the $set->module. It
        # uses the ->actions method on that to set up
        # VertRes::PipelineManager::Action objects and uses those to track state
        # as it does ->run_action($name) etc. Note that this process is only
        # creating Action and Submission objects; we are not calling ->submit
        # on the Submissions, so run() doesn't itself directly cause anything
        # to be executed on the farm.
        #
        # As it loops through everything it keeps track of various things like
        # number of data_source elements, actions completed and still to go per
        # element and overall, and number of jobs
        # completed/running/to-be-run/failed per action, per element and
        # overall. It returns the total number of jobs still to go (those
        # currently running, or scheduled to run later), so if this is 0 then
        # this pipeline is finished for now at least.
        #
        # If run($pip, stdout => 1) is called, then it outputs this kind of
        # stuff in a pretty way to STDOUT as it goes along. Otherwise this info
        # is queryable after run($pip) returns:
        my %status = $man->status($set);
        
        # calling run() should be very cheap, since for the most part it will
        # just be doing lookups into the Persistent db where our state is
        # stored. However, if an action method has to actually run, then it may
        # not be fast enough for constructing a live-updated status webpage. For
        # that purpose, ->status($set) can be called without calling ->run($set)
        # first, in which case it will get all the stats it can by going through
        # the same loops as run(), except that it won't call ->run_action (so
        # will be able to say that an action has not yet started, but not how
        # many jobs there are to go for that action). ->status($set,
        # overall_jobs_only => 1) takes a shortcut and just directly grabs the
        # Submission objects from $set->children and generates a simple set of
        # completed/running/to-be-run/failed overall stats. NB: we never know
        # about actions that we haven't gotten to, so the total number of jobs
        # is never accurate until the last action has been run on all data
        # elements.
        
        # a simple hash from status isn't powerful/easy enough to allow things
        # like showing failed jobs and error messages in a web-frontend. For
        # that purpose, we have a set of methods that would let a user drill
        # down and discover what's going wrong with their pipeline:
        my @failed_actions = $man->failed_actions($set);
        # a 'failed' action is one that has submissions that have finished but
        # not ok.
        foreach my $action (@failed_actions) {
            # meta-information (like the dataelement_key) can be extracted in
            # the normal way from an Action for display
            
            my %outputs = $man->failed_outputs($action);
            # this gets the last STDOUT and STDERR from every failed Submission
            # associated with $action and returns a hash with Submission->id
            # keys and {stdout => 'string', stderr => 'string'} values.
            
            # you could imagine that a web-frontend would have a button that
            # called:
            $action->reset;
            # that the user might use if they looked at the outputs and fixed
            # the problem.
        }
        
        # perhaps the user fixed a problem and now wants to reset all the
        # actions that failed in one go:
        $man->reset(@failed_actions);
        
        # perhaps something really stupid and wrong happened with a pipeline and
        # you just want to start everything over completely from scratch:
        $man->reset($set);
        # this grabs all the Action objects from $set->children and does ->reset
        # on them.
        
        # run-time/memory-usage summary stats can be found on a per-action
        # basis, averaged over every pipeline that ran that action, and also
        # per pipeline:
        my %stats = $man->action_stats(); # averaged over all pipelines
        %stats = $man->action_stats($set); # for just this pipeline config
        # %stats == (action_name => { walltime => $seconds, memory => $mb });
    }
    elsif (time() - $set->last_updated > 7776000) {
        # it's been inactive for ~3months; perhaps we want to trash old stuff?
        $man->expunge($set);
        # expunges all Persistent entries in the Peristent db that were created
        # during the course of this pipeline setup, that are not used by another
        # pipeline setup
    }
 }

 # Having ->run at least 1 pipeline setup, we'll have a bunch of Submission
 # objects in the Persistent db, but none of them will have been submitted. We
 # deal with submitting, coping with failed Submissions, retrying Submissions
 # etc. in the following few methods, which don't need to know anything about
 # the pipelines that spawned them. This is completly generic, and you could
 # imagine running a script that calls these methods every 10mins in a cron job
 # or something.

 # get all Submissions that have not completed ok (used internally by the
 # subsequent methods, answer cached once per instance of Manager):
 my @incomplete_submissions = $man->submissions();
 # again, this is all incomplete Submission objects, regardless of what pipeline
 # or action or process spawned them. submissions() will block if the running()
 # boolean is true, to give a run() time to finish, so we'll see all of its
 # Submission objects. After a reasonable amount of time it will stop blocking
 # though, incase a run() went bad and failed to turn off running().

 # deal with failed Submissions:
 $man->resubmit_failures(max_retries => 3);
 # this checks all the ->submissions for ones that failed, and uses the standard
 # Submission and Job methods to work out why and resubmit them as appropriate,
 # potentially with updated Requirements. max_retries defaults to 3.

 $man->check_running;
 # this does the dance of checking if any of the currently running ->submissions
 # are approaching their time limit, and switching queues as appropriate. It
 # also checks that all running Jobs have had a recent heartbeat, and if not
 # will do the kill dance and resubmit.

 # batch Submissions into arrays and submit them:
 $man->submit;
 # this groups all the non-permanently-failed ->submissions into lists based
 # on shared Requirements objects, creates [PersistentArray, Requirements] for
 # each group, and uses those to do a Scheduler->submit for each group. If a
 # group consists of only 1 Submission, then it won't bother with the creating
 # a PersistentArray step. Suppling ->submit($percent) makes it consider only
 # that percentage of incomplete Submissions. This would be useful in a script
 # that would be run by, say, 4 team members: the script would say ->submit(25),
 # and the workload would be shared amongst the members, helpful when dealing
 # with fair-share systems in schedulers.
 
 # expunge all the PersistentArray objects that submit() has made in the past,
 # but are now defunct (they reference Submissions that have all finished OK):
 $man->cleanup;

=cut

package pipeline_ideas;
use strict;
use warnings;

1;
