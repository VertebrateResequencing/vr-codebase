#!/usr/bin/perl -w
use strict;
use warnings;
use Parallel::ForkManager;

BEGIN {
    use Test::Most;
    eval {
        require VRTrack::Testconfig;
    };
    if ($@) {
        plan skip_all => "Skipping all tests because VRTrack tests have not been configured";
    }
    else {
        plan tests => 14;
    }
    
    use_ok('VRTrack::VRTrack');
    use_ok('VRTrack::Lane');
    use_ok('VRTrack::File');
}

my $connection_details = { database => VRTrack::Testconfig->config('test_db'),
                           host     => VRTrack::Testconfig->config('host'),
                           port     => VRTrack::Testconfig->config('port'),
                           user     => VRTrack::Testconfig->config('user'),
                           password => VRTrack::Testconfig->config('password'), 
                        };

# access the db as normal with new()
ok my $vrtrack = VRTrack::VRTrack->new($connection_details), 'VRTrack new() returned something';
isa_ok $vrtrack, 'VRTrack::VRTrack';

# create 1000 lanes with 1 file each
my @lane_ids;
my @file_ids;
foreach my $lid (1..1000) {
    my $lane = VRTrack::Lane->create($vrtrack, 'lane_'.$lid);
    my $file = $lane->add_file('file_'.$lid.'.bam');
    $file->type(4);
    $file->update;
    $lane->update;
    push(@file_ids, $file->id);
    push(@lane_ids, $lane->id);
}
is scalar(grep { defined $_ } @lane_ids, @file_ids), 2000, 'created 1000 lanes and files to test with';

# fork to 15 processes and try to update all the file table entries: this should
# cause deadlocks due to the simultaneous accesses
my $fm = Parallel::ForkManager->new(15);
my %retrieved_responses = ();
$fm->run_on_finish (
    sub {
        my $data_structure_reference = $_[5];
        if (defined($data_structure_reference)) {
            $retrieved_responses{$$data_structure_reference} = 1;
        }
    }
);

foreach my $file_id (@file_ids) {
    $fm->start and next; # fork
    
    $vrtrack = VRTrack::VRTrack->new($connection_details);
    
    # also test that we can nest transactions (the update() call also calls
    # transaction())
    my $worked = $vrtrack->transaction(sub {
        my $file = VRTrack::File->new($vrtrack, $file_id);
        $file->is_processed('import' => 1);
        $file->update;
        my $lane = VRTrack::Lane->new($vrtrack, $file->lane_id);
        $lane->is_processed('import' => 1);
        $lane->update;
    });
    
    $fm->finish(0, \$worked); # exit in the child process
}
$fm->wait_all_children;

is check_updated('import'), 2000, 'managed to update all 1000 lanes and files, having done 15 updates simultaneously in a nested transaction';
my @results = keys %retrieved_responses;
my $all_worked = 0;
if (@results == 1 && $results[0] == 1) {
    $all_worked = 1;
}
is $all_worked, 1, 'transaction() claimed that all 1000 worked';

# test that nested transactions work as expected - dieing midway through the
# top transaction should result in updates in the first half that were
# successful not actually happening
%retrieved_responses = ();
foreach my $file_id (@file_ids) {
    $fm->start and next;
    
    $vrtrack = VRTrack::VRTrack->new($connection_details);
    
    my $worked = $vrtrack->transaction(sub {
        my $file = VRTrack::File->new($vrtrack, $file_id);
        $file->is_processed('mapped' => 1);
        $file->update;
        die "mid transaction\n";
        my $lane = VRTrack::Lane->new($vrtrack, $file->lane_id);
        $lane->is_processed('mapped' => 1);
        $lane->update;
    });
    
    my $failed_with_expected_error = 0;
    unless ($worked) {
        $failed_with_expected_error = 1 if $vrtrack->{transaction_error} eq "Transaction failed, rolled back. Error was: mid transaction";
    }
    
    $fm->finish(0, \$failed_with_expected_error);
}
$fm->wait_all_children;

is check_updated('mapped'), 0, 'deliberately dieing during a nested transaction results in nothing getting updated';
@results = keys %retrieved_responses;
my $all_failed_with_expected_error = 0;
if (@results == 1 && $results[0] == 1) {
    $all_failed_with_expected_error = 1;
}
is $all_failed_with_expected_error, 1, 'transaction() claimed that none of them worked, and we got a suitable error message for them';

# make sure we're compatable with the old transaction_start and _commit methods
foreach my $file_id (@file_ids) {
    $fm->start and next;
    
    $vrtrack = VRTrack::VRTrack->new($connection_details);
    
    $vrtrack->transaction_start();
    
    my $sig_warn = $SIG{'__WARN__'};
    $SIG{'__WARN__'} = sub { };
    
    eval {
        my $file = VRTrack::File->new($vrtrack, $file_id);
        $file->is_processed('mapped' => 1);
        $file->update;
        my $lane = VRTrack::Lane->new($vrtrack, $file->lane_id);
        $lane->is_processed('mapped' => 1);
        $lane->update;
        
        $vrtrack->transaction_commit();
    };
    if ($@) {
        eval { $vrtrack->transaction_rollback(); };
    }
    
    $SIG{'__WARN__'} = $sig_warn;
    
    $fm->finish(0);
}
$fm->wait_all_children;
cmp_ok check_updated('mapped'), '>', 10, 'we were able to get at least some updates to work when using transaction_start/commit with transaction() nested; deadlocks can still occur this way';

foreach my $file_id (@file_ids) {
    $fm->start and next;
    
    $vrtrack = VRTrack::VRTrack->new($connection_details);
    
    $vrtrack->transaction_start();
    
    my $sig_warn = $SIG{'__WARN__'};
    $SIG{'__WARN__'} = sub { };
    
    eval {
        my $file = VRTrack::File->new($vrtrack, $file_id);
        $file->is_processed('improved' => 1);
        $file->update;
        die "mid transaction\n";
        my $lane = VRTrack::Lane->new($vrtrack, $file->lane_id);
        $lane->is_processed('improved' => 1);
        $lane->update;
        
        $vrtrack->transaction_commit();
    };
    if ($@) {
        eval { $vrtrack->transaction_rollback(); };
    }
    
    $SIG{'__WARN__'} = $sig_warn;
    
    $fm->finish(0);
}
$fm->wait_all_children;
is check_updated('improved'), 0, 'we got no updates to work when using transaction_start/commit with a die inside, and transaction() nested - transactional safety is still intact this way';

foreach my $file_id (@file_ids) {
    $fm->start and next;
    
    $vrtrack = VRTrack::VRTrack->new($connection_details);
    
    $vrtrack->transaction(sub {
        my $file = VRTrack::File->new($vrtrack, $file_id);
        $file->is_processed('improved' => 1);
        $file->update;
        
        $vrtrack->transaction_start();
        my $lane = VRTrack::Lane->new($vrtrack, $file->lane_id);
        $lane->is_processed('improved' => 1);
        $lane->update;
        $vrtrack->transaction_commit();
    });
    
    $fm->finish(0);
}
$fm->wait_all_children;
is check_updated('improved'), 2000, 'we got all updates to work when using transaction_start/commit nested inside a transaction()';

foreach my $file_id (@file_ids) {
    $fm->start and next;
    
    $vrtrack = VRTrack::VRTrack->new($connection_details);
    
    $vrtrack->transaction(sub {
        my $file = VRTrack::File->new($vrtrack, $file_id);
        $file->is_processed('snp_called' => 1);
        $file->update;
        
        $vrtrack->transaction_start();
        my $lane = VRTrack::Lane->new($vrtrack, $file->lane_id);
        $lane->is_processed('snp_called' => 1);
        $lane->update;
        $vrtrack->transaction_commit();
        
        die "force fail\n";
    });
    
    $fm->finish(0);
}
$fm->wait_all_children;
is check_updated('snp_called'), 0, 'we got no updates to work when using transaction_start/commit nested inside a transaction() that was forced to fail';


# cleanup
my $dbh = $vrtrack->{_dbh};
foreach ($dbh->tables()){
    next if /`latest_/;
    next if /schema_version/;
    $dbh->do("TRUNCATE TABLE $_");
}
exit;


sub check_updated {
    my $processed = shift;
    my $vrtrack = VRTrack::VRTrack->new($connection_details);
    my $num_updated = 0;
    foreach my $file_id (@file_ids) {
        $num_updated++ if VRTrack::File->new($vrtrack, $file_id)->is_processed($processed);
    }
    foreach my $lane_id (@lane_ids) {
        $num_updated++ if VRTrack::Lane->new($vrtrack, $lane_id)->is_processed($processed);
    }
    return $num_updated;
}
