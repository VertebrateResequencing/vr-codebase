#!/usr/bin/env perl
use strict;
use warnings;

my $auth = $ENV{ASANA_API};
my $base_url = 'https://app.asana.com/api/1.0';
my $workspace_id = 870050815366;
my @fields = qw(project_name task_name assignee created_at modified_at due_on assignee_status completed_at completed notes);

# get the list of projects
my $projects = asana("workspaces/$workspace_id/projects");

# for each project, get the tasks
print join("\t", @fields), "\n";
foreach my $p_ref (@$projects) {
    my $project_name = $p_ref->{name};
    my $tasks = asana("projects/$p_ref->{id}/tasks");
    
    # output details of each task
    foreach my $t_ref (@$tasks) {
        my $task_name = $t_ref->{name};
        next unless $task_name;
        next if $task_name =~ /:$/; # sub sections
        my $d = asana("tasks/$t_ref->{id}");
        my $assignee = $d->{assignee} ? $d->{assignee}->{name} : '';
        my $notes = $d->{notes} || '';
        $notes =~ s/\n/ /g;
        print join("\t", $project_name, $task_name, $assignee, map { $d->{$_} || '' } @fields[3..8], $notes), "\n";
    }
}

exit;

sub asana {
    my $command = shift;
    my $curl = qq[curl -s -u "$auth:" "$base_url/$command"];
    
    my $data;
    my $max_retries = 4;
    while (! $data) {
        my $return = `$curl`;
        my $hash = string_to_hash($return);
        $data = $hash->{data};
        unless ($data) {
            $max_retries--;
            last if $max_retries == 0;
            sleep(1);
        }
    }
    
    unless ($data) {
        warn "no data retrieved for [$curl]; is asano down or overloaded?\n";
    }
    return $data;
}

sub string_to_hash {
    my $str = shift;
    $str =~ s/:null/:undef/g;
    $str =~ s/:true/:1/g;
    $str =~ s/:false/:0/g;
    $str =~ s/":/" => /g;
    my $hash = eval $str;
    return $hash;
}