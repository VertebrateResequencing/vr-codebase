#!/usr/bin/env perl
use strict;
use warnings;
use Time::Format;
use Getopt::Long;
use CGI;

# get user input
my ($help, $dir);
my $mode = 'users';
GetOptions('mode=s' => \$mode,
           'dir=s' => \$dir,
           'h|help' => \$help);

unless ($mode eq 'users' || $mode eq 'projects') {
    warn "bad --mode; must be users|projects, not '$mode'\n";
    $help = 1;
}

($help) and die <<USAGE;
Gets tasks out of Asana. 

Usage: $0 --mode projects
        --mode [users|projects] In users mode (default), shows tasks by user,
                                and so does not show tasks they have archived.
                                Tasks marked for 'later' are also ignored.
                                It outputs html with headers, suitable for a web
                                server. If the user_name variable is posted to
                                the script, limits to that user.
                                
                                In projects mode, shows tasks by project, not
                                affected by what users have archived. It outputs
                                in tsv format.
        --dir [path]            If supplied, and if in projects mode, instead of
                                outputting to stdout, outputs to a compressed
                                file in the supplied directory, name after
                                today's date.
        --help                  (This message)
        
The environment variable ASANA_API must be set, providing the authorization
credentials. Note that this is currently hard-coded for the Vertebrate
Resequencing workspace id.

USAGE

my $auth = $ENV{ASANA_API} || die "ASANA_API environment variable not set\n";
my $base_url = 'https://app.asana.com/api/1.0';
my $workspace_id = 870050815366;

if ($mode eq 'projects') {
    my $ofh;
    if ($dir && -d $dir) {
        open($ofh, "| gzip -9 > $dir/$time{'yyyy-mm-dd'}.gz");
        select($ofh);
    }
    
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
    
    if ($ofh) {
        close($ofh);
    }
}
else {
    my $q = CGI->new;
    #print $q->header();
    my $desired_user_name = $q->param('user_name');
    
    my $css = <<'CSS';
body { 
    font-size: 12pt; 
    color: black; 
    background: white; 
    margin: 10px; 
    padding: 0px;
    font-family: sans-serif;
}

h1, h2, h3 {
    margin: 0px;
    padding: 0px;
}
h1 {
    margin-top: 40px;
}
h2 {
    margin-top: 30px;
    margin-left: 20px;
}
h3 {
    margin-top: 20px;
    margin-left: 40px;
}

.task_details {
    padding: 10px;
    margin-left: 60px;
    margin-bottom: 5px;
    border-bottom: 1px solid #C0C0C0;
}

.project {
    font-size: 10pt;
    color: #555753;
    border-right: 1px solid grey;
    padding-right: 6px;
    margin-right: 3px;
}

.notes {
    font-size: 10pt;
    color: #909090;
    margin-top: 2px;
}
CSS
    
    print $q->start_html(-title => 'Asana tasks '.($desired_user_name ? "for $desired_user_name" : 'by user'),
                         -style => { -code => $css });
    
    # we go by tasks of each user instead of all tasks, so that we only see what
    # that user hasn't archived
    my $users = asana("workspaces/$workspace_id/users");
    
    my %tasks;
    foreach my $user_ref (@$users) {
        my $user_name = $user_ref->{name};
        if ($desired_user_name) {
            next unless $user_name eq $desired_user_name;
        }
        
        # get all this user's tasks
        my $tasks = asana("tasks?workspace=$workspace_id&assignee=$user_ref->{id}");
        
        # get details of each task
        foreach my $t_ref (@$tasks) {
            my $task_name = $t_ref->{name};
            next unless $task_name;
            next if $task_name =~ /:$/; # sub sections
            my $d = asana("tasks/$t_ref->{id}");
            my $assignee = $d->{assignee} ? $d->{assignee}->{name} : '';
            $assignee || next;
            $assignee eq $user_name || next;
            my $project_name = $d->{projects}->[0]->{name};
            my $notes = $d->{notes} || '';
            
            push(@{$tasks{$assignee}->{$d->{completed} ? 'completed' : 'todo'}->{$d->{assignee_status} || 'inbox'}}, [$project_name, $task_name, $notes]);
        }
    }
    
    foreach my $user (sort keys %tasks) {
        print $q->h1($user), "\n";
        foreach my $status (sort keys %{$tasks{$user}}) {
            my $hash = $tasks{$user}->{$status};
            print "\t", $q->h2($status), "\n";
            
            if ($status eq 'completed') {
                my @completed_tasks;
                foreach my $array (values %$hash) {
                    push(@completed_tasks, @$array);
                }
                foreach my $task_details (sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] } @completed_tasks) {
                    html_task_output($q, $task_details);
                }
            }
            else {
                foreach my $as_status ('inbox', 'today', 'upcoming') {
                    defined $hash->{$as_status} || next;
                    print "\t\t", $q->h3($as_status), "\n";
                    
                    foreach my $task_details (sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1] } @{$hash->{$as_status}}) {
                        html_task_output($q, $task_details);
                    }
                }
            }
        }
    }
    
    print $q->end_html, "\n";
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
        warn "no data retrieved for [$curl]; is asana down or overloaded?\n";
    }
    return $data;
}

sub string_to_hash {
    my $str = shift;
    $str =~ s/:null/:undef/g;
    $str =~ s/:true/:1/g;
    $str =~ s/:false/:0/g;
    $str =~ s/":/" => /g;
    $str =~ s/\B\$(\w)/\\\$$1/g;
    my $hash = eval $str;
    return $hash;
}

sub html_task_output {
    my $q = shift;
    my ($project_name, $task_name, $notes) = @{$_[0]};
    print "\t\t\t", $q->div({class => 'task_details'}, $q->span({class => 'project'}, $project_name), $task_name, $q->div({class => 'notes'}, $notes)), "\n";
}