
=pod 

=head1 NAME

    Bio::EnsEMBL::Hive::PipeConfig::SimpleExample_conf;

=head1 SYNOPSIS

init_pipeline.pl Bio::EnsEMBL::Hive::PipeConfig::SimpleExample_conf -input_fofn my.fofn -append foo -chunk_size 10

=head1 DESCRIPTION

    This is the PipeConfig file acting as a simple example.
    It takes all the files in an input fofn, splits them into chunks of up to
    -chunk_size lines long each, alters each line in the chunks by appending
    -append, then cats the chunks together again.

    The SimpleExample pipeline consists of four "analyses" (types of tasks):
    'start', 'split', 'append' and 'cat'.

    * A 'start' job takes in a fofn and creates several jobs of the 'split'
      analysis (one for each file mentioned in the fofn).

    * A 'split' job takes in a file and splits it into -chunk_size chunks,
      recording state and stats to the db, creating an append job for each
      chunk it makes, planning to run a 'cat' job after those have completed.

    * An 'append' job takes a chunk file, appends -append onto the end of every
      line, checking for truncation vs db stats.

    * A 'cat' job waits for the appends of a split to complete before
      concatenating the chunks and checking for truncation, then deleting the
      chunk files.

=cut

package Bio::EnsEMBL::Hive::PipeConfig::SimpleExample_conf;

use strict;
use warnings;

# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

=head2 default_options

    Description : Implements default_options() interface method of
                  Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used
                  to initialize default options.
                  In addition to the standard things it defines two options,
                  'append' and 'chunk_size', leaving input_fofn as required
                  input from the user.

=cut

sub default_options {
    my ($self) = @_;
    
    my $opts = $self->SUPER::default_options(@_);
    $opts->{pipeline_name} = 'simple_example';
    $opts->{append} = 'default_append';
    $opts->{chunk_size} = 10;
    
    return $opts;
}

=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of
                  Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists
                  the commands that will create and set up the Hive database.
                  In addition to the standard creation of the database and
                  populating it with Hive tables and procedures it also creates
                  two pipeline-specific tables used by Runnables to communicate.

=cut

sub pipeline_create_commands {
    my ($self) = @_;
    return [
        # inherit database and hive tables' creation
        @{$self->SUPER::pipeline_create_commands},
        # additional tables needed for simple example pipeline's operation:
        'mysql '.$self->dbconn_2_mysql('pipeline_db', 1)." -e 'CREATE TABLE append_files (orig_file VARCHAR(200), max_lines INT, append VARCHAR(200), chunk INT, append_file VARCHAR(200), actual_lines INT, PRIMARY KEY (orig_file, max_lines, append, chunk))'",
        'mysql '.$self->dbconn_2_mysql('pipeline_db', 1)." -e 'CREATE TABLE completed (orig_file VARCHAR(200), max_lines INT, append VARCHAR(200), result_file VARCHAR(200), actual_lines INT, PRIMARY KEY (orig_file, max_lines, append))'"
    ];
}

=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of
                  Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines
                  the structure of the pipeline: analyses, jobs, rules, etc.
                  Here it defines four analyses:

                  * 'start' with one job (parse the input_fofn)
                    Will dataflow (create more jobs) via branch #1 into 'split'

                  * 'split' initially without jobs (they will flow from 'start')
                    Will dataflow via branch #2 into 'append' and via branch
                    #1 into 'cat'.

                  * 'append' initially without jobs (they will flow from
                    'split')

                  * 'cat' initially without jobs (they will flow from 'split').
                    All 'cat' jobs will wait for completion of *all* 'append'
                    jobs for the corresponding split before their own execution
                    (to ensure all data is available).

=cut

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {   -logic_name => 'start',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Start',
            -parameters => {},
            -input_ids => [
                { 'input_file_list' => $self->o('input_fofn'),
                  'lines_per_split' => $self->o('chunk_size'),
                  'string_to_append' => $self->o('append') },
            ],
            -flow_into => {
                1 => [ 'split' ], # will create a fan of jobs
            },
        },
        
        {   -logic_name    => 'split',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Split',
            -parameters    => {},
            -input_ids     => [
                # (jobs for this analysis will be flown_into via branch-1 from
                #  'start' job)
            ],
            -flow_into => {
                2 => [ 'append' ],
                1 => [ 'cat' ]
            },
        },
        
        {   -logic_name => 'append',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Append',
            -parameters => {},
            -input_ids => [
                # (jobs for this analysis will be flown_into via branch-2 from
                #  'split' jobs)
            ],
        },
        
        {   -logic_name => 'cat',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SimpleExample::Cat',
            -parameters => {},
            -input_ids => [
                # (jobs for this analysis will be flown_into via branch-1 from
                #  'split' jobs)
            ],
            # we can only start catting when all splits chunks have been
            # appended
            -wait_for => [ 'append' ],
        },
    ];
}

1;

