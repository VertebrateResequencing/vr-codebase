# vr-codebase
A collection of modules and scripts for processing and analysing large quantities of sequencing data.

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/vr-codebase/blob/master/LICENSE)

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Required dependencies](#required-dependencies)
    * [Optional dependencies](#optional-dependencies)
    * [Vrtrack database](#vrtrack-database)
  * [Usage](#usage)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)

## Introduction
This repository contains the internally-developed software used by the Vertebrate Resequencing group at the Sanger Institute. It comprises mostly self-documented Perl code. There are both scripts and modules within subfolders. Each module has its own POD, so please use perldoc or similar for further help.

For example:
```
perldoc VertRes::Utils::Sam
```
An overview, guides and how-tos are available on the [GitHub wiki](https://github.com/VertebrateResequencing/vr-codebase/wiki).

## Installation
vr-codebase has the following dependencies:

### Required dependencies

 * samtools
 * Filesys::DfPortable
 * Filesys::DiskUsage
 * File::Fetch
 * File::Rsync
 * File::Temp
 * Net::FTP::Robust
 * Time::Format
 * IO::Capture::Stderr
 * Math::Random
 * Bio::Assembly::Improvement
 * Bio::AutomatedAnnotation
 * Bio::RNASeq
 * Bio::VertRes::Permissions

### Optional dependencies

 * GATK
 * picard
 * beagle

For samtools you need the source version compiled with -fPIC and -m64 in the CFLAGS, and the environment variable SAMTOOLs pointing to that source directory (which should now contain bam.h and libbam.a).

It is recommended that you set PERL_INLINE_DIRECTORY to ~/.Inline

Some dependencies need environment variables setup (using setenv in csh or export in bash). The following table shows the name of the software, the environment variable you need to set, and the value you should set it to.

|Software|Environment variable|Value                                     |
|--------|--------------------|------------------------------------------|
|samtools|SAMTOOLS            |/path/to/samtools/source_directory        |
|GATK    |GATK                |/path/to/gatk_jar_files                   |
|GATK    |STING_DIR           |/path/to/gatk_source_code_checkout        |
|GATK    |GATK_RESOURCES      |/path/to/resource_files_like_reference_etc|
|picard  |PICARD              |/path/to/picard_jar_files                 |
|beagle  |BEAGLE              |/path/to/beagle_install_directory         |

Details for the installation are provided below. If you encounter an issue when installing vr-codebase please contact your local system administrator. If you encounter a bug please email us at path-help@sanger.ac.uk.

To build vr-codebase:
```
perl Build.PL
```
If this says you have "ERRORS/WARNINGS FOUND IN PREREQUISITES", install missing prerequisites from CPAN:
```
./Build installdeps
```
To test the code prior to use:
```
perl Build.PL

./Build test
```
To install:
```
./Build install
```
Most likely some of the tests will fail due to you not having certain external software installed. If you don't plan on making use of that software, just ignore it when a test script for that software fails.

### Vrtrack database
Tracking of meta-data, required for many of our pipelines, occurs in a mysql database. We use environement variables to define how to access the mysql database:

 * VRTRACK_HOST
 * VRTRACK_PORT
 * VRTRACK_RO_USER
 * VRTRACK_RW_USER
 * VRTRACK_PASSWORD

The RW user you setup should have permissions to create and alter databases. We assume that the RO user does not require a password.

## Usage
```
Usage: run-pipeline [OPTIONS]
Options:
   -c, --config <file>             The file containing the projects.
   -d, --done                      List of completed jobs (i.e. jobs with the status 'done')
   -f, --failed                    List jobs which failed repeatedly (i.e. jobs with the status 'failed').
   -l, --logfile <file>            Where to put logging stuff.
   -L, --lockfile <file>           Run only if there is no other process running.
   -m, --max-lsf-jobs <int>        Allow at most this many LSF jobs in queue.
   -o, --once                      Run all jobs only once and then exit.
   -s, --sleep-time <int>          Check the projects every N minutes.
   -t, --todo                      List jobs which should be run (i.e. jobs with the status 0,1,2).
   -a, --dont_check_load           Dont halt script if the load is higher than the number of CPUs.
   -b, --same_database_for_all     All config files refer to the same DB so reuse connection.
   -v, --verbose                   Be verbose. (Given twice increases the verbosity level.)
```

For more information, please see the [GitHub wiki](https://github.com/VertebrateResequencing/vr-codebase/wiki).

## License
Vr-codebase is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/vr-codebase/blob/master/LICENSE)

## Feedback/Issues
Please report any issues to path-help@sanger.ac.uk.
