#!/bin/sh
umask 002
ROOT="/lustre/scratch102"
PROJ_ROOT="$ROOT/projects/vrtrack_mouse"
export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

export PATH=/software/bin:/usr/local/bin:/usr/bin:/bin
export PERL5LIB=/software/vertres/modules:/software/vertres/lib/all

date="`date +'%y%m%d'`"
mysqldump -u vreseq_rw -pt3aml3ss -P3306 -hmcs4a vrtrack_mouse | gzip -c > "$PROJ_ROOT/sql_dumps/vrtrack_mouse_$date.sql.gz"

#update the individuals and samples files
/software/vertres/bin/generateIndividualSamplesMap.pl -s $PROJ_ROOT/meta-data/studies.fofn -m $PROJ_ROOT/meta-data/individuals2Samples.tab -spe Mus_musculus -t 10090 -nm $PROJ_ROOT/meta-data/new_individuals2Samples.tab -ns $PROJ_ROOT/meta-data/samples.tab 2> "$ROOT/log/sample_map_vrtrack_mouse.out"

mv $PROJ_ROOT/meta-data/individuals2Samples.tab $PROJ_ROOT/meta-data/individuals2Samples.tab_old
mv $PROJ_ROOT/meta-data/new_individuals2Samples.tab $PROJ_ROOT/meta-data/individuals2Samples.tab

#load the individuals (if there are new ones)
/software/vertres/bin/load_vrtrack_individuals.pl --indiv $PROJ_ROOT/meta-data/samples.tab -db vrtrack_mouse 2> "$ROOT/log/load_ind_vrtrack_mouse.out"

#update the database
/software/vertres/bin/update_vrtrack.pl --projects $PROJ_ROOT/meta-data/studies.fofn --database vrtrack_mouse --sample_map $PROJ_ROOT/meta-data/individuals2Samples.tab 2> "$ROOT/log/update_vrtrack_mouse.out"

#/software/vertres/bin/update_vrtrack.pl --spp mouse --projects $DISK_ROOT/conf/mouse_projects  > "$DISK_ROOT/log/vrtrack_update_mouse.log" 2> "$DISK_ROOT/log/vrtrack_update_mouse.err"
