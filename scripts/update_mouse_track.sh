#!/bin/sh
umask 002
ROOT="/lustre/scratch102"
WGS_PROJ_ROOT="$ROOT/projects/vrtrack_mouse"
EXOME_PROJ_ROOT="$ROOT/projects/vrtrack_mouse_exome"
export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

export PATH=/software/bin:/usr/local/bin:/usr/bin:/bin
export PERL5LIB=/software/vertres/modules:/software/vertres/lib/all

#backup databases
date="`date +'%y%m%d'`"
mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST  vrtrack_mouse | gzip -c > "$WGS_PROJ_ROOT/sql_dumps/vrtrack_mouse_$date.sql.gz"
mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST  vrtrack_mouse_exome | gzip -c > "$EXOME_PROJ_ROOT/sql_dumps/vrtrack_mouse_exome_$date.sql.gz"

#WGS update the individuals and samples files
/software/vertres/scripts/generateIndividualSamplesMap.pl -s $WGS_PROJ_ROOT/meta-data/studies.fofn -m $WGS_PROJ_ROOT/meta-data/individuals2Samples.tab -spe Mus_musculus -t 10090 -nm $WGS_PROJ_ROOT/meta-data/new_individuals2Samples.tab -ns $WGS_PROJ_ROOT/meta-data/samples.tab 2> "$ROOT/log/sample_map_vrtrack_mouse.out"
mv $WGS_PROJ_ROOT/meta-data/individuals2Samples.tab $WGS_PROJ_ROOT/meta-data/individuals2Samples.tab_old
mv $WGS_PROJ_ROOT/meta-data/new_individuals2Samples.tab $WGS_PROJ_ROOT/meta-data/individuals2Samples.tab

#EXOME update the individuals and samples files
/software/vertres/scripts/generateIndividualSamplesMap.pl -s $EXOME_PROJ_ROOT/meta-data/studies.fofn -m $EXOME_PROJ_ROOT/meta-data/individuals2Samples.tab -spe Mus_musculus -t 10090 -nm $EXOME_PROJ_ROOT/meta-data/new_individuals2Samples.tab -ns $EXOME_PROJ_ROOT/meta-data/samples.tab 2> "$ROOT/log/sample_map_vrtrack_mouse_exome.out"
mv $EXOME_PROJ_ROOT/meta-data/individuals2Samples.tab $EXOME_PROJ_ROOT/meta-data/individuals2Samples.tab_old
mv $EXOME_PROJ_ROOT/meta-data/new_individuals2Samples.tab $EXOME_PROJ_ROOT/meta-data/individuals2Samples.tab

#load the individuals (if there are new ones)
/software/vertres/scripts/load_vrtrack_individuals.pl --indiv $WGS_PROJ_ROOT/meta-data/samples.tab -db vrtrack_mouse 2> "$ROOT/log/load_ind_vrtrack_mouse.out"
/software/vertres/scripts/load_vrtrack_individuals.pl --indiv $EXOME_PROJ_ROOT/meta-data/samples.tab -db vrtrack_mouse_exome 2> "$ROOT/log/load_ind_vrtrack_mouse_exome.out"

#update the database
/software/vertres/scripts/update_vrtrack.pl --projects $WGS_PROJ_ROOT/meta-data/studies.fofn --database vrtrack_mouse --sample_map $WGS_PROJ_ROOT/meta-data/individuals2Samples.tab 2> "$ROOT/log/update_vrtrack_mouse.out"
/software/vertres/scripts/update_vrtrack.pl --projects $EXOME_PROJ_ROOT/meta-data/studies.fofn --database vrtrack_mouse_exome --sample_map $EXOME_PROJ_ROOT/meta-data/individuals2Samples.tab 2> "$ROOT/log/update_vrtrack_mouse_exome.out"

#/software/vertres/scripts/update_vrtrack.pl --spp mouse --projects $DISK_ROOT/conf/mouse_projects  > "$DISK_ROOT/log/vrtrack_update_mouse.log" 2> "$DISK_ROOT/log/vrtrack_update_mouse.err"
