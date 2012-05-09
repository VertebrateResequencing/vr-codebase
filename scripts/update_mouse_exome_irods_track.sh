#!/bin/sh
umask 002
ROOT="/lustre/scratch102/projects/vrtrack_mouse_exome_irods"
export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0
export PATH=/software/bin:/usr/local/bin:/usr/bin:/bin
export PERL5LIB=/software/vertres/modules:/software/vertres/lib/all

#backup databases
date="`date +'%y%m%d'`"
mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST  vrtrack_mouse_exome_irods | gzip -c > "$ROOT/sql_dumps/vrtrack_mouse_exome_irods_$date.sql.gz"

#EXOME update the individuals and samples files
/software/vertres/scripts/generateIndividualSamplesMap.pl -s $ROOT/meta-data/studies.fofn -m $ROOT/meta-data/individuals2Samples.tab -spe Mus_musculus -t 10090 -nm $ROOT/meta-data/new_individuals2Samples.tab -ns $ROOT/meta-data/samples.tab 2> "$ROOT/log/sample_map_vrtrack_mouse_exome_irods.out"
mv $ROOT/meta-data/individuals2Samples.tab $ROOT/meta-data/individuals2Samples.tab_old
mv $ROOT/meta-data/new_individuals2Samples.tab $ROOT/meta-data/individuals2Samples.tab

#load the individuals (if there are new ones)
/software/vertres/scripts/load_vrtrack_individuals.pl --indiv $ROOT/meta-data/samples.tab -db vrtrack_mouse_exome_irods 2> "$ROOT/log/load_ind_vrtrack_mouse_exome_irods.out"

#update the database
/software/vertres/scripts/update_vrtrack.pl --projects $ROOT/meta-data/studies.fofn --database vrtrack_mouse_exome_irods --sample_map $ROOT/meta-data/individuals2Samples.tab 2> "$ROOT/log/update_vrtrack_mouse_exome_irods.out"