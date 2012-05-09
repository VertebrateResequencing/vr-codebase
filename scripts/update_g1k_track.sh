#!/bin/sh
umask 002
ROOT="/lustre/scratch105"
CONF="/nfs/vertres01/conf"
DUMPS="/warehouse/g1k-04/sql_dumps"
export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

date="`date +'%y%m%d'`"
mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST g1k_track_phase2 > "$DUMPS/g1k_track_phase2_$date.sql"

update_vrtrack.pl --database g1k_track_phase2 --projects $CONF/g1k_phase2_studies --create_individuals --no_fastq > "$ROOT/log/vrtrack_update_g1k_phase2.log" 2> "$ROOT/log/vrtrack_update_g1k_phase2.err"

