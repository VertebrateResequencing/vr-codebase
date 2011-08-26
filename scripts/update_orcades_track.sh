#!/bin/sh
umask 002
ROOT="/lustre/scratch105"
CONF="/nfs/vertres01/conf"
DUMPS="/warehouse/g1k-01/sql_dumps"
export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

date="`date +'%y%m%d'`"
mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST vrtrack_orcades > "$DUMPS/vrtrack_orcades_$date.sql"

update_vrtrack.pl --database vrtrack_orcades --projects SEQCAP_WGS_ORCADES  --create_individuals --no_fastq > "$ROOT/log/vrtrack_update_orcades.log" 2> "$ROOT/log/vrtrack_update_orcades.err"

