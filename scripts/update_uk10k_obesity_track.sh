#!/bin/sh
umask 002
ROOT="/lustre/scratch105"
CONF="/nfs/vertreseq01/conf"
DUMPS="/warehouse/g1k-01/sql_dumps"
export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

date="`date +'%y%m%d'`"
mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST vrtrack_uk10k_obesity > "$DUMPS/vrtrack_uk10k_obesity_$date.sql"

update_vrtrack.pl --database vrtrack_uk10k_obesity --projects $CONF/uk10k_obesity_studies --create_individuals --no_fastq > "$ROOT/log/vrtrack_update_uk10k_obesity.log" 2> "$ROOT/log/vrtrack_update_uk10k_obesity.err"

