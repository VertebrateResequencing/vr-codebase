#!/bin/sh
umask 002
ROOT="/lustre/scratch105"
CONF="/nfs/vertres01/conf"
DUMPS="/warehouse/g1k-04/sql_dumps"
COMPRESS_CMD="/usr/bin/lzma --best --force "
STUDY="uk10k_rare"
DB="vrtrack_$STUDY"
export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

date="`date +'%y%m%d'`"
mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST $DB > $DUMPS/$DB"_"$date.sql

$COMPRESS_CMD $DUMPS/$DB"_"$date.sql

update_vrtrack.pl --database $DB --projects $CONF/$STUDY"_studies" --create_individuals --no_fastq > "$ROOT/log/update_$DB.log" 2> "$ROOT/log/update_$DB.err"
