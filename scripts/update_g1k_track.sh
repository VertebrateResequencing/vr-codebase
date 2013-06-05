#!/bin/sh
umask 002
date="`date +'%y%m%d'`"
ROOT="/lustre/scratch105"
CONF="/nfs/vertres01/conf"
COMPRESS_CMD="/usr/bin/lzma --best --force "
STUDY="g1k_phase2"
DB="g1k_track_phase2"
DUMPS="/warehouse/g1k-04/sql_dumps/"$DB"_"$date".sql"

export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST $DB > $DUMPS

bsub -M8000 -R'select[mem>8000] rusage[mem=8000]' -o "$ROOT/log/compress_$DB.log" -e "$ROOT/log/compress_$DB.err" "$COMPRESS_CMD $DUMPS"

update_vrtrack.pl --database $DB --projects $CONF/$STUDY"_studies" --create_individuals --no_fastq > "$ROOT/log/update_$DB.log" 2> "$ROOT/log/update_$DB.err"
