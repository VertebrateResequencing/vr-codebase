#!/bin/sh
umask 002
date="`date +'%y%m%d'`"
ROOT="/lustre/scratch105"
CONF="/nfs/vertres01/conf"
COMPRESS_CMD="/usr/bin/lzma --best --force "
SCRIPTS="/software/vertres/scripts"
BIN_EXT="/software/vertres/bin-external/update_pipeline"
DB="vrtrack_human_wes"
DUMPS="/warehouse/g1k-04/sql_dumps/"$DB"_"$date".sql"

export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST $DB > $DUMPS

bsub -M8000000 -R'select[mem>8000] rusage[mem=8000]' -o "$ROOT/log/compress_$DB.log" -e "$ROOT/log/compress_$DB.err" "$COMPRESS_CMD $DUMPS"

$BIN_EXT/update_pipeline.pl -s $CONF/$DB"_studies" -d $DB -v

$SCRIPTS/vrtrack_individual_supplier_name -d $DB
