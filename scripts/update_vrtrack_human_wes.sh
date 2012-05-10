#!/bin/sh
umask 002
DUMPS="/warehouse/g1k-04/sql_dumps"
CONF="/nfs/vertres01/conf"
SCRIPTS="/software/vertres/scripts"
BIN_EXT="/software/vertres/bin-external/update_pipeline"
COMPRESS_CMD="/usr/bin/lzma --best --force "
DB="vrtrack_human_wes"

export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

date="`date +'%y%m%d'`"
mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST $DB > $DUMPS/$DB"_"$date.sql

$COMPRESS_CMD $DUMPS/$DB"_"$date.sql

$BIN_EXT/update_pipeline.pl -s $CONF/$DB"_studies" -d $DB -v

$SCRIPTS/vrtrack_individual_supplier_name -d $DB
