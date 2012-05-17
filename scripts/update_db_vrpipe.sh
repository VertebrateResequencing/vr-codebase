#!/bin/sh
umask 002
if [ $# -ne 1 ]
then
    echo "Error in $0 - Invalid Argument Count"
    echo "Syntax: $0 vr-pipe_db_name (vrtrack_mouse_wgs, etc.)"
    exit
fi

DB="$1"
DBEXISTS=$(mysql -u vreseq_ro -hmcs4a --batch --skip-column-names -e "SHOW DATABASES LIKE '$DB'" | grep $DB > /dev/null; echo "$?")
if [ $DBEXISTS -eq 1 ];then
    echo "A database with the name $DB does not exist."
    exit
fi

date="`date +'%y%m%d'`"
ROOT="/lustre/scratch105"
CONF="/nfs/vertres01/conf"
SCRIPTS="/software/vertres/scripts"
BIN_EXT="/software/vertres/bin-external/update_pipeline"
DUMPS="/warehouse/g1k-04/sql_dumps/"$DB"_"$date".sql"

export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST $DB > $DUMPS

$BIN_EXT/update_pipeline.pl -s $CONF/$DB"_studies" -d $DB -v

$SCRIPTS/vrtrack_individual_supplier_name -d $DB
