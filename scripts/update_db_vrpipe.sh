#!/bin/sh
umask 002
if [ $# -lt 1 -o $# -gt 2  ]
then
    echo "Error in $0 - Invalid Argument Count"
    echo "Syntax: $0 vr-pipe_db_name(vrtrack_mouse_wgs, etc.) -OR-"
    echo "Syntax: $0 vr-pipe_db_name taxon_id (9606, 10090....)"
    exit
fi

DB="$1"
DBEXISTS=$(mysql -u $VRTRACK_RO_USER -h$VRTRACK_HOST --batch --skip-column-names -e "SHOW DATABASES LIKE '$DB'" | grep $DB > /dev/null; echo "$?")
if [ $DBEXISTS -eq 1 ];then
    echo "A database with the name $DB does not exist."
    exit
fi

TAX=""
if [ $# -eq 2  ]
then
    TAX="-tax $2"
fi

ARG_UP=""
if [[ $DB =~ "g1k_track_phase3" ]]
then
    ARG_UP="-u -sup -v"
elif [[ $DB =~ "g1k_track_phase2_vrpipe" ]] 
then
	ARG_UP="-u -sup -nop -v -min 9000"
else
    ARG_UP="-u -sup -nop -v"
fi

ROOT="/lustre/scratch105"
CONF="/nfs/vertres01/conf"
SCRIPTS="/software/vertres/scripts"
BIN_EXT="/software/vertres/bin-external/update_pipeline"
DUMPS="/warehouse/g1k-04/sql_dumps/$DB.sql"

export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST $DB > $DUMPS

$BIN_EXT/update_pipeline.pl -s $CONF/$DB"_studies" -d $DB $TAX $ARG_UP
