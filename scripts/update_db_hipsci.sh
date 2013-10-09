#!/bin/bash
umask 002

function usage
{
    echo "usage: update_db_hipsci.sh [[-d vrtrack_db_name -f file_type -g output_location (optional: -t tax_id ) ]  | [-h]]"
}

FILE=""
DB=""
TAX=""
GSF=""
while [ "$1" != "" ]; do
    case $1 in
        -d | --db )             shift
                                DB=$1
                                ;;
        -t | --tax )            shift
        						TAX="-tax $1"
                                ;;
        -f | --file )           shift
        						FILE="-f $1"
                                ;;
        -g | --gsf )            shift
        						GSF="-gsf $1"
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

if [ "$DB" = "" ];then
    echo "A valid tracking database name must be entered."
    exit
fi

DBEXISTS=$(mysql -u $VRTRACK_RO_USER -h$VRTRACK_HOST --batch --skip-column-names -e "SHOW DATABASES LIKE '$DB'" | grep $DB > /dev/null; echo "$?")
if [ $DBEXISTS -eq 1 ];then
    echo "A database with the name $DB does not exist."
    exit
fi

ARG_UP="-u -v"

CONF="/nfs/vertres01/conf"
BIN_EXT="/software/vertres/update_pipeline_hipsci"
DUMPS="/warehouse/g1k-04/sql_dumps/$DB.sql"

export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST $DB > $DUMPS

$BIN_EXT/update_pipeline.pl -s $CONF/$DB"_studies" -d $DB $TAX $FILE $GSF
