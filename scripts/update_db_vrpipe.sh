#!/bin/bash
umask 002

function usage
{
    echo "usage: update_db_vrpipe.sh [[-d vrtrack_db_name (optional: -t tax_id -m min_run_id -s study) ]  | [-h]]"
}

TAX=""
DB=""
MIN_RUN=""
STUDY=""
SPE=""
FORMAT=""
while [ "$1" != "" ]; do
    case $1 in
        -d | --db )             shift
                                DB=$1
                                ;;
        -t | --tax )            shift
        						TAX="-tax $1"
                                ;;
        -n | --species )        shift
                                                        SPE="-spe '$1'"
                                ;;
        -m | --min )            shift
        						MIN_RUN="-min $1"
                                ;;
        -s | --study )          shift
        						STUDY=$1
                                ;;
	-f | --format )		shift
							FORMAT="--file_type $1"
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

ARG_UP="-u -sup -nop -md5 -wdr -trd -v"

ROOT="/lustre/scratch105"
SCRIPTS="/software/vertres/scripts"
BIN_EXT="/software/vertres/update_pipeline"
DUMPS="/warehouse/g1k-04/sql_dumps/$DB.sql"

export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST $DB > $DUMPS

if [ "$STUDY" = "" ]; then
	cd $BIN_EXT && perl bin/update_pipeline.pl -s $CONF/$DB"_studies" -d $DB $TAX $SPE $MIN_RUN $FORMAT $ARG_UP
else
	cd $BIN_EXT && perl bin/update_pipeline.pl -n $STUDY -d $DB $TAX $SPE $MIN_RUN $FORMAT $ARG_UP
fi
