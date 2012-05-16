#!/bin/sh
umask 002
if [ $# -ne 1 ]
then
    echo "Error in $0 - Invalid Argument Count"
    echo "Syntax: $0 vrtrack_db_name (vrtrack_uk10k_obesity, vrtrack_ctcf, etc.)"
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
DUMPS="/warehouse/g1k-04/sql_dumps/"$DB"_"$date".sql"

export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST $DB > $DUMPS

update_vrtrack.pl --database $DB --projects $CONF/$DB"_studies" --create_individuals --no_fastq > "$ROOT/log/update_$DB.log" 2> "$ROOT/log/update_$DB.err"
