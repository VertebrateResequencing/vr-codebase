#!/bin/sh
umask 002
ROOT="/lustre/scratch105"
CONF="/nfs/vertreseq01/conf"
DUMPS="/warehouse/g1k-01/sql_dumps"
export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

#export PATH=/software/bin:/usr/local/bin:/usr/bin:/bin
#export PERL5LIB=/software/vertres/modules:/software/vertres/lib/all

date="`date +'%y%m%d'`"
mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST g1k_track > "$DUMPS/g1k_track_$date.sql"

update_vrtrack.pl --spp g1k --projects $CONF/g1k_projects > "$ROOT/log/vrtrack_update_g1k.log" 2> "$ROOT/log/vrtrack_update_g1k.err"

#grep -v 'updating' "$ROOT/log/g1k_update_$date"
