#!/bin/sh
umask 002
ROOT="/lustre/scratch105"
CONF="/nfs/vertres01/conf"
DUMPS="/warehouse/g1k-04/sql_dumps"
export LD_LIBRARY_PATH=/software/badger/lib:/software/oracle_client-10.2.0/lib
export ORACLE_HOME=/software/oracle_client-10.2.0

date="`date +'%y%m%d'`"
mysqldump -u $VRTRACK_RW_USER -p$VRTRACK_PASSWORD -P$VRTRACK_PORT -h$VRTRACK_HOST vrtrack_vbseq > "$DUMPS/vrtrack_vbseq_$date.sql"

update_vrtrack.pl --database vrtrack_vbseq --projects SEQCAP_WGS_High-powered_complex_trait_association_mapping_through_whole_genome_sequencing_of_a_selected_subpopulation_of_the_INGI-Val_Borbera_genetic_isolate   --create_individuals --no_fastq > "$ROOT/log/vrtrack_update_vbseq.log" 2> "$ROOT/log/vrtrack_update_vbseq.err"

