ALTER TABLE `mapstats` ADD `percentage_reads_with_transposon` FLOAT unsigned DEFAULT NULL;
DROP VIEW if EXISTS `latest_mapstats`;
create view latest_mapstats as select * from mapstats where latest=true;
ALTER TABLE seq_request MODIFY seq_type ENUM('Single ended sequencing','Paired end sequencing','HiSeq Paired end sequencing','MiSeq sequencing')  NOT NULL; 
DROP VIEW if EXISTS `latest_seq_request`;
create view latest_seq_request as select * from seq_request where latest=true;
ALTER TABLE `assembly` ADD `taxon_id` mediumint(8) unsigned DEFAULT NULL;
ALTER TABLE `assembly` ADD `translation_table` smallint(5) unsigned DEFAULT NULL;
update schema_version set schema_version=16;