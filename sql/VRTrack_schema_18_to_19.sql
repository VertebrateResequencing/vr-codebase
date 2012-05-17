ALTER TABLE `mapstats` ADD `is_qc` tinyint(1) DEFAULT '0';
ALTER TABLE `mapstats` ADD `prefix` varchar(40) DEFAULT '_';
DROP VIEW if EXISTS `latest_mapstats`;
create view latest_mapstats as select * from mapstats where latest=true;
update schema_version set schema_version=19;
update  `mapstats`  set is_qc = 1 where rmdup_bases_mapped is NULL;
