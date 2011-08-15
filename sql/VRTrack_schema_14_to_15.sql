ALTER TABLE `mapstats` ADD `percentage_reads_with_transposon` FLOAT unsigned DEFAULT NULL;
DROP VIEW if EXISTS `latest_mapstats`;
create view latest_mapstats as select * from mapstats where latest=true;
update schema_version set schema_version=15;
