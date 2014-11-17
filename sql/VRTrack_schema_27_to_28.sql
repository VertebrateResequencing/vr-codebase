ALTER TABLE `file` ADD `reference` varchar(255) DEFAULT NULL;
DROP VIEW if EXISTS `latest_file`;
create view latest_lane as select * from file where latest=true;
update schema_version set schema_version=28;