ALTER TABLE `project` ADD `data_access_group` varchar(255) DEFAULT NULL;
DROP VIEW if EXISTS `latest_project`;
create view latest_project as select * from project where latest=true;
update schema_version set schema_version=29;
