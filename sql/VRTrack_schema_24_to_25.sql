ALTER TABLE `lane` ADD `manually_withdrawn` tinyint(1) DEFAULT NULL;
DROP VIEW if EXISTS `latest_lane`;
create view latest_lane as select * from lane where latest=true;
update schema_version set schema_version=25;
