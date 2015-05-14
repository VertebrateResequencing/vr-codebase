ALTER TABLE species MODIFY species_id mediumint(8) unsigned NOT NULL auto_increment ;
ALTER TABLE individual MODIFY species_id mediumint(8) unsigned DEFAULT NULL ;

ALTER TABLE `lane` ADD `manually_withdrawn` tinyint(1) DEFAULT NULL;
DROP VIEW if EXISTS `latest_lane`;
create view latest_lane as select * from lane where latest=true;

update schema_version set schema_version=25;
