ALTER TABLE `assembly` MODIFY `md5` char(30);
ALTER TABLE `image` MODIFY `name` varchar(255);
ALTER TABLE `library_type_id` MODIFY `name` varchar(255);
ALTER TABLE `mapper` MODIFY `name` varchar(255);
ALTER TABLE `population` MODIFY `name` varchar(255);
ALTER TABLE `study` MODIFY `name` varchar(255);
ALTER TABLE `sample` MODIFY `name` varchar(255);
ALTER TABLE `seq_centre` MODIFY `name` varchar(255);
ALTER TABLE `seq_tech` MODIFY `name` varchar(255);
ALTER TABLE `submission` MODIFY `name` varchar(255);

DROP TABLE IF EXISTS `autoqc`;
CREATE TABLE `autoqc`
(
  `autoqc_id` mediumint(8) unsigned NOT NULL auto_increment,
   mapstats_id mediumint(8) unsigned NOT NULL DEFAULT 0,
   test varchar(50) NOT NULL default '',
   result smallint(5) unsigned NOT NULL DEFAULT 0,
   reason varchar(200) NOT NULL default '',
   PRIMARY KEY (`autoqc_id`),
  KEY  `mapstats_id` (`mapstats_id`),
   UNIQUE KEY `mapstats_test` (`mapstats_id`, `test`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

update schema_version set schema_version=20;
