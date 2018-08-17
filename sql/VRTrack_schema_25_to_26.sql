ALTER TABLE `assembly` MODIFY `taxon_id` mediumint(8) unsigned DEFAULT 0;
ALTER TABLE `species` MODIFY `taxon_id` mediumint(8) unsigned DEFAULT 0;
ALTER TABLE `allocation` MODIFY `study_id` smallint(5) unsigned NOT NULL;
ALTER TABLE `allocation` MODIFY `individual_id` smallint(5) unsigned NOT NULL;
ALTER TABLE `allocation` MODIFY `seq_centre_id` smallint(5) unsigned NOT NULL;
update schema_version set schema_version=26;
