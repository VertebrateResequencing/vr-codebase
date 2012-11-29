ALTER TABLE species MODIFY species_id mediumint(8) unsigned NOT NULL auto_increment ;
ALTER TABLE individual MODIFY species_id mediumint(8) unsigned DEFAULT NULL ;
update schema_version set schema_version=25;
