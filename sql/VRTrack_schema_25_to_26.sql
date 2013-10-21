ALTER TABLE individual MODIFY individual_id int(10) unsigned NOT NULL auto_increment ;
ALTER TABLE allocation MODIFY individual_id int(10) unsigned NOT NULL DEFAULT '0';
ALTER TABLE sample MODIFY individual_id int(10) unsigned DEFAULT NULL;
ALTER TABLE assembly MODIFY assembly_id int(10) unsigned NOT NULL auto_increment ;
ALTER TABLE mapstats MODIFY assembly_id int(10) unsigned DEFAULT NULL;

ALTER TABLE library MODIFY sample_id int(10)  UNSIGNED NOT NULL;
ALTER TABLE library_request MODIFY sample_id int(10)  UNSIGNED NOT NULL;
ALTER TABLE sample MODIFY sample_id int(10)  UNSIGNED NOT NULL;

ALTER TABLE library MODIFY library_id int(10)  UNSIGNED NOT NULL;
ALTER TABLE library_multiplex_pool MODIFY library_id int(10)  UNSIGNED NOT NULL;
ALTER TABLE seq_request MODIFY library_id int(10)  UNSIGNED NOT NULL;
ALTER TABLE lane MODIFY library_id int(10)  UNSIGNED NOT NULL;

ALTER TABLE lane MODIFY lane_id int(10)  UNSIGNED NOT NULL;
ALTER TABLE mapstats MODIFY lane_id int(10)  UNSIGNED NOT NULL;
ALTER TABLE file MODIFY lane_id int(10)  UNSIGNED NOT NULL;

ALTER TABLE image MODIFY image_id int(10) unsigned NOT NULL auto_increment ;
update schema_version set schema_version=26;
