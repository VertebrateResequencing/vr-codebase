ALTER TABLE assembly MODIFY COLUMN taxon_id INT(10) unsigned DEFAULT NULL;
ALTER TABLE file MODIFY COLUMN file_id INT(10)  unsigned NOT NULL DEFAULT '0';
ALTER TABLE image MODIFY COLUMN mapstats_id INT(10)  unsigned NOT NULL DEFAULT '0';
ALTER TABLE individual MODIFY COLUMN species_id INT(10) unsigned DEFAULT NULL;
ALTER TABLE lane MODIFY COLUMN seq_request_id INT(10)  unsigned NOT NULL DEFAULT '0';
ALTER TABLE library MODIFY COLUMN ssid INT(10) unsigned DEFAULT NULL;
ALTER TABLE mapstats MODIFY COLUMN mapstats_id INT(10)  unsigned NOT NULL DEFAULT '0';
ALTER TABLE sample MODIFY COLUMN ssid INT(10) unsigned DEFAULT NULL;
ALTER TABLE species MODIFY COLUMN species_id INT(10) unsigned NOT NULL AUTO_INCREMENT;
ALTER TABLE species MODIFY COLUMN taxon_id INT(10) NOT NULL DEFAULT '0';

update schema_version set schema_version=32;
