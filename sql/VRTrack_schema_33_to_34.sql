ALTER TABLE `individual` ADD INDEX `acc` (`acc`);
ALTER TABLE library ADD INDEX library_id_sample_id_latest (library_id,sample_id,latest);
ALTER TABLE lane ADD INDEX library_id_latest (library_id,latest);
update schema_version set schema_version=34;