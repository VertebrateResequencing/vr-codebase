ALTER TABLE `individual` ADD INDEX `acc` (`acc`);
update schema_version set schema_version=34;