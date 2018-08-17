ALTER TABLE `library` MODIFY `ssid` int(8) unsigned DEFAULT NULL;
ALTER TABLE `multiplex_pool` MODIFY `ssid` int(8) unsigned DEFAULT NULL;
ALTER TABLE `library_request` MODIFY `ssid` int(8) unsigned DEFAULT NULL;
ALTER TABLE `seq_request` MODIFY `ssid` int(8) unsigned DEFAULT NULL;
ALTER TABLE `project` MODIFY `ssid` int(8) unsigned DEFAULT NULL;
ALTER TABLE `study` MODIFY `ssid` int(8) unsigned DEFAULT NULL;
update schema_version set schema_version=27;
