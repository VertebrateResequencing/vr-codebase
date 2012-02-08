ALTER TABLE species ADD KEY `name` (`name`) ;
UPDATE schema_version SET schema_version=18 ;
