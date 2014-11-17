DROP TABLE IF EXISTS db_projects; 
DROP TABLE IF EXISTS tracking_database;
DROP TABLE IF EXISTS schema_version; 
DROP TABLE IF EXISTS sample_id_mapping;
DROP TABLE IF EXISTS vrpipe_pipelinesetup;
DROP TABLE IF EXISTS vrpipe_file_info;
DROP TABLE IF EXISTS vrpipe_usage_top_level;
DROP TABLE IF EXISTS vrpipe_usage_total;
DROP TABLE IF EXISTS vrpipe_root_top_level;
DROP TABLE IF EXISTS pipeline_file;
DROP TABLE IF EXISTS file_info;

CREATE TABLE tracking_database(
  db_id INT AUTO_INCREMENT NOT NULL,
  db_name VARCHAR(100) NOT NULL,
  imported datetime NOT NULL,
  PRIMARY KEY (db_id),
  UNIQUE (db_name)
)ENGINE=InnoDB DEFAULT CHARSET=latin1;
ALTER TABLE tracking_database ADD INDEX tracking_db_id (db_id); 
ALTER TABLE tracking_database ADD INDEX tracking_db_name (db_name);

CREATE TABLE schema_version (
  schema_version mediumint(8) unsigned NOT NULL,
  imported datetime NOT NULL,
  PRIMARY KEY  (schema_version)
)ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE db_projects (
  id INT AUTO_INCREMENT NOT NULL,
  db_id INT NOT NULL,
  project_id smallint(5) unsigned NOT NULL,
  project_name varchar(255) NOT NULL,
  ssid mediumint(8) unsigned DEFAULT NULL,
  imported datetime NOT NULL,
  PRIMARY KEY (id),
  KEY db_id (db_id),
  KEY project_id (project_id),
  FOREIGN KEY(db_id) REFERENCES tracking_database(db_id) ON DELETE CASCADE
)ENGINE=InnoDB DEFAULT CHARSET=latin1;

insert into schema_version values (24, NOW());

CREATE TABLE sample_id_mapping(
  id INT AUTO_INCREMENT NOT NULL,
  db_id INT NOT NULL,
  db_name VARCHAR(100) NOT NULL,  
  project_id smallint(5) unsigned NOT NULL,
  project_name varchar(255) NOT NULL,
  supplier_name varchar(255) DEFAULT NULL,
  accession_number varchar(50) DEFAULT NULL,
  sanger_sample_name varchar(40) NOT NULL,
  PRIMARY KEY (id)
)ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE vrpipe_pipelinesetup (
  ps_id int(9) NOT NULL,
  ps_name varchar(64) NOT NULL,
  ps_user varchar(64) NOT NULL,
  ps_type varchar(64) NOT NULL,
  pipeline int(9) NOT NULL,
  output_root text NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE vrpipe_file_info (
  ps_id int(9) NOT NULL,
  ps_name varchar(64) NOT NULL,
  ps_user varchar(64) NOT NULL,  
  ps_type varchar(64) NOT NULL,
  step_number smallint(4) NOT NULL,
  file_id int(9) NOT NULL,
  s bigint(20) DEFAULT NULL,
  type VARCHAR(4) DEFAULT NULL,
  path VARCHAR(255) NOT NULL,
  path_root VARCHAR(128) DEFAULT NULL,
  mtime datetime DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE vrpipe_root_top_level (
  path_root VARCHAR(128) DEFAULT NULL,
  file_type VARCHAR(32) DEFAULT NULL,
  s bigint(40) DEFAULT NULL,
  s_gb bigint(30) DEFAULT NULL,
  file_count int(9) DEFAULT NULL,
  top_level_display tinyint(1) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE vrpipe_usage_top_level (
  ps_id int(9) NOT NULL,
  ps_name varchar(64) NOT NULL,
  ps_user varchar(64) NOT NULL,
  ps_type varchar(64) NOT NULL,
  path_root VARCHAR(128) DEFAULT NULL,
  file_type VARCHAR(4) DEFAULT NULL,
  s bigint(30) DEFAULT NULL,
  s_gb bigint(22) DEFAULT NULL,
  file_count int(9) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE vrpipe_usage_total (
  ps_id int(9) NOT NULL,
  ps_name varchar(64) NOT NULL,
  ps_user varchar(64) NOT NULL,
  ps_type varchar(64) NOT NULL,
  total_s bigint(30) DEFAULT NULL,
  total_s_gb bigint(22) DEFAULT NULL,
  file_count int(9) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE file_info (
  id int(9) NOT NULL,
  s bigint(20) DEFAULT NULL,
  type VARCHAR(4) DEFAULT NULL,
  path VARCHAR(255) NOT NULL,
  moved_to int(9) DEFAULT NULL,
  path_root VARCHAR(128) DEFAULT NULL  
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE pipeline_file (
  pipelinesetup int(9) NOT NULL,
  file int(9) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

ALTER TABLE file_info ADD INDEX file_info_id_idx (id);
ALTER TABLE pipeline_file ADD INDEX pipeline_file_pipelinesetup_file_idx (pipelinesetup,file);
