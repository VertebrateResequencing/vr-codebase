ALTER TABLE seq_request MODIFY seq_type ENUM(
'HiSeq Paired end sequencing'  ,
'Illumina-A HiSeq Paired end sequencing'  ,
'Illumina-A Paired end sequencing'  ,
'Illumina-A Pulldown ISC' ,
'Illumina-A Pulldown SC'  ,
'Illumina-A Pulldown WGS' ,
'Illumina-A Single ended hi seq sequencing' ,
'Illumina-A Single ended sequencing'  ,
'Illumina-B HiSeq Paired end sequencing'  ,
'Illumina-B Paired end sequencing'  ,
'Illumina-B Single ended hi seq sequencing' ,
'Illumina-B Single ended sequencing'  ,
'Illumina-C HiSeq Paired end sequencing'  ,
'Illumina-C MiSeq sequencing' ,
'Illumina-C Paired end sequencing'  ,
'Illumina-C Single ended hi seq sequencing' ,
'Illumina-C Single ended sequencing'  ,
'MiSeq sequencing'  ,
'Paired end sequencing' ,
'Single ended hi seq sequencing'  ,
'Single Ended Hi Seq Sequencing Control'  ,
'Single ended sequencing' 
)  NOT NULL; 

DROP VIEW if EXISTS `latest_seq_request`;
create view latest_seq_request as select * from seq_request where latest=true;
update schema_version set schema_version=24;