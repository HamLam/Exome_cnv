/* DROP TABLE IF EXISTS `cnv_control_name_bwa_pileup_no_dup`;*/
CREATE TABLE `cnv_control_name_bwa_pileup_no_dup`(
chr VARCHAR(32),
pos INT,
coverage INT,
mapped_reads INT
);
/* LOAD DATA LOCAL INFILE 'cnv_control_name_bwa_pileup_no_dup.txt' INTO TABLE `cnv_control_name_bwa_pileup_no_dup` FIELDS TERMINATED BY '\t';
CREATE INDEX `cnv_control_name_pileup_nd` ON `cnv_control_name_bwa_pileup_no_dup`(chr,pos);*/

/* DROP TABLE IF EXISTS `cnv_control_name_bwa_pileup`; */
CREATE TABLE `cnv_control_name_bwa_pileup`(
chr VARCHAR(32),
pos INT,
coverage INT
);
/* LOAD DATA LOCAL INFILE 'cnv_control_name_bwa_pileup.txt' INTO TABLE `cnv_control_name_bwa_pileup` FIELDS TERMINATED BY '\t';
CREATE INDEX `cnv_control_name_bwa_pileup_i1` ON `cnv_control_name_bwa_pileup`(chr,pos);*/

/* DROP TABLE IF EXISTS `cnv_control_name_bowtie_pileup`;*/
CREATE TABLE `cnv_control_name_bowtie_pileup`(
chr VARCHAR(32),
pos INT,
coverage INT
);
/* LOAD DATA LOCAL INFILE 'cnv_control_name_bowtie2_pileup.txt' INTO TABLE `cnv_control_name_bowtie_pileup` FIELDS TERMINATED BY '\t';
CREATE INDEX `cnv_control_name_bowtie_pileup_i1` ON `cnv_control_name_bowtie_pileup`(chr,pos);*/

