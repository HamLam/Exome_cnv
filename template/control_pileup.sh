#!/bin/bash

module load parallel

c_S1_R1=c_s1r1Fastq
c_S1_R2=c_s1r2Fastq

# Check to see if fastq files are compressed. If they are
# uncompress them into the working directory
#
# NOTE: The copying in the ELSE clause is not necessary. The files should be readable from data release. However, 
# there are instances where files permission are not set properly and user is unable to read files from data release. 
# This copying is a precautionary measure to make sure the program does not break if that happens. 

WORKING_PATH=working_dir
BWA_DB=bwa_db_value
BOWTIE2_DB=bowtie2_db_value
S_DB=seq_db

bwacommand="/home/msistaff/lamx0031/my_software/BWA/bwa-0.7.15/bwa mem -M -t 24 $BWA_DB $c_S1_R1 $c_S1_R2 | samtools view -q 10 -bS - > c_bwa.bam"
btcommand="bowtie2 -p 24 -k 5 -x $BOWTIE2_DB -1 $c_S1_R1 -2 $c_S1_R2 | samtools view -q 10 -bS - > c_bowtie2.bam"

echo ${bwacommand} > $WORKING_PATH/aligncommands
echo ${btcommand} >> $WORKING_PATH/aligncommands
cat ${WORKING_PATH}/aligncommands | parallel -P 2 $1

java -Xmx4g -jar  $CLASSPATH/picard.jar FixMateInformation SORT_ORDER=coordinate INPUT=c_bwa.bam OUTPUT=c_bwa.fixed.bam
picard1="java -Xmx4g -jar  $CLASSPATH/picard.jar MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true METRICS_FILE=c_bwa_duplicate_stats.txt INPUT=c_bwa.fixed.bam OUTPUT=c_bwa.fixed_nodup.bam"
picard2="java -Xmx4g -jar  $CLASSPATH/picard.jar FixMateInformation SORT_ORDER=coordinate INPUT=c_bowtie2.bam OUTPUT=c_bowtie2.fixed.bam"
echo ${picard1} > $WORKING_PATH/cpicardcommands
echo ${picard2} >> $WORKING_PATH/cpicardcommands
cat ${WORKING_PATH}/cpicardcommands | parallel -P 2 $1

indexcomm1="samtools index c_bwa.fixed.bam"
indexcomm2="samtools index c_bwa.fixed_nodup.bam"
indexcomm3="samtools index c_bowtie2.fixed.bam"
echo ${indexcomm1} > $WORKING_PATH/indexcommands
echo ${indexcomm2} >> $WORKING_PATH/indexcommands
echo ${indexcomm3} >> $WORKING_PATH/indexcommands
cat ${WORKING_PATH}/indexcommands | parallel -P 3 $1

ref=/panfs/roc/rissdb/genomes/Homo_sapiens/hg19_canonical/seq/hg19_canonical.fa

samtools view -H c_bwa.fixed.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f1 |  xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d10000000 -f $ref  -r \"{}\" c_bwa.fixed.bam > \"{}\".c_bwa_count"
cat *.c_bwa.count | awk '{FS=" ";print $1,"\t",$2,"\t",$4}' - >> cnv_control_name_bwa_pileup.txt
samtools view -H c_bwa.fixed_nodup.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f1 |  xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d10000000 -f $ref  -r \"{}\" c_bwa.fixed_nodup.bam > \"{}\".c_bwa_nodup_count"
cat *.c_bwa_nodup_count | awk '{FS=" ";print $1,"\t",$2,"\t",$4}' - >> cnv_control_name_bwa_pileup_no_dup.txt
samtools view -H c_bowtie2.fixed.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f1 |  xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d10000000 -f $ref  -r \"{}\" c_bowtie2.fixed.bam > \"{}\".c_bowtie2_count"
cat *.c_bowtie2_count | awk '{FS=" ";print $1,"\t",$2,"\t",$4}' - >> cnv_control_name_bowtie2_pileup.txt

