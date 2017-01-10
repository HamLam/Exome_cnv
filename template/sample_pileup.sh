#!/bin/bash

module load parallel

s_S1_R1=s_s1r1Fastq
s_S1_R2=s_s1r2Fastq

# Check to see if fastq files are compressed. If they are
# uncompress them into the working directory
#
# NOTE: The copying in the ELSE clause is not necessary. The files should be readable from data release. However, 
# there are instances where files permission are not set properly and user is unable to read files from data release. 
# This copying is a precautionary measure to make sure the program does not break if that happens. 

BWA_DB=bwa_db_value
BOWTIE2_DB=bowtie2_db_value
S_DB=seq_db

bwacommand="/home/msistaff/lamx0031/my_software/BWA/bwa-0.7.15/bwa mem -M -t 24 $BWA_DB $s_S1_R1 $s_S1_R2 | samtools view -q 10 -bS - > s_bwa_s1.bam"
btcommand="bowtie2 -p 24 -k 5 -x $BOWTIE2_DB -1 $s_S1_R1 -2 $s_S1_R2 | samtools view -q 10 -bS - > s_bowtie2_s1.bam"

echo ${bwacommand} > $WORKING_PATH/saligncommands
echo ${btcommand} >> $WORKING_PATH/saligncommands
cat ${WORKING_PATH}/saligncommands | parallel -P 2 $1

java -Xmx4g -jar  $CLASSPATH/picard.jar FixMateInformation SORT_ORDER=coordinate INPUT=s_bwa_s1.bam OUTPUT=s_bwa.fixed.bam

picard1="java -Xmx4g -jar  $CLASSPATH/picard.jar MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true METRICS_FILE=s_bwa_duplicate_stats.txt INPUT=s_bwa.fixed.bam OUTPUT=s_bwa.fixed_nodup.bam"
picard2="java -Xmx4g -jar  $CLASSPATH/picard.jar FixMateInformation SORT_ORDER=coordinate INPUT=s_bowtie2_s1.bam OUTPUT=s_bowtie2.fixed.bam"

echo ${picard1} > $WORKING_PATH/spicardcommands
echo ${picard2} >> $WORKING_PATH/spicardcommands
cat ${WORKING_PATH}/spicardcommands | parallel -P 2 $1


indexcomm1="samtools index s_bwa.fixed.bam"
indexcomm2="samtools index s_bwa.fixed_nodup.bam"
indexcomm3="samtools index s_bowtie2.fixed.bam"

echo ${indexcomm1} > $WORKING_PATH/sindexcommands
echo ${indexcomm2} >> $WORKING_PATH/sindexcommands
echo ${indexcomm3} >> $WORKING_PATH/sindexcommands
cat ${WORKING_PATH}/sindexcommands | parallel -P 3 $1


ref=/panfs/roc/rissdb/genomes/Homo_sapiens/hg19_canonical/seq/hg19_canonical.fa

samtools view -H s_bwa.fixed.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f1 |  xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d10000000 -f $ref  -r \"{}\" s_bwa.fixed.bam > \"{}\".s_bwa_count"
cat *.s_bwa_count | awk '{FS=" ";print $1,"\t",$2,"\t",$4}' - >> cnv_sample_name_bwa_pileup.txt
samtools view -H s_bwa.fixed_nodup.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f1 |  xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d10000000 -f $ref  -r \"{}\" s_bwa.fixed_nodup.bam > \"{}\".s_bwa_nodup_count"
cat *.s_bwa_nodup_count | awk '{FS=" ";print $1,"\t",$2,"\t",$4}' - >> cnv_sample_name_bwa_pileup_no_dup.txt
samtools view -H s_bowtie2.fixed.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f1 |  xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d10000000 -f $ref  -r \"{}\" s_bowtie2.fixed.bam > \"{}\".s_bowtie2_count"
cat *.s_bowtie2_count | awk '{FS=" ";print $1,"\t",$2,"\t",$4}' - >> cnv_sample_name_bowtie2_pileup.txt


#samtools mpileup -BQ0 -d10000000 -f $S_DB -q 1 s_bwa.fixed.bam | cut -f 1,2,4 > cnv_sample_name_bwa_pileup.txt
#samtools mpileup -BQ0 -d10000000 -f $S_DB -q 1 s_bwa.fixed_nodup.bam | cut -f 1,2,4 > cnv_sample_name_bwa_pileup_no_dup.txt
#samtools mpileup -BQ0 -d10000000 -f $S_DB -q 1 s_bowtie2.fixed.bam | cut -f 1,2,4 > cnv_sample_name_bowtie2_pileup.txt
