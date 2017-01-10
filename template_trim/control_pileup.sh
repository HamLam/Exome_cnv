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

chrfiles_path=exo_chrfiles
script_path=scripts_location
WORKING_PATH=working_dir
BWA_DB=bwa_db_value
BOWTIE2_DB=bowtie2_db_value
S_DB=seq_db
ref=/panfs/roc/rissdb/genomes/Homo_sapiens/hg19_canonical/seq/hg19_canonical.fa

bwacommand="/home/msistaff/lamx0031/my_software/BWA/bwa-0.7.15/bwa mem -M -t 24 $BWA_DB $c_S1_R1 $c_S1_R2 | samtools view -q 10 -bS - > c_bwa.bam"
btcommand="bowtie2 -p 24 -k 5 -x $BOWTIE2_DB -1 $c_S1_R1 -2 $c_S1_R2 | samtools view -q 10 -bS - > c_bowtie2.bam"

echo ${bwacommand} > $WORKING_PATH/aligncommands
echo ${btcommand} >> $WORKING_PATH/aligncommands
cat ${WORKING_PATH}/aligncommands | parallel -j +0 $1

java -Xmx4g -jar  $CLASSPATH/picard.jar FixMateInformation SORT_ORDER=coordinate INPUT=c_bwa.bam OUTPUT=c_bwa.fixed.bam
picard1="java -Xmx4g -jar  $CLASSPATH/picard.jar MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true METRICS_FILE=c_bwa_duplicate_stats.txt INPUT=c_bwa.fixed.bam OUTPUT=c_bwa.fixed_nodup.bam"
picard2="java -Xmx4g -jar  $CLASSPATH/picard.jar FixMateInformation SORT_ORDER=coordinate INPUT=c_bowtie2.bam OUTPUT=c_bowtie2.fixed.bam"
echo ${picard1} > $WORKING_PATH/cpicardcommands
echo ${picard2} >> $WORKING_PATH/cpicardcommands
cat ${WORKING_PATH}/cpicardcommands | parallel -j +0 $1

indexcomm1="samtools index c_bwa.fixed.bam"
indexcomm2="samtools index c_bwa.fixed_nodup.bam"
indexcomm3="samtools index c_bowtie2.fixed.bam"
echo ${indexcomm1} > $WORKING_PATH/indexcommands
echo ${indexcomm2} >> $WORKING_PATH/indexcommands
echo ${indexcomm3} >> $WORKING_PATH/indexcommands
cat ${WORKING_PATH}/indexcommands | parallel -j +0 $1


samtools view -H c_bwa.fixed.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f1 |  xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d10000000 -f $ref  -r \"{}\" c_bwa.fixed.bam | cut -f 1,2,4 > cnv_control_name_bwa_pileup.\"{}\""

chr1_files=($WORKING_PATH/cnv_control_name_bwa_pileup.chr*)
chr2_files=($chrfiles_path/file.*)
for ((i=0;i<${#chr1_files[@]};i++)); do
echo "perl $script_path/trimfile_t.pl "${chr1_files[i]}" "${chr2_files[i]}" " >> chopcommands
done
cat ${WORKING_PATH}/chopcommands | parallel -j +0

#cat *.c_bwa.count | awk '{FS=" ";print $1,"\t",$2,"\t",$4}' - >> cnv_control_name_bwa_pileup.txt

samtools view -H c_bwa.fixed_nodup.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f1 |  xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d10000000 -f $ref  -r \"{}\" c_bwa.fixed_nodup.bam | cut -f 1,2,4 > cnv_control_name_bwa_pileup_no_dup.\"{}\""
chr1a_files=($WORKING_PATH/cnv_control_name_bwa_pileup_no_dup.chr*)
chr2a_files=($chrfiles_path/file.*)
for ((i=0;i<${#chr1a_files[@]};i++)); do
echo "perl $script_path/trimfile_t.pl "${chr1a_files[i]}" "${chr2a_files[i]}" " >> chopcommands2
done
cat ${WORKING_PATH}/chopcommands2 | parallel -j +0

#cat *.c_bwa_nodup_count | awk '{FS=" ";print $1,"\t",$2,"\t",$4}' - >> cnv_control_name_bwa_pileup_no_dup.txt

samtools view -H c_bowtie2.fixed.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f1 |  xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -BQ0 -d10000000 -f $ref  -r \"{}\" c_bowtie2.fixed.bam | cut -f 1,2,4 > cnv_control_name_bowtie_pileup.\"{}\""
chr1b_files=($WORKING_PATH/cnv_control_name_bowtie_pileup.chr*)
chr2b_files=($chrfiles_path/file.*)
for ((i=0;i<${#chr1b_files[@]};i++)); do
echo "perl $script_path/trimfile_t.pl "${chr1b_files[i]}" "${chr2b_files[i]}" " >> chopcommands3
done
cat ${WORKING_PATH}/chopcommands3 | parallel -j +0



