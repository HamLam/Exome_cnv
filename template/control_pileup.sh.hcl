#!/bin/bash

c_S1_R1=c_s1r1Fastq
c_S1_R2=c_s1r2Fastq
c_S2_R1=c_s2r1Fastq
c_S2_R2=c_s2r2Fastq

# Check to see if fastq files are compressed. If they are
# uncompress them into the working directory
#
# NOTE: The copying in the ELSE clause is not necessary. The files should be readable from data release. However, 
# there are instances where files permission are not set properly and user is unable to read files from data release. 
# This copying is a precautionary measure to make sure the program does not break if that happens. 

if [[ $c_S1_R1 = *.gz ]] ; then
    gunzip -c $c_S1_R1 > c_S1_R1.fastq
    c_S1_R1=c_S1_R1.fastq
else
    cp $c_S1_R1 c_S1_R1.fastq
    c_S1_R1=c_S1_R1.fastq
fi	

if [[ $c_S1_R2 = *.gz ]] ; then
    gunzip -c $c_S1_R2 > c_S1_R2.fastq
    c_S1_R2=c_S1_R2.fastq
else 
    cp $c_S1_R2 c_S1_R2.fastq
    c_S1_R2=c_S1_R2.fastq
fi

if [[ $c_S2_R1 = *.gz ]] ; then
    gunzip -c $c_S2_R1 > c_S2_R1.fastq
    c_S2_R1=c_S2_R1.fastq
else
    cp $c_S2_R1 c_S2_R1.fastq
    c_S2_R1=c_S2_R1.fastq
fi	

if [[ $c_S2_R2 = *.gz ]] ; then
    gunzip -c $c_S2_R2 > c_S2_R2.fastq
    c_S2_R2=c_S2_R2.fastq
else
    cp $c_S2_R2 c_S2_R2.fastq
    c_S2_R2=c_S2_R2.fastq
fi   

BWA_DB=bwa_db_value
BOWTIE2_DB=bowtie2_db_value
S_DB=seq_db

# bwa mem -M -t 24 $BWA_DB $c_S1_R1 $c_S1_R2 > c_bwa_s1.sam
# bwa mem -M -t 24 $BWA_DB $c_S2_R1 $c_S2_R2 > c_bwa_s2.sam

bwa mem -M -t 24 $BWA_DB $c_S1_R1 $c_S1_R2 | samtools view -q 10 -bS - > c_bwa_s1.bam
bwa mem -M -t 24 $BWA_DB $c_S2_R1 $c_S2_R2 | samtools view -q 10 -bS - > c_bwa_s2.bam

# bowtie2 -p 24 -k 5 -x $BOWTIE2_DB -1 $c_S1_R1 -2 $c_S1_R2 -S c_bowtie2_s1.sam
# bowtie2 -p 24 -k 5 -x $BOWTIE2_DB -1 $c_S2_R1 -2 $c_S2_R2 -S c_bowtie2_s2.sam

bowtie2 -p 24 -k 5 -x $BOWTIE2_DB -1 $c_S1_R1 -2 $c_S1_R2 | samtools view -q 10 -bS - > c_bowtie2_s1.bam
bowtie2 -p 24 -k 5 -x $BOWTIE2_DB -1 $c_S2_R1 -2 $c_S2_R2 | samtools view -q 10 -bS - > c_bowtie2_s2.bam

# samtools view -q 10 -bS c_bwa_s1.sam > c_bwa_s1.bam
# samtools view -q 10 -bS c_bwa_s2.sam > c_bwa_s2.bam
# samtools view -q 10 -bS c_bowtie2_s1.sam > c_bowtie2_s1.bam
# samtools view -q 10 -bS c_bowtie2_s2.sam > c_bowtie2_s2.bam

samtools merge c_bwa.bam c_bwa_s1.bam c_bwa_s2.bam
samtools merge c_bowtie2.bam c_bowtie2_s1.bam c_bowtie2_s2.bam

samtools sort c_bwa.bam c_bwa.sort

## suggested by Christy to replace FixMateInformation
## $PTOOL/picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT

java -Xmx4g -jar  $CLASSPATH/picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true ASSUME_SORTED=true METRICS_FILE=c_bwa_duplicate_stats.txt INPUT=c_bwa.sort.bam OUTPUT=c_bwa.fixed.bam


# bedtools bamtobed -i c_bwa.fixed.bam > c_bwa.bed
# bedtools coverage -sorted -a c_bwa.bed -b c_bwa.fixed.bam > c_bwa_coverage # get back same number of rows as in c_bwa.bed


samtools mpileup -BQ0 -d10000000 -f $S_DB -q 1 c_bwa.fixed.bam | cut -f 1,2,4 > cnv_control_name_bwa_pileup.txt

#samtools mpileup -BQ0 -d10000000 -f $S_DB -q 1 c_bwa.fixed_nodup.bam | cut -f 1,2,4 > cnv_control_name_bwa_pileup_no_dup.txt
#samtools mpileup -BQ0 -d10000000 -f $S_DB -q 1 c_bowtie2.fixed.bam | cut -f 1,2,4 > cnv_control_name_bowtie2_pileup.txt
