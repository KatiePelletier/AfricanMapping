#!/bin/bash

#raw_dir=/home/katie/africans/zi192ef96
trim=/usr/local/trimmomatic/trimmomatic-0.36.jar
adapt=/home/katie/adapterSequences/TruSeq3-PE.fa:2:30:10
#trim_dir=/2/scratch/Katie/africans/trim
trim_dir=$1

files=(*_R1.fastq.gz)
for file in ${files[@]} 
do
name=${file}
base=`basename ${name} _R1.fastq.gz`
java -jar ${trim} PE -threads 5 -phred33 -trimlog ${trim_dir}/trimlog.txt ${base}_R1.fastq.gz ${base}_R2.fastq.gz ${trim_dir}/${base}_R1_PE.fastq.gz ${trim_dir}/${base}_R1_SE.fastq.gz ${trim_dir}/${base}_R2_PE.fastq.gz ${trim_dir}/${base}_R2_SE.fastq.gz ILLUMINACLIP:${adapt} LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36
done
