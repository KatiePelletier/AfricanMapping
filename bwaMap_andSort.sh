#!/bin/bash

# directory of processed sequences with trimmomatic 
#trimdir=/2/scratch/Katie/africans/trim/zi192ef81

#variable for output/run
set=$1

# variable for the reference genome
refGenome=/home/katie/flyGenome/dmel_r6.23/bwa/dmel-all-chromosome-r6.23.fasta 

# make output directory from mapping outputs
output=/2/scratch/Katie/africans/map/${set}

# make BWA directory path
bwa_dir=/usr/local/bwa/0.7.8

#list all files to be read (this selects the left end from each PE pair)
files=(*_R1_PE.fastq.gz)

#echo ${files[@]}

#For loop over every file
for file in ${files[@]} 
do
name=${file}

base=`basename ${name} _R1_PE.fastq.gz`

#echo ${base} 

bwa mem -t 12 -M ${refGenome} ${base}_R1_PE.fastq.gz ${base}_R2_PE.fastq.gz | samtools view -b -q 20 -@5 | samtools sort -o ${output}/${base}_bwa_PE.bam

done
