#!/bin/bash

# Threads
THREADS=8

# Create folders
mkdir -p data qc alignment variants expression epigenomics motifs reference

#QC 
wget -P data \
https://www.encodeproject.org/files/ENCFF002DKA/@@download/ENCFF002DKA.fastq.gz \
https://www.encodeproject.org/files/ENCFF002DKE/@@download/ENCFF002DKE.fastq.gz

fastqc data/*.fastq.gz -o qc -t $THREADS
multiqc qc -o qc

#Alignment 
wget -P reference https://genome-idx.s3.amazonaws.com/bt/hg19.zip
unzip reference/hg19.zip -d reference/

bowtie2 -p $THREADS -x reference/hg19 \
-1 data/ENCFF002DKA.fastq.gz \
-2 data/ENCFF002DKE.fastq.gz \
-S alignment/aligned.sam

samtools view -bS alignment/aligned.sam > alignment/aligned.bam
samtools sort alignment/aligned.bam -o alignment/aligned_sorted.bam
samtools index alignment/aligned_sorted.bam

#Variant Calling 
wget -P reference https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip -k reference/hg19.fa.gz

bcftools mpileup -f reference/hg19.fa alignment/aligned_sorted.bam > variants/raw.vcf
bcftools call -mv variants/raw.vcf > variants/final.vcf

#Expression (Cufflinks) 
cufflinks -G hg19.gtf alignment/aligned_sorted.bam -o expression/

cuffdiff -o expression/diff hg19.gtf \
alignment/aligned_sorted.bam alignment/aligned_sorted.bam

#ChIP-seq 
bowtie2 -x reference/hg19 \
-U chip.fastq.gz -S epigenomics/chip.sam

bowtie2 -x reference/hg19 \
-U control.fastq.gz -S epigenomics/control.sam

samtools view -bS epigenomics/chip.sam > epigenomics/chip.bam
samtools sort epigenomics/chip.bam -o epigenomics/chip_sorted.bam

samtools view -bS epigenomics/control.sam > epigenomics/control.bam
samtools sort epigenomics/control.bam -o epigenomics/control_sorted.bam

macs2 callpeak \
-t epigenomics/chip_sorted.bam \
-c epigenomics/control_sorted.bam \
-f BAM -g hs -n peaks --outdir epigenomics/

#Motif 
bedtools getfasta -fi reference/hg19.fa \
-bed epigenomics/peaks_peaks.narrowPeak \
-fo motifs/sequences.fasta

meme motifs/sequences.fasta -dna -oc motifs/ -nmotifs 5
