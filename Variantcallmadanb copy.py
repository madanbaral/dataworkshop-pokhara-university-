#! /usr/bin/env python

import sys
import os

def main():
    
    # Preprocessing and QC with FastQC V-11.7
    os.system("fastqc -t 4 -o qc/raw *.fastq.gz")
    
    #Trimming reads to remove bad reads and improve quality of reads with Trimmomatic V-0.36
    os.system("sh trimmomatic.sh")

    #Rerunning QC to evaluate the quality of reads with FastQC V-11.7
    os.system("fastqc -t 4 -o qc/qctrimmed qc/trimmomatic/trimmed/*trimmed.fastq.gz")

    #Aligning the Parent reads with HISAT2 V-2.1.0 to reference genome

    #Creating HISAT2 V-2.1.0  Index before Aligning in three steps
    os.system("hisat2_extract_splice_sites.py /Users/madanbaral/GCF_000005845.2_ASM584v2_genomic.gff >ecoliparent.ss")

    os.system("hisat2_extract_exons.py /Users/madanbaral/parent.gtf >ecoliparent.exon")

    os.system("hisat2-build --ss ecoliparent.ss --exon ecoliparent.exon /Users/madanbaral/GCF_000005845.2_ASM584v2_genomic.fna ecoliparent_tran")

    #Aligning RNAseq Parent sample reads to Reference Genome
    os.system("hisat2 -p 8 --dta -x /Users/madanbaral/ecoliparent_tran -1 /Users/madanbaral/qc/trimmomatic/trimmed/sample-1_S1_L001_R1_001_trimmed.fastq.gz -2 /Users/madanbaral/qc/trimmomatic/trimmed/sample-1_S1_L001_R2_001_trimmed.fastq.gz -S parent.sam")

    #Sorting and converting SAM to BAM file usung SAMTOOLS V-1.7
    os.system("samtools sort -@ 8 -o parent.bam parent.sam")

    #Assembling and Quantifying expressed genes and transcripts usung STRINGTIE V-1.3.4c
    os.system("stringtie -p 8 -G /Users/madanbaral/GCF_000005845.2_ASM584v2_genomic.gff -o parent.gtf parent.bam")

    #Extracting consensus transcript sequences using BEDTOOLS V-2.27.1
    os.system("bedtools getfasta -fi /Users/madanbaral/GCF_000005845.2_ASM584v2_genomic.fasta -bed parent.gtf -fo parent.fasta")

    #Aligning the Daughter reads with HISAT2 V-2.1.0 to parent consensus genome
    
    #Creating HISAT2 V-2.1.0 Index before Aligning in three steps
    os.system("hisat2_extract_splice_sites.py /Users/madanbaral/parent.gtf >ecolidaughter.ss")
    
    os.system("hisat2_extract_exons.py /Users/madanbaral/parent.gtf >ecolidaughter.exon")
    
    os.system("hisat2-build --ss ecolidaughter.ss --exon ecolidaughter.exon /Users/madanbaral/parent.fasta ecolidaughter_tran")
    
    #Aligning RNAseq Parent sample reads to Reference Genome with HISAT2 V-2.1.0
    os.system("hisat2 -p 8 --dta -x /Users/madanbaral/ecolidaughter_tran -1 /Users/madanbaral/qc/trimmomatic/trimmed/sample-2_S6_L001_R1_001_trimmed.fastq.gz -2 /Users/madanbaral/qc/trimmomatic/trimmed/sample-2_S6_L001_R2_001_trimmed.fastq.gz -S daughter.sam")
    
    #Sorting and converting Daughter SAM to BAM file using SAMTOOLS V-1.7
    os.system("samtools sort -@ 8 -o daughter.bam daughter.sam")
    
    #Creating transcript sequences from the reference genome files using BEDTOOLS V-2.27.1
    os.system("bedtools getfasta -fi /Users/madanbaral/GCF_000005845.2_ASM584v2_genomic.fasta -bed GCF_000005845.2_ASM584v2_genomic.gff -fo ref.fasta")
    
    #Creating VCF files from daughter BAM file using SAMTOOLS V-1.7 and BCFTOOLS V-1.7
    os.system("samtools mpileup -E -uf /Users/madanbaral/ref.fasta daughter.bam >daughter.mpileup")
    os.system("bcftools view -vcg daughter.mpileup >ann.vcf")
    os.system("vcfutils.pl varFilter -D100 > ann.vcf")

if __name__ == "__main__":
    main()

