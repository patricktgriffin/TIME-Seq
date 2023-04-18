#!/usr/bin/env Rscript

# This file takes samplesheet info and user supplied pool IDs and directory
# and writes scripts for each demultiplexed TIME-Seq sample 

# read in dplyr
suppressMessages(library("dplyr"))

# get the user supplied info for scripts writing
args = commandArgs(trailingOnly=TRUE)
dir <- args[1]
samplesheet <- as.character(args[2])
sampsheet <- read.csv(paste(samplesheet), header=T)
pools <- as.list(args[3])
intFile <- as.character(args[4])
intFileName <- substr(intFile, 1,nchar(intFile)-4)
sequencing <- as.character(args[5])

# filtering the samplesheet by pools that you want to map
s <- sampsheet %>%
  dplyr::filter(pool_ID %in% pools)

# write one script for each TIME-Seq sample

for (i in unique(s$sample_ID)) {
  a <- dplyr::filter(s, sample_ID == i) #selecting just the row for each unique sampleID
  setwd(dir)
  dirslist <- list.dirs(recursive = F)
  sampDir <- dirslist[grep(a$pool_ID,dirslist)]
  setwd(sampDir)
  
  script <- file(paste(i,"_analyzeTimeSeq.sh",sep="")) #making a script
  
  
  # this is the pipeline if sequenced with PE configuration
  if (sequencing == "PE") {
    writeLines(c("#!/bin/bash", 
                 "#SBATCH -t 2-00:00",
                 "#SBATCH -p medium # Partition to run in",
                 "#SBATCH --mem 20G",
                 paste("#SBATCH -e ",dir,"/",i,"_mapping_stderr.txt",sep=''),
                 paste("#SBATCH -o ",dir,"/",i,"_mapping_stdout.txt",sep=''),
                 "",
                 "module load gcc",
                 "module load fastqc conda2/4.2.13",
                 "module list",
                 "",
                 paste("cd ",dir,sep=''),
                 paste("cd ",sampDir,sep=''),
                 paste("mkdir -p mapped_",i,sep=""),
                 #paste("if [ -f ",i,"_*R*.fastq ]; then mv ",i,"_*R*.fastq mapped_",i,"; fi",sep=""),
                 paste("mv ",i,"_*R*.fastq mapped_",i,sep=""),
                 
                 paste("cd mapped_",i,sep=""),
                 # #fastqc
                 paste("mkdir -p fastqc", sep = ""),
                 "",
                 paste("fastqc -o fastqc ",i,"*R1*fastq ",i,"*R2*fastq ", sep = ""),
                 "",
                 #load conda virtual env with updated cutadapt and bowtie2
                 "source activate py39",
                 "env",
                 "",
                  # remove adaptor sequence from reads
                  paste("cutadapt -a CTNTCTCTTATACACATCT -A CTNTCTCTTATACACATCT -o ",i,"_R1_t1.fastq -p ",i,"_R2_t1.fastq ", i,"_R1.fastq ",i,"_R2.fastq",sep=""),
                  paste("cutadapt -U 19 -o ",i,"_R1_t2.fastq -p ",i,"_R2_t2.fastq ", i,"_R1_t1.fastq ",i,"_R2_t1.fastq",sep=""),
                  # "",
                  # map with bismark
                  paste("bismark -p 4 -N 1 --unmapped --non_directional --basename ",i,"_bismark --genome /n/data2/hms/genetics/sinclair/Patrick/genomes/iGenome/",a$genome,"/iGenome/",a$genome,"/Sequence/WholeGenomeFasta -1 ",i,"_R1_t2.fastq  -2 ",i,"_R2_t2.fastq",sep = ""),
                  paste("gunzip -f ",i,"_bismark_unmapped_reads_*.fq.gz",sep = ""),
                  paste("bismark -p 4 -N 1 --non_directional --basename ",i,"_unmapped_R1 --genome /n/data2/hms/genetics/sinclair/Patrick/genomes/iGenome/",a$genome,"/iGenome/",a$genome,"/Sequence/WholeGenomeFasta -SE ",i,"_bismark_unmapped_reads_1.fq",sep = ""),
                  paste("bismark -p 4 -N 1 --non_directional --basename ",i,"_unmapped_R2 --genome /n/data2/hms/genetics/sinclair/Patrick/genomes/iGenome/",a$genome,"/iGenome/",a$genome,"/Sequence/WholeGenomeFasta -SE ",i,"_bismark_unmapped_reads_2.fq",sep = ""),

                  # filter strands that are fully methylated
                  paste("filter_non_conversion --threshold 11 ",i,"_bismark_pe.bam",sep=""),
                  paste("filter_non_conversion --threshold 11 ",i,"_unmapped_R1.bam",sep=""),
                  paste("filter_non_conversion --threshold 11 ",i,"_unmapped_R2.bam",sep=""),
                 
                  # #sort, index, convert to sam
                  paste("samtools sort ",i,"_bismark_pe.nonCG_filtered.bam > ",i,".bam",sep=""),
                  paste("samtools index ",i,".bam",sep=""),
                  paste("samtools sort ",i,"_unmapped_R1.bam > ",i,"_unmapped_R1.sorted.bam",sep=""),
                  paste("samtools index ",i,"_unmapped_R1.sorted.bam",sep=""),
                  paste("samtools sort ",i,"_unmapped_R2.bam > ",i,"_unmapped_R2.sorted.bam",sep=""),
                  paste("samtools index ",i,"_unmapped_R2.sorted.bam",sep=""),                 

                  # extract methylation levels from deduplicated file
                  paste("bismark_methylation_extractor --parallel 4 --ignore 10 --ignore_r2 10 --ignore_3prime 10 --ignore_3prime_r2 10 --no_overlap --genome_folder /n/data2/hms/genetics/sinclair/Patrick/genomes/iGenome/",a$genome,"/iGenome/",a$genome,"/Sequence/WholeGenomeFasta/ ",i,"_bismark_pe.nonCG_filtered.bam",sep=""),
                  paste("bismark_methylation_extractor  --parallel 4 --ignore 10 --ignore_3prime 10 --genome_folder /n/data2/hms/genetics/sinclair/Patrick/genomes/iGenome/",a$genome,"/iGenome/",a$genome,"/Sequence/WholeGenomeFasta/ ",i,"_unmapped_R1.nonCG_filtered.bam",sep=""),
                  paste("bismark_methylation_extractor  --parallel 4 --ignore 10 --ignore_3prime 10 --genome_folder /n/data2/hms/genetics/sinclair/Patrick/genomes/iGenome/",a$genome,"/iGenome/",a$genome,"/Sequence/WholeGenomeFasta/ ",i,"_unmapped_R2.nonCG_filtered.bam",sep=""),
                  paste("bismark2bedGraph -o ",i,"_bismark_merged CpG_*_",i,"_bismark_pe.nonCG_filtered.txt CpG_*_",i,"_unmapped_R1.nonCG_filtered.txt CpG_*_",i,"_unmapped_R2.nonCG_filtered.txt",sep=""),
                  
                # cleanup
                 paste("cp ",dir,"/",i,"_mapping_stderr.txt ",i,"_mapping_stderr.txt", sep = ""),
                 paste("cp ",dir,"/",i,"_mapping_stdout.txt ",i,"_mapping_stdout.txt", sep = ""),
                 paste("rm ",dir,"/",i,"_mapping_stderr.txt", sep = ""),
                 paste("rm ",dir,"/",i,"_mapping_stdout.txt", sep = ""),
                 paste("mv ../",i,"_analyzeTimeSeq.sh ./",sep="")), script)
  }
  
# this is the pipeline if sequenced with SE configuration

  if (sequencing == "SE") {
    writeLines(c("#!/bin/bash", 
                 "#SBATCH -t 0-02:00",
                 "#SBATCH -p short # Partition to run in",
                 "#SBATCH --mem 20G",
                 "#SBATCH -c 4",
                 paste("#SBATCH -e ",dir,"/",i,"_mapping_stderr.txt",sep=''),
                 paste("#SBATCH -o ",dir,"/",i,"_mapping_stdout.txt",sep=''),
                 
                 #load softwares and make directories for the sample
                 "module load gcc",
                 "module load fastqc conda2/4.2.13",
                 paste("cd ",dir,sep=''),
                 paste("cd ",sampDir,sep=''),
                 paste("mkdir -p mapped_",i,sep=""),
                 paste("mv ",i,"_*R*.fastq mapped_",i,sep=""),
                 paste("cd mapped_",i,sep=""),
                 
                 #run fastqc for each sample
                 paste("mkdir fastqc", sep = ""),
                 paste("fastqc -o fastqc ",i,"*R1*fastq ", sep = ""),
                 
                 #load conda virtual env with updated cutadapt and bowtie2
                 "source activate py39",
                 "env",
                 
                 #Map and get methylation data with bismark
                 paste("cutadapt --minimum-length 40 -a CTNTCTCTTATACACATCT -o ",i,"_R1_t1.fastq ", i,"_R1.fastq",sep=""),
                 paste("bismark -p 4 -N 1 --non_directional --basename ",i," --genome /n/data2/hms/genetics/sinclair/Patrick/genomes/iGenome/",a$genome,"/iGenome/",a$genome,"/Sequence/WholeGenomeFasta ",i,"_R1_t1.fastq ",sep = ""),
                 paste("filter_non_conversion --threshold 11 ",i,".bam",sep=""),
                 
                 #sort, index, sam file
                 paste("samtools sort ",i,".nonCG_filtered.bam > ",i,".sorted.bam",sep=""),
                 paste("samtools index ",i,".sorted.bam",sep=""),
                 paste("samtools view -Sh ",i,".nonCG_filtered.bam > ",i,".sam", sep = ""),
                 
                 # #extract methylation levels
                 paste("bismark_methylation_extractor  --parallel 4 --ignore 10 --ignore_3prime 10 --bedGraph --genome_folder /n/data2/hms/genetics/sinclair/Patrick/genomes/iGenome/",a$genome,"/iGenome/",a$genome,"/Sequence/WholeGenomeFasta/ ",i,".nonCG_filtered.bam",sep=""),
                 
                 # bismark reports 
                 paste("bismark2report ",i,".nonCG_filtered.bam",sep=""),
                 paste("bismark2summary ",i,".nonCG_filtered.bam",sep=""),
                 paste("module load python/2.7.12",sep=""),
                 paste("source ~/python2.7/bin/activate", sep="")), script)
    
  }
  

  close(script)
}
