#!/usr/bin/env Rscript

#load dplyr R library
library("dplyr")

# collecting input from analyzeTimeSeq.sh script.
args = commandArgs(trailingOnly=TRUE)
dir <- args[1] 
samplesheet <- as.character(args[2]) 
sampsheet <- read.csv(paste(samplesheet), header=T)
pools <- as.list(args[3])
intFile <- as.character(args[4])
intFileName <- substr(intFile, 1,nchar(intFile)-4)
sequencing <- as.character(args[5])
genDir <- as.character(args[6])

# filtering samplesheet by selected pools
s <- sampsheet %>%
  dplyr::filter(pool_ID %in% pools)

# writing scripts near 
for (i in unique(s$sample_ID)) {
  a <- dplyr::filter(s, sample_ID == i) #selecting just the row for each unique sampleID
  setwd(dir)
  dirslist <- list.dirs(recursive = F)
  sampDir <- dirslist[grep(a$pool_ID,dirslist)]
  setwd(sampDir)
  
  script <- file(paste(i,"_analyzeTimeSeq.sh",sep="")) #making a script
  
  
  if (sequencing == "PE") {
    writeLines(c("#!/bin/bash", 
                 "#SBATCH -t 2-00:00", # time to run pipeline, for PE 1m reads this will overshoot
                 "#SBATCH -p medium", # partition to run in
                 "#SBATCH --mem 20G", # need at least 15GB for Bismark mapping
                 paste("#SBATCH -e ",dir,"/",i,"_mapping_stderr.txt",sep=''),
                 paste("#SBATCH -o ",dir,"/",i,"_mapping_stdout.txt",sep=''),
                 "",
                 #load required modules
                 "module load gcc",
                 "module load fastqc conda2/4.2.13",
                 
                 #change and make 
                 paste("cd ",dir,sep=''),
                 paste("cd ",sampDir,sep=''),
                 paste("mkdir -p mapped_",i,sep=""),
                 paste("mv ",i,"_*R*.fastq mapped_",i,sep=""),
                 paste("cd mapped_",i,sep=""),
                 
                 #fastqc
                 paste("mkdir -p fastqc", sep = ""),
                 paste("fastqc -o fastqc ",i,"*R1*fastq ",i,"*R2*fastq ", sep = ""),
                 
                 #load conda virtual env with updated cutadapt and bowtie2
                 "source activate py39",
                 "env",
                  
                  # remove adaptor sequence from reads
                  paste("cutadapt -a CTNTCTCTTATACACATCT -A CTNTCTCTTATACACATCT -o ",i,"_R1_t1.fastq -p ",i,"_R2_t1.fastq ", i,"_R1.fastq ",i,"_R2.fastq",sep=""),
                  paste("cutadapt -U 19 -o ",i,"_R1_t2.fastq -p ",i,"_R2_t2.fastq ", i,"_R1_t1.fastq ",i,"_R2_t1.fastq",sep=""),
                  
                  # map with bismark
                  paste("bismark -p 4 -N 1 --unmapped --non_directional --basename ",i,"_bismark --genome ",genDir," -1 ",i,"_R1_t2.fastq  -2 ",i,"_R2_t2.fastq",sep = ""),
                  paste("gunzip -f ",i,"_bismark_unmapped_reads_*.fq.gz",sep = ""),
                  paste("bismark -p 4 -N 1 --non_directional --basename ",i,"_unmapped_R1 --genome ",genDir," -SE ",i,"_bismark_unmapped_reads_1.fq",sep = ""),
                  paste("bismark -p 4 -N 1 --non_directional --basename ",i,"_unmapped_R2 --genome ",genDir," -SE ",i,"_bismark_unmapped_reads_2.fq",sep = ""),
                  
                  # filter strands that are fully methylated (artificially)
                  paste("filter_non_conversion --threshold 11 ",i,"_bismark_pe.bam",sep=""),
                  paste("filter_non_conversion --threshold 11 ",i,"_unmapped_R1.bam",sep=""),
                  paste("filter_non_conversion --threshold 11 ",i,"_unmapped_R2.bam",sep=""),
                  
                  # sort, index
                  paste("samtools sort ",i,"_bismark_pe.nonCG_filtered.bam > ",i,".bam",sep=""),
                  paste("samtools index ",i,".bam",sep=""),
                  paste("samtools sort ",i,"_unmapped_R1.bam > ",i,"_unmapped_R1.sorted.bam",sep=""),
                  paste("samtools index ",i,"_unmapped_R1.sorted.bam",sep=""),
                  paste("samtools sort ",i,"_unmapped_R2.bam > ",i,"_unmapped_R2.sorted.bam",sep=""),
                  paste("samtools index ",i,"_unmapped_R2.sorted.bam",sep=""),                 
                 
                  # extract DNAme levels 
                  paste("bismark_methylation_extractor --parallel 4 --ignore 10 --ignore_r2 10 --ignore_3prime 10 --ignore_3prime_r2 10 --no_overlap --genome_folder ",genDir," ",i,"_bismark_pe.nonCG_filtered.bam",sep=""),
                  paste("bismark_methylation_extractor  --parallel 4 --ignore 10 --ignore_3prime 10 --genome_folder ",genDir," ",i,"_unmapped_R1.nonCG_filtered.bam",sep=""),
                  paste("bismark_methylation_extractor  --parallel 4 --ignore 10 --ignore_3prime 10 --genome_folder ",genDir," ",i,"_unmapped_R2.nonCG_filtered.bam",sep=""),
                  
                  # combining the DNAme data from pair-mapped reads and the individually-mapped R1 and R2 
                  paste("bismark2bedGraph -o ",i,"_bismark_merged CpG_*_",i,"_bismark_pe.nonCG_filtered.txt CpG_*_",i,"_unmapped_R1.nonCG_filtered.txt CpG_*_",i,"_unmapped_R2.nonCG_filtered.txt",sep=""),
                  
                  ## collecting simple stats for each sample / pool
                  # mapping
                  paste("PAIRS=$(grep 'Sequence pairs analysed in total:' ",i,"_bismark_PE_report.txt | cut -f 2)", sep=""),
                  paste("MAPEFF=$(grep 'Mapping efficiency:' ",i,"_bismark_PE_report.txt | cut -f 2 | tr -d '%')", sep=""),
                  paste("MAPEFFR1=$(grep 'Mapping efficiency:' ",i,"_unmapped_R1_SE_report.txt | cut -f 2 | tr -d '%')", sep=""),
                  paste("MAPEFFR2=$(grep 'Mapping efficiency:' ",i,"_unmapped_R2_SE_report.txt | cut -f 2 | tr -d '%')", sep=""),
                  
                  # reads removed due to being fully methylated (by the non-con function)
                  paste("NONCON=$(grep 'Sequences removed' ",i,"_bismark_pe.non-conversion_filtering.txt | cut -d ' ' -f 18 | tr -d '(' | tr -d ')' |tr -d '%' )", sep=""),
                  paste("NONCONR1=$(grep 'Sequences removed' ",i,"_unmapped_R1.non-conversion_filtering.txt | cut -d ' ' -f 15 | tr -d '(' | tr -d ')' |tr -d '%' )", sep=""),
                  paste("NONCONR2=$(grep 'Sequences removed' ",i,"_unmapped_R2.non-conversion_filtering.txt | cut -d ' ' -f 15 | tr -d '(' | tr -d ')' |tr -d '%' )", sep=""),
                  
                  # methylation percent
                  paste("MCPG=$(grep 'C methylated in CpG context:' ",i,"_bismark_pe.nonCG_filtered_splitting_report.txt | cut -f 2 | tr -d '%')", sep=""),
                  paste("MCHH=$(grep 'C methylated in CHH context:' ",i,"_bismark_pe.nonCG_filtered_splitting_report.txt | cut -f 2 | tr -d '%')", sep=""),
                  # ontarget analysis
                  paste("MAPPED=$(samtools view -b -F 4 ",i,"_bismark_pe.bam | samtools view -c)",sep=""),
                  paste("MAPPEDR1=$(samtools view -b -F 4 ",i,"_unmapped_R1.sorted.bam | samtools view -c)",sep=""),
                  paste("MAPPEDR2=$(samtools view -b -F 4 ",i,"_unmapped_R2.sorted.bam | samtools view -c)",sep=""),
                  
                  paste("ONTARG=$(bedtools intersect -u -a ",i,"_bismark_pe.bam -b ",dir,"/captureAnalysis/",intFile," | samtools view -c)", sep=""), 
                  paste("ONTARGR1=$(bedtools intersect -u -a ",i,"_unmapped_R1.bam -b ",dir,"/methylationClock/captureAnalysis/",intFile," | samtools view -c)", sep=""), 
                  paste("ONTARGR2=$(bedtools intersect -u -a ",i,"_unmapped_R2.bam -b ",dir,"/captureAnalysis/",intFile," | samtools view -c)", sep=""),
                  
                  'PRCNT=$(echo "scale = 5; ${ONTARG} / ${MAPPED} * 100" | bc -l)',
                  'PRCNTR1=$(echo "scale = 5; ${ONTARGR1} / ${MAPPEDR1} * 100" | bc -l)',
                  'PRCNTR2=$(echo "scale = 5; ${ONTARGR2} / ${MAPPEDR2} * 100" | bc -l)',
                  
                  #printing to stats output file
                  paste("echo ",a$pool_ID,",",i,",$PAIRS,$MAPEFF,$MAPEFFR1,$MAPEFFR2,$NONCON,$NONCONR1,$NONCONR2,$MCPG,$MCHH,$MAPPED,$MAPPEDR1,$MAPPEDR2,$ONTARG,$ONTARGR1,$ONTARGR2,$PRCNT,$PRCNTR1,$PRCNTR2 >> ",dir,"/runStats_",a$pool_ID,"_",intFileName,".bed.csv", sep=""),
                  
                  # cleanup
                  paste("cp ",dir,"/",i,"_mapping_stderr.txt ",i,"_mapping_stderr.txt", sep = ""),
                  paste("cp ",dir,"/",i,"_mapping_stdout.txt ",i,"_mapping_stdout.txt", sep = ""),
                  paste("rm ",dir,"/",i,"_mapping_stderr.txt", sep = ""),
                  paste("rm ",dir,"/",i,"_mapping_stdout.txt", sep = ""),
                  paste("mv ../",i,"_analyzeTimeSeq.sh ./",sep="")), script)
    }
  
  
  if (sequencing == "SE") {
    writeLines(c("#!/bin/bash", 
                 "#SBATCH -t 0-02:00", 
                 "#SBATCH -p short # Partition to run in",
                 "#SBATCH --mem 20G",
                 "#SBATCH -c 4",
                 paste("#SBATCH -e ",dir,"/",i,"_mapping_stderr.txt",sep=''),
                 paste("#SBATCH -o ",dir,"/",i,"_mapping_stdout.txt",sep=''),
                 "",
                 "module load gcc",
                 "module load fastqc conda2/4.2.13",
                 paste("cd ",dir,sep=''),
                 paste("cd ",sampDir,sep=''),
                 paste("mkdir -p mapped_",i,sep=""),
                 paste("mv ",i,"_*R*.fastq mapped_",i,sep=""),
                 paste("cd mapped_",i,sep=""),
                 paste("#FASTQC", sep = ""),
                 paste("mkdir fastqc", sep = ""),
                 paste("fastqc -o fastqc ",i,"*R1*fastq ", sep = ""),
                 
                 #load conda virtual env with updated cutadapt and bowtie2
                 "source activate py39",
                 "env",

                 # map and get methylation data with bismark
                 paste("cutadapt --minimum-length 40 -a CTNTCTCTTATACACATCT -o ",i,"_R1_t1.fastq ", i,"_R1.fastq",sep=""),
                 paste("bismark -p 4 -N 1 --non_directional --basename ",i," --genome ",genDir," ",i,"_R1_t1.fastq ",sep = ""),
                 paste("filter_non_conversion --threshold 11 ",i,".bam",sep=""),
                 
                 #sort, index
                 paste("samtools sort ",i,".nonCG_filtered.bam > ",i,".sorted.bam",sep=""),
                 paste("samtools index ",i,".sorted.bam",sep=""),
                 
                 # extract methylation levels
                 paste("bismark_methylation_extractor  --parallel 4 --ignore 10 --ignore_3prime 10 --bedGraph --genome_folder ",genDir," ",i,".nonCG_filtered.bam",sep=""),
                 # bismark reports 
                 paste("bismark2report ",i,".nonCG_filtered.bam",sep=""),
                 paste("bismark2summary ",i,".nonCG_filtered.bam",sep=""),
  
                 #collecting stats
                 #mapping
                 paste("PAIRS=$(grep 'Sequences analysed in total:' ",i,"_SE_report.txt | cut -f 2)", sep=""),
                 paste("MAPEFF=$(grep 'Mapping efficiency:' ",i,"_SE_report.txt | cut -f 2 | tr -d '%')", sep=""),
                 
                 # percent removed due to fully methylated artifact
                 paste("NONCON=$(grep 'Sequences removed' ",i,".non-conversion_filtering.txt | cut -d ' ' -f 15 | tr -d '(' | tr -d ')' |tr -d '%' )", sep=""),
                 
                 # methylation percent
                 paste("MCPG=$(grep 'C methylated in CpG context:' ",i,".nonCG_filtered_splitting_report.txt | cut -f 2 | tr -d '%')", sep=""),
                 paste("MCHH=$(grep 'C methylated in CHH context:' ",i,".nonCG_filtered_splitting_report.txt | cut -f 2 | tr -d '%')", sep=""),
                 
                 # ontarget analysis
                 paste("MAPPED=$(samtools view -b -F 4 ",i,".bam | samtools view -c)",sep=""),
                 paste("ONTARG=$(bedtools intersect -u -a ",i,".bam -b ",dir,"/captureAnalysis/",intFile," | samtools view -c)", sep=""), 
                 'PRCNT=$(echo "scale = 5; ${ONTARG} / ${MAPPED} * 100" | bc -l)',
                 paste("echo ",a$pool_ID,",",i,",$PAIRS,$MAPEFF,$NONCON,$MCPG,$MCHH,$MAPPED,$ONTARG,$PRCNT >> ",dir,"/runStats_",a$pool_ID,"_",intFileName,".bed.csv", sep=""),
                 
                 #cleanup
                 paste("cp ",dir,"/",i,"_mapping_stderr.txt ",i,"_mapping_stderr.txt", sep = ""),
                 paste("cp ",dir,"/",i,"_mapping_stdout.txt ",i,"_mapping_stdout.txt", sep = ""),
                 paste("rm ",i,"_mappingstats.txt ", sep = ""),
                 paste("rm ",dir,"/",i,"_mapping_stderr.txt", sep = ""),
                 paste("rm ",dir,"/",i,"_mapping_stdout.txt", sep = ""),
                 paste("mv ../",i,"_analyzeTimeSeq.sh ./",sep="")), script)
  }
  
  close(script)
}

