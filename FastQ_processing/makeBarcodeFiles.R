#!/usr/bin/env Rscript

## This script makes barcode files formatted for demultiplexing with sabre using the provided samplesheet

suppressMessages(library("dplyr"))

args = commandArgs(trailingOnly=TRUE)
dir <- args[1]
samplesheet <- as.character(args[2])
pools <- as.list(args[3])
s <- read.csv(samplesheet, header=T) %>%
  dplyr::filter(pool_ID %in% pools)

for (i in unique(s$pool_ID)) {
  a <- dplyr::filter(s, pool_ID == i) #selecting just the row for each unique sampleID
  setwd(dir)
  dirslist <- list.dirs(recursive = F)
  sampDir <- dirslist[grep(a$pool_ID[1],dirslist)]
  setwd(sampDir)
 
  b <- a %>% #making files for sabre demultiplexing to read
    mutate(read1_fastq =  paste(sample_ID,"_R2.fastq",sep="")) %>%
    mutate(read2_fastq =  paste(sample_ID,"_R1.fastq",sep="")) %>%
    dplyr::select(barcode,read1_fastq,read2_fastq)
  wd <- getwd()
  pool <- as.character(a$pool_ID[1])
  write.table(b,file =paste(wd,"/barcodes_",pool,".txt", sep = ""),quote=F,sep="\t", row.names = F, col.names = F)
  
  }
  
