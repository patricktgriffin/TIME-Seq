# TIME-Seq Sample Processing and Clock Analysis

This github contains example data and code used in Griffin et al. "TIME-Seq Enables Scalable and Inexpensive Epigenetic Age Predictions" (pre-print: https://doi.org/10.1101/2021.10.25.465725).

TIME-Seq is a method for highly-scalable  DNA methylation sequencing compatible with efficient in-solution hybridization enrichment. We have leveraged TIME-Seq to make efficient epigenetic clocks, which are up to 100-fold less expensive and considerably less laborious that BeadChip-based clocks. This method has enabled us to profile epigenetic age in thousands of human and mouse samples from a variety of tissues, cell-types, and intervnetion contexts. 

TIME-Seq relies on Tn5 transposition of a specially designed barcoded adaptor set comptible with bisulfite conversion and high-efficiency hybridization enrichment. 

<img width="388" alt="image" src="https://user-images.githubusercontent.com/94640617/232860614-6e23660e-c1a9-4e44-a9c7-526bf1351902.png">

Here we provide:

(1) The sample processing pipeline for demultiplexing TIME-Seq data from fastq, mapping data, and calling methylation.

(2) Current TIME-Seq clock loci and coefficients.

(3) R code that can be used to analyze TIME-Seq-based epigenetic clocks from bismark-based DNAme data. 

