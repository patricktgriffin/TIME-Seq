# TIME-Seq Sample Processing and Clock Analysis

This github contains example data and code used in Griffin et al. "TIME-Seq Enables Scalable and Inexpensive Epigenetic Age Predictions" (pre-print: https://doi.org/10.1101/2021.10.25.465725).

**TIME-Seq is a method for highly-scalable  DNA methylation sequencing compatible with efficient in-solution hybridization enrichment.** We have leveraged TIME-Seq to make efficient epigenetic clocks, which are up to 100-fold less expensive and considerably less laborious that BeadChip-based clocks. This method has enabled us to profile epigenetic age in thousands of human and mouse samples from a variety of tissues, cell-types, and intervnetion contexts. 

TIME-Seq relies on Tn5 transposition of a specially designed barcoded adaptor set comptible with bisulfite conversion and high-efficiency hybridization enrichment. 

<img width="388" alt="image" src="https://user-images.githubusercontent.com/94640617/232860614-6e23660e-c1a9-4e44-a9c7-526bf1351902.png">

_Here we the following :_

(1) The sample processing pipeline for demultiplexing TIME-Seq data from fastq, mapping data, and calling methylation. This pipeline uses a sample sheet with barcode identifyiers for each sample (example provided) to demultiplex raw fastq files based on the TIME-Seq barcode that is contained in Read 2. Once demultiplexed, samples are processed with a relatively standard pipeline using bismark to map reads (using bowtie2) and call methylation status. 

(2) R code that can be used to analyze TIME-Seq-based epigenetic clocks from bismark-based DNAme data. This code multiplies coefficients by methylation percentages reported by Bismark (0-100), sums the weighted methylation, adds the intercept and then applies model adjustments coefficients a and c.

(3) Current TIME-Seq clock loci and coefficients, including:
    - Mouse Multi-tissue Clock
    - Mouse Blood Clock (version 1.1 and 1.2)
    - Mouse Skin Clock
    - Mouse Liver Clock
    - Human Blood Clock

(4) Example TIME-Seq data, samplesheet, processed data, and epigenetic age predictions.


