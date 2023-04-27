# TIME-Seq Analysis 

This github contains example data and code used in Griffin et al. "TIME-Seq Enables Scalable and Inexpensive Epigenetic Age Predictions" (pre-print: https://doi.org/10.1101/2021.10.25.465725).

**TIME-Seq is a method for highly-scalable  DNA methylation sequencing compatible with efficient in-solution hybridization enrichment.** We have leveraged TIME-Seq to make efficient epigenetic clocks, which are up to 100-fold less expensive and considerably less laborious that BeadChip-based clocks. This method has enabled us to profile epigenetic age in thousands of human and mouse samples from a variety of tissues, cell-types, and intervnetion contexts. 

TIME-Seq relies on Tn5 transposition of a specially designed barcoded adaptor set comptible with bisulfite conversion and high-efficiency hybridization enrichment. 

<img width="388" alt="image" src="https://user-images.githubusercontent.com/94640617/232860614-6e23660e-c1a9-4e44-a9c7-526bf1351902.png">

_Here we include the following :_

(1) The sample processing pipeline for demultiplexing TIME-Seq data from fastq, mapping data, and calling methylation. 

This pipeline uses a sample sheet with barcode identifyiers for each sample (example provided) to demultiplex raw fastq files based on the TIME-Seq barcode that is contained in Read 2. Once demultiplexed, samples are processed with a relatively standard pipeline using bismark to map reads (using bowtie2) and call methylation status. 

(2) R code that can be used to analyze TIME-Seq-based epigenetic clocks from bismark-based DNAme data. 

This code multiplies coefficients by methylation percentages reported by Bismark (0-100), sums the weighted methylation, adds the intercept and then applies model adjustments coefficients a and c.

(3) Current TIME-Seq clock loci and coefficients, including:

    - Mouse Multi-tissue Clock
    
    - Mouse Blood Clock (version 1.1 and 1.2)
    
    - Mouse Skin Clock
    
    - Mouse Liver Clock
    
    - Human Blood Clock

(4) Example TIME-Seq data, samplesheet, processed data, and epigenetic age predictions.

_____________

## Usage
The BASH script “analyzeTimeSeq.sh” has been used for sample demultiplexing and subsequently writing mapping / methylation calling scripts for each sample on a computing cluster like the O2 cluster at Harvard Medical School. 

Please note: this pipeline is included as it was used in the TIME-Seq manuscript. You will likely need to adapt the pipeline or individual commands to fit with the computational platform you are working on.

The submission of the BASH script requires several inputs that should be included in quotes, described in the usage output at the top of the script. 

```
analyzeTimeSeq.sh <options> -d "<directory>" -s "<sample_sheet>" -p "<pool_IDs>" -b "<clock_BedFile>" -g "<genDir>"

option: -d     directory (required) this is the directory with the folders with fastq files for each pool
option: -s     format [ .csv ] (required) this is the name of the sample sheet. An examples samplesheet is provided: “example_demultiplexing_samplesheet.csv”
option: -p     parenthesized list (required) this is a list of the pool IDs. Separate pool IDs with a space.
option: -b     bed (required) this is a bed file that must be in this directory "../sinclair/Patrick/methylationClock/captureAnalysis" 
option: -e     Either PE for paired-end sequencing or SE for single-end sequencing
option: -g     directory where genome is contained
```

For each pool name that you provide as input, the script writes a barcode file that works with the sabre software.

Next, sabre (https://github.com/najoshi/sabre) demultiplexes samples into individual FASTQ files based on the barcode file. 
Example barcode file provided “barcodes_pool.txt”

After demultiplexing is done, the samples in each pool have a script written for mapping / methylation calling using the “writeScripts_analyzeTimeSeq.R” script. 
This R script uses the provided info and the sample sheet and writes 1 script for each sample in the pool. 

In addition to mapping, processing, and methylation calling (Bismark), the script written by “writeScripts_analyzeTimeSeq.R” also collects and print sample processing stats and overlap with target loci (provided as a .bed file).

The genomes that we used were downloaded from iGenome (https://support.illumina.com/sequencing/sequencing_software/igenome.html). If using the script without altering it, you need to reference your own genome directory that contains the genome fasta that then has been BS converted (bismark_genome_preparation) and prepared into a bowtie2 index that is automatically named: "Bisulfite_Genome".

The bed files that you want to intersect for on-target analysis of hybrid capture need to be contained in the base directory in a subdirectory called captureAnalysis.

For any help or questions please contact ptgriffin {at} g.harvard.edu 
