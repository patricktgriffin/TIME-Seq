#!/bin/bash
#SBATCH -n 1 # Request one core
#SBATCH -t 0-11:59                             
#SBATCH -p short # Partition to run in
#SBATCH --mem 30G

# This bash script was used to demultiplex TIME-Seq pools then write mapping / analysis scripts
# and submit them to the O2 computing cluster at HMS. 
# Note: You may need to adapt for your computing cluster and slurm version.

# Dependencies:
# https://github.com/najoshi/sabre
# R

usage() {
  echo "" > /dev/stderr
  echo "Usage: analyzeTimeSeq.sh <options> -d <directory> -s <sample_sheet> -p <pool_IDs> -b <clock_BedFile> -g <genomeDir>"  > /dev/stderr
  echo "" > /dev/stderr
  echo " required arguments:" > /dev/stderr
  echo "" > /dev/stderr
  echo "  -d     directory (required) this is the directory with the folders with fastq files for each pool" > /dev/stderr
  echo "  -s     format [ .csv ] (required) this is the sample sheet" > /dev/stderr
  echo "  -p     parenthesized list (required) this is a list of the pool IDs" > /dev/stderr 
  echo "  -b     bed (required) this is a bed file that must be in the directory ${DIR}/captureAnalysis" > /dev/stderr 
  echo "  -e     Either PE for paired-end sequencing or SE for single-end sequencing" > /dev/stderr 
  echo "  -g     directory where genome is contained" > /dev/stderr 
  echo "" > /dev/stderr
  exit 1
}

# getting options
while getopts d:s:p:b:e:g: option
do
case "${option}"
in
d) DIR=(${OPTARG});;    
s) SHEET=(${OPTARG});;   
p) POOLS=(${OPTARG});;   
b) BED=(${OPTARG});;   
e) SEQ=(${OPTARG});;   
g) GDIR=(${OPTARG});;   

esac
done

#DEMULTIPLEX
cd ${DIR}
gunzip ./*/*fastq.gz
module load gcc R/3.6.1
source ~/R-VersionSelected/bin/activate # this is a persional R library for HMS O2 cluster
pwd

# making barcode file for sabre demultiplexing from provided samplesheet
echo "Making Barcodes File"
for i in "${!POOLS[@]}"
do
	Rscript ~/scripts/makeBarcodeFiles.R ${DIR} ${SHEET} ${POOLS[i]}
done

# combining FASTQ files for pools if run across different lanes. 
## note: we have never sequenced across more than 2 lanes but this code would need to be updates if >2 lanes per pool. 
echo 'DEMULTIPLEXING!'
for i in "${!POOLS[@]}"
do
    cd ${DIR} 
    cd ${POOLS[i]}_*  
    if [ -f ${POOLS[i]}_combined_R1.fastq ] && [ -f ${POOLS[i]}_combined_R2.fastq ]; then
     	sabre pe -f ${POOLS[i]}_combined_R2.fastq -r ${POOLS[i]}_combined_R1.fastq -b barcodes_${POOLS[i]}.txt -u no_bc_match_${POOLS[i]}_R2.fq -w no_bc_match_${POOLS[i]}_R1.fq | tee sabre_report.txt 
    elif [ -f ${POOLS[i]}_*L002*_R1_*.fastq ] && [ -f ${POOLS[i]}_*L001*_R1_*.fastq ]; then
    	cat ${POOLS[i]}_*L001*_R1_*.fastq ${POOLS[i]}*L002*_R1_*.fastq > ${POOLS[i]}_combined_R1.fastq 
    	cat ${POOLS[i]}_*L001*_R2_*.fastq ${POOLS[i]}*L002*_R2_*.fastq > ${POOLS[i]}_combined_R2.fastq
    	sabre pe -f ${POOLS[i]}_combined_R2.fastq -r ${POOLS[i]}_combined_R1.fastq -b barcodes_${POOLS[i]}.txt -u no_bc_match_${POOLS[i]}_R2.fq -w no_bc_match_${POOLS[i]}_R1.fq | tee sabre_report.txt 
    else
    	sabre pe -f ${POOLS[i]}_*_R2_*.fastq -r ${POOLS[i]}_*_R1_*.fastq -b barcodes_${POOLS[i]}.txt -u no_bc_match_${POOLS[i]}_R2.fq -w no_bc_match_${POOLS[i]}_R1.fq | tee sabre_report.txt
	fi 	
done
echo 'DONE DEMULTIPLUEXING!'
echo 'WRITIING MAPING SCRIPTS!!'
#MAP
cd ${DIR}

# writing the demultiplexing and mapping scripts to be submitted on the cluster, also starting the runstats files
for i in "${!POOLS[@]}"
do
    cd ${DIR}
        
    if [ $SEQ == "PE" ]; then
    	echo "POOL,SAMPLE,READS_INITIAL,MAPPING_EFFICIENCY,MAPPING_EFFICIENCY_R1,MAPPING_EFFICIENCY_R2,PERCENT_NONCON,PERCENT_NONCON_R1,PERCENT_NONCON_R2,PERCENT_mCPG,PERCENT_mCHH,READS_MAPPED_FINAL,READS_MAPPED_FINAL_R1,READS_MAPPED_FINAL_R2,READS_ONTARG,READS_ONTARG_R1,READS_ONTARG_2,PERCENT_ONTARG,PERCENT_ONTARG_R1,PERCENT_ONTARG_R2" > runStats_${POOLS[i]}_${BED[i]}.csv
    else
    	echo "POOL,SAMPLE,READS_INITIAL,MAPPING_EFFICIENCY,PERCENT_NONCON,PERCENT_mCPG,PERCENT_mCHH,READS_MAPPED_FINAL,READS_ONTARG,PERCENT_ONTARG" > runStats_${POOLS[i]}_${BED[i]}.csv
    fi
    
    Rscript ~/scripts/writeScripts_analyzeTimeSeq.R ${DIR} ${SHEET} ${POOLS[i]} ${BED[i]} ${SEQ} ${GDIR}

done
echo 'DONE WRITING SCRIPTS!'

#submitting scripts to cluster
for i in "${!POOLS[@]}"
do
	DIRLIST=(./${POOLS[i]}_*/*_analyzeTimeSeq.sh)
	for i in "${!DIRLIST[@]}"
	do
    	chmod u+x ${DIRLIST[i]}
    	echo ${DIRLIST[i]}
    	sbatch ${DIRLIST[i]}
	done
done
