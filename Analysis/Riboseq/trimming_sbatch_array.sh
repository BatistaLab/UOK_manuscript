#!/bin/bash
set -e

declare -A INPUTS
INPUTS[0]=2
INPUTS[1]=4
INPUTS[2]=5
INPUTS[3]=6
INPUTS[4]=7
INPUTS[5]=12
INPUTS[6]=13
INPUTS[7]=14
INPUTS[8]=15
INPUTS[9]=16
INPUTS[10]=18
INPUTS[11]=19
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}

perl /data/BatistaLab_NGS/Scripts/icSHAPE-master/scripts/trimming.pl \
  -U /data/BatistaLab_NGS/Riboseq3/Nextseq_2000/splitFastq/$INPUT.collapsed.fastq \
  -o /data/BatistaLab_NGS/Riboseq3/Nextseq_2000/splitFastq/$INPUT.trimm.collapsed.fastq \
  -l 17 \
  -t 0 \
  -c phred33 \
  -a /data/BatistaLab_NGS/Scripts/icSHAPE-master/data/adapter/TruSeq2-PE.fa \
  -m 15


# command line input
# sbatch --array=11 --partition=ccr --cpus-per-task=8 --mem=20g trimming_sbatch_array.sh

# make sure the .sh file is in the same folder as the files, the files are numbered and then run using the command below
# sbatch --array=3,4,5 --partition=ccr --cpus-per-task=4 --mem=10g trimming_sbatch_array.sh
