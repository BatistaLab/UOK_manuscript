#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=12
#SBATCH --gres=lscratch:32
#SBATCH --job-name=ADAPTOR_TRIMMING
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

cd $SLURM_SUBMIT_DIR

module load cutadapt
module load umitools

directory=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" samples_to_map.txt)
read1=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$2}" samples_to_map.txt)
sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$4}" samples_to_map.txt)



cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 20 -q 10 -j 12 --discard-untrimmed -o step1_${sample}.fastq.gz --quiet ${directory}/${read1}
umi_tools extract --3prime --stdin=step1_${sample}.fastq.gz --bc-pattern=NNNNNNNNNNNN --stdout=step2_${sample}.fastq.gz
cutadapt -g GTTCAGAGTTCTACAGTCCGACGATC -m 20 -j 12 -o step3_${sample}.fastq.gz --quiet step2_${sample}.fastq.gz
cutadapt -a AACTGTAGGCACCATCAAT -m 10 -q 10 -j 12 --discard-untrimmed -o Processed_${sample}.fastq.gz --quiet step3_${sample}.fastq.gz

