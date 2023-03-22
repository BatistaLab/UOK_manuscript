#! /bin/bash
set -e

module load samtools

cd /data/maligireddyss/18s_genomic_map
mkdir -p pileup_files

declare -A INPUTS
INPUTS[0]=18s_BC8_1Aligned.sortedByCoord.out.bam 18s_BC8_2Aligned.sortedByCoord.out.bam
INPUTS[1]=18s_BC9_1Aligned.sortedByCoord.out.bam 18s_BC9_2Aligned.sortedByCoord.out.bam
INPUTS[2]=18s_BC10_1Aligned.sortedByCoord.out.bam 18s_BC10_2Aligned.sortedByCoord.out.bam 
INPUTS[3]=18s_BC11_1Aligned.sortedByCoord.out.bam 18s_BC11_2Aligned.sortedByCoord.out.bam
INPUTS[4]=18s_BC12_1Aligned.sortedByCoord.out.bam 18s_BC12_2Aligned.sortedByCoord.out.bam
INPUTS[5]=18s_BC13_1Aligned.sortedByCoord.out.bam 18s_BC13_2Aligned.sortedByCoord.out.bam
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}

declare -A OUTFILES
OUTFILES[0]=18s_BC8.pileup
OUTFILES[1]=18s_BC9.pileup
OUTFILES[2]=18s_BC10.pileup
OUTFILES[3]=18s_BC11.pileup
OUTFILES[4]=18s_BC12.pileup
OUTFILES[5]=18s_BC13.pileup
OUTPUT=${OUTFILES[$SLURM_ARRAY_TASK_ID]}

bcftools mpileup -d 100000000 $INPUT -f /data/maligireddyss/repeat_index_star2.7.6a/repeatRNA.fa -a INFO/AD > $OUTPUT