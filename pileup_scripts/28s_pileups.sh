#! /bin/bash
set -e

module load samtools

cd /data/maligireddyss/28s_repeat_map
mkdir -p pileup_files

declare -A INPUTS
INPUTS[0]="28s_BC14_1Aligned.sortedByCoord.out.bam 28s_BC14_2Aligned.sortedByCoord.out.bam"
INPUTS[1]="28s_BC15_1Aligned.sortedByCoord.out.bam 28s_BC15_2Aligned.sortedByCoord.out.bam"
INPUTS[2]="28s_BC16_1Aligned.sortedByCoord.out.bam 28s_BC16_2Aligned.sortedByCoord.out.bam"
INPUTS[3]="28s_BC18_1Aligned.sortedByCoord.out.bam 28s_BC18_2Aligned.sortedByCoord.out.bam"
INPUTS[4]="28s_BC19_1Aligned.sortedByCoord.out.bam 28s_BC19_2Aligned.sortedByCoord.out.bam"
INPUTS[5]="28s_BC20_1Aligned.sortedByCoord.out.bam 28s_BC20_2Aligned.sortedByCoord.out.bam"
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}

declare -A OUTFILES
OUTFILES[0]=28s_BC14.pileup
OUTFILES[1]=28s_BC15.pileup
OUTFILES[2]=28s_BC16.pileup
OUTFILES[3]=28s_BC18.pileup
OUTFILES[4]=28s_BC19.pileup
OUTFILES[5]=28s_BC20.pileup
OUTPUT=${OUTFILES[$SLURM_ARRAY_TASK_ID]}

bcftools mpileup -d 100000000 $INPUT --fasta-ref /data/maligireddyss/repeat_index_star2.7.6a/repeatRNA.fa  -a INFO/AD > $OUTPUT
