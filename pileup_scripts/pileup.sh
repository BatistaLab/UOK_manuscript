#! /bin/bash
set -e

module load samtools

cd /data/maligireddyss/ND5_pileup_bams
mkdir -p pileup_files

declare -A INPUTS
INPUTS[0]=ND5_BC1_1Aligned.sortedByCoord.out.bam
INPUTS[1]=ND5_BC1_2Aligned.sortedByCoord.out.bam
INPUTS[2]=ND5_BC2_1Aligned.sortedByCoord.out.bam
INPUTS[3]=ND5_BC2_2Aligned.sortedByCoord.out.bam
INPUTS[4]=ND5_BC3_1Aligned.sortedByCoord.out.bam
INPUTS[5]=ND5_BC3_2Aligned.sortedByCoord.out.bam
INPUTS[6]=ND5_BC4_1Aligned.sortedByCoord.out.bam
INPUTS[7]=ND5_BC4_2Aligned.sortedByCoord.out.bam
INPUTS[8]=ND5_BC5_1Aligned.sortedByCoord.out.bam
INPUTS[9]=ND5_BC5_2Aligned.sortedByCoord.out.bam
INPUTS[10]=ND5_BC6_1Aligned.sortedByCoord.out.bam
INPUTS[11]=ND5_BC6_2Aligned.sortedByCoord.out.bam
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}

declare -A OUTFILES
OUTFILES[0]=ND5_BC1_1.pileup
OUTFILES[1]=ND5_BC1_2.pileup
OUTFILES[2]=ND5_BC2_1.pileup
OUTFILES[3]=ND5_BC2_2.pileup
OUTFILES[4]=ND5_BC3_1.pileup
OUTFILES[5]=ND5_BC3_2.pileup
OUTFILES[6]=ND5_BC4_1.pileup
OUTFILES[7]=ND5_BC4_2.pileup
OUTFILES[8]=ND5_BC5_1.pileup
OUTFILES[9]=ND5_BC5_2.pileup
OUTFILES[10]=ND5_BC6_1.pileup
OUTFILES[11]=ND5_BC6_2.pileup
OUTPUT=${OUTFILES[$SLURM_ARRAY_TASK_ID]}

bcftools mpileup -d 100000000 $INPUT -f /data/maligireddyss/pileup/GRCh38.primary_assembly.genome.fa  -r chrM:13532-13807 -a INFO/AD > $OUTPUT

