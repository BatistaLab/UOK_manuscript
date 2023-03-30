#!/bin/bash
# Christina Fitzsimmons
set -e
set -o pipefail

source /data/fitzsimmonscm/conda/etc/profile.d/conda.sh

declare -A INPUTS
INPUTS[0]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC02Aligned.sortedByCoord.out.bam"
INPUTS[1]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC04Aligned.sortedByCoord.out.bam"
INPUTS[2]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC05Aligned.sortedByCoord.out.bam"
INPUTS[3]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC06Aligned.sortedByCoord.out.bam"
INPUTS[4]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC07Aligned.sortedByCoord.out.bam"
INPUTS[5]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC12Aligned.sortedByCoord.out.bam"
INPUTS[6]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC13Aligned.sortedByCoord.out.bam"
INPUTS[7]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC14Aligned.sortedByCoord.out.bam"
INPUTS[8]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC15Aligned.sortedByCoord.out.bam"
INPUTS[9]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC16Aligned.sortedByCoord.out.bam"
INPUTS[10]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC18Aligned.sortedByCoord.out.bam"
INPUTS[11]="/data/BatistaLab_NGS/Riboseq3/Nextseq_2000/mapped/riboseq3_gencode/OnePass/BC19Aligned.sortedByCoord.out.bam"
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}

declare -A OUTNAMES
OUTNAMES[0]=BC02MANEcount
OUTNAMES[1]=BC04MANEcount
OUTNAMES[2]=BC05MANEcount
OUTNAMES[3]=BC06MANEcount
OUTNAMES[4]=BC07MANEcount
OUTNAMES[5]=BC12MANEcount
OUTNAMES[6]=BC13MANEcount
OUTNAMES[7]=BC14MANEcount
OUTNAMES[8]=BC15MANEcount
OUTNAMES[9]=BC16MANEcount
OUTNAMES[10]=BC18MANEcount
OUTNAMES[11]=BC19MANEcount
OUTNAME=${OUTNAMES[$SLURM_ARRAY_TASK_ID]}

conda activate plastid

cs count /data/BatistaLab_NGS/Riboseq3/Nextseq_2000/plastid_metagene_counts/MANE/MANE_position_gene.positions \
--count_files $INPUT \
--fiveprime \
--countfile_format BAM \
--min_length 25 \
--max_length 34 \
$OUTNAME
