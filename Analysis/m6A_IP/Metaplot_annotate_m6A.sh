#! /bin/bash
# author: Christina Fitzsimmons
# last updated: 2022-08-12
# The purpose of this file is to annotate the m6A peaklist with metagene regions (TSS / 5UTR / Exon / 3UTR)
# the m6A bed files are approx 100 bp. Consider re-running using single nt bed files generated upstream of merge

set -e


declare -A INPUTS
INPUTS[0]='/data/BatistaLab_NGS/UOK_manuscript/m6AIP/metagene_analysis/2022.08/Mut.sorted.bed'
INPUTS[1]='/data/BatistaLab_NGS/UOK_manuscript/m6AIP/metagene_analysis/2022.08/WT.sorted.bed'
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}


declare -A OUTNAMES
OUTNAMES[0]='mutant'
OUTNAMES[1]='wildtype'
OUTNAME=${OUTNAMES[$SLURM_ARRAY_TASK_ID]}

module load perl
module load bedtools
perl /data/BatistaLab_NGS/UOK_manuscript/m6AIP/metagene_analysis/metaPlotR/annotate_bed_file.pl --bed $INPUT --bed2 /data/BatistaLab_NGS/UOK_manuscript/m6AIP/metagene_analysis/hg38_sorted_metagene_annot.bed > annot_m6A_sorted_$OUTNAME.bed
