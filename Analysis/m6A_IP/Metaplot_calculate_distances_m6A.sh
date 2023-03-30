#! /bin/bash
# author: Christina Fitzsimmons
# last updated: 2022-08-12
# The purpose of this file is to calculate distances on the m6A metaplot

set -e


declare -A INPUTS
INPUTS[0]='/data/BatistaLab_NGS/UOK_manuscript/m6AIP/metagene_analysis/2022.08/annot_m6A_sorted_mutant.bed'
INPUTS[1]='/data/BatistaLab_NGS/UOK_manuscript/m6AIP/metagene_analysis/2022.08/annot_m6A_sorted_wildtype.bed'
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}


declare -A OUTNAMES
OUTNAMES[0]='mutant'
OUTNAMES[1]='wildtype'
OUTNAME=${OUTNAMES[$SLURM_ARRAY_TASK_ID]}

module load perl
module load bedtools
perl /data/BatistaLab_NGS/UOK_manuscript/m6AIP/metagene_analysis/metaPlotR/rel_and_abs_dist_calc.pl --bed $INPUT --regions /data/BatistaLab_NGS/UOK_manuscript/m6AIP/metagene_analysis/region_sizes.txt > $OUTNAME.m6A-distance-measures.txt
