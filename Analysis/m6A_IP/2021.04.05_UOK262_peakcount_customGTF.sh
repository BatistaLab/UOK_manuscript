#! /bin/bash
# author = fitzsimmonscm
# last updated 2021-04-02
set -e

declare -A INPUTS
INPUTS[0]=BC2Aligned
INPUTS[1]=BC4Aligned
INPUTS[2]=BC5Aligned
INPUTS[3]=BC6Aligned
INPUTS[4]=BC7Aligned
INPUTS[5]=BC12Aligned
INPUTS[6]=BC13Aligned
INPUTS[7]=BC14Aligned
INPUTS[8]=BC15Aligned
INPUTS[9]=BC16Aligned
INPUTS[10]=BC18Aligned
INPUTS[11]=BC19Aligned
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}

declare -A OUTNAMES
OUTNAMES[0]=BC2
OUTNAMES[1]=BC4
OUTNAMES[2]=BC5
OUTNAMES[3]=BC6
OUTNAMES[4]=BC7
OUTNAMES[5]=BC12
OUTNAMES[6]=BC13
OUTNAMES[7]=BC14
OUTNAMES[8]=BC15
OUTNAMES[9]=BC16
OUTNAMES[10]=BC18
OUTNAMES[11]=BC19
OUTNAME=${OUTNAMES[$SLURM_ARRAY_TASK_ID]}

module load htseq
htseq-count -m union -s yes -f bam -r pos -t exon -i gene_id /data/BatistaLab_NGS/UOK_manuscript/m6AIP/262_genome_map/STAR2.7.6a/$INPUT.sortedByCoord.out.bam /data/BatistaLab_NGS/UOK_manuscript/m6AIP/262_readcounts/2021.04.05_m6AIP_merged-mutant-wildtype_customGFF.gff > /data/BatistaLab_NGS/UOK_manuscript/m6AIP/262_readcounts/peak_counts/2021.04.05/$OUTNAME.strand-yes.count


# FLAGS
# -m = mode (options = union / intersection-strict / intersection-nonempty)
# -s = strand-specific (options = yes / no / reverse)
# -r = how data are sorted (position / name)
# -f = format of data (BAM / SAM)
# -t = feature type (3rd column of GFF file; e.g. exons)
# -i = GFF atribute to use as feature ID (default = gene_id)
