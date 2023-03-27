#! /bin/bash
# author = fitzsimmonscm
# This script is designed to count RNAseq libraries prepared with the TruSeq library method
set -e

declare -A OUTNAMES
OUTNAMES[0]=HCA
OUTNAMES[1]=HCB
OUTNAMES[2]=HCC
OUTNAMES[3]=HFA
OUTNAMES[4]=HFB
OUTNAMES[5]=HFC
OUTNAMES[6]=HWA
OUTNAMES[7]=HWB
OUTNAMES[8]=HWC
OUTNAMES[9]=NCA
OUTNAMES[10]=NCB
OUTNAMES[11]=NCC
OUTNAMES[12]=NFA
OUTNAMES[13]=NFB
OUTNAMES[14]=NFC
OUTNAMES[15]=NWA
OUTNAMES[16]=NWB
OUTNAMES[17]=NWC
OUTNAME=${OUTNAMES[$SLURM_ARRAY_TASK_ID]}

module load htseq
htseq-count -m intersection-nonempty -s reverse -f bam -r pos -t exon -i gene_id /data/BatistaLab_NGS/UOK_manuscript/RNAseq/gencodeV35/picard_remove_dup/$OUTNAME.remove.dup.bam /data/BatistaLab_NGS/annotations/v35/gencode.v35.annotation.gtf > /data/BatistaLab_NGS/UOK_manuscript/RNAseq/gencodeV35/262_counts/$OUTNAME.picard.v35.count


# FLAGS
# -m = mode (options = union / intersection-strict / intersection-nonempty)
# -s = strand-specific (options = yes / no / reverse)
# -r = how data are sorted (position / name)
# -f = format of data (BAM / SAM)
# -t = feature type (3rd column of GFF file; e.g. exons)
# -i = GFF atribute to use as feature ID (default = gene_id)
