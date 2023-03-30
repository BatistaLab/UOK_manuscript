#!/bin/bash
set -e

module load singularity
  export SINGULARITY_BINDPATH="/data/BatistaLab_NGS/"
singularity exec docker://tobneu/slamdunk alleyoop utrrates -r /data/BatistaLab_NGS/annotations/v35/GRCh38.p13.genome.fa \
  -b /data/BatistaLab_NGS/UOK_manuscript/SLAMseq/01_3utr_annotation/2021.03.20_gencodeV35_3UTR_SLAMseq_ENSG_ids.bed \
  -o /data/BatistaLab_NGS/UOK_manuscript/SLAMseq/02_SLAMDUNK_pipeline/alleyoop \
  -t 4 \
  /data/BatistaLab_NGS/UOK_manuscript/SLAMseq/02_SLAMDUNK_pipeline/gencodeV35_multimap_WT/filter/*.bam
