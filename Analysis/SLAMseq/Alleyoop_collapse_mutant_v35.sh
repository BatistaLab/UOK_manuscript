#!/bin/bash
set -e

module load singularity
  export SINGULARITY_BINDPATH="/data/BatistaLab_NGS/"
singularity exec docker://tobneu/slamdunk alleyoop collapse /data/BatistaLab_NGS/UOK_manuscript/SLAMseq/02_SLAMDUNK_pipeline/gencodeV35_multimap_mut/count/*.tsv \
-t $SLURM_CPUS_PER_TASK \
-o /data/BatistaLab_NGS/UOK_manuscript/SLAMseq/03_alleyoop/read_collapse/collapse_mutant
