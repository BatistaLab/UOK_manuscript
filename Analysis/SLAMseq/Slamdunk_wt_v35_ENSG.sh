#!/bin/bash
# SLAMseq analysis of CMF library #2 using gencodeV35
# optional multimapper flag included
set -e

module load singularity
#SINGULARITY_BINDPATH needs the location of the reference geneome (fasta) and the raw data (fastq)
export SINGULARITY_BINDPATH="/data/BatistaLab_NGS/"
singularity exec docker://tobneu/slamdunk:latest slamdunk all /data/BatistaLab_NGS/UOK_manuscript/SLAMseq/02_SLAMDUNK_pipeline/library2_manifest_WT_v27.tsv \
  -r /data/BatistaLab_NGS/annotations/v35/GRCh38.p13.genome.fa \
  -b /data/BatistaLab_NGS/UOK_manuscript/SLAMseq/01_3utr_annotation/2021.03.20_gencodeV35_3UTR_SLAMseq_ENSG_ids.bed \
  -o /data/BatistaLab_NGS/UOK_manuscript/SLAMseq/02_SLAMDUNK_pipeline/gencodeV35_multimap_WT \
  -t $SLURM_CPUS_PER_TASK \
  -5 12 \
  -a 4 \
  -rl 500 \
  -m

  # required inputs
    # -r /path/to/reference/genome.fa
    # -o /path/to/output/directory (create before running the script)
    # -t number of threads you need. Typical starting amount are 8 threads and 20-40g memory
    # -b /path/to/3'UTRcoordinates/bedfile.bed [count]
    # -m		Use 3’UTR annotation to filter multimappers [filter]
    # A .tsv file containing the info for 1 sample per row. Do NOT include a header row. List information in the exact order below:
