#! /bin/bash
set -e

module load bcl2fastq/2.20.0 || exit 1
bcl2fastq --runfolder-dir /data/BatistaLab_NGS/Riboseq3/220425_VH00687_49_AAAYTCLM5 \
          -r 4 -w 4 -p 14 \
          --barcode-mismatches 0

# sbatch --cpus-per-task=20 --partition=ccr

