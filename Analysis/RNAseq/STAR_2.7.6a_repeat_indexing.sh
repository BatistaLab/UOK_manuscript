#! /bin/bash
# author: Christina Fitzsimmons
# This script generates a repeat index for STAR/2.7.6a using hg19 and PB repeat RNA fasta file
set -o pipefail
set -e

function fail {
    echo "$@" >&2
    exit 1
    }

module load samtools/1.11         || fail "could not load samtools module"
module load STAR/2.7.6a          || fail "could not load STAR module"

cd /data/BatistaLab_NGS/annotations/
mkdir -p repeat_Index_STAR2.7.6a

STAR \
  --runThreadN $SLURM_CPUS_PER_TASK \
  --runMode genomeGenerate \
  --genomeDir repeat_Index_STAR2.7.6a \
  --genomeFastaFiles /data/BatistaLab_NGS/annotations/repeat_Index/repeatRNA.fa \

# make sure all gtf and fasta files are unzipped or this mode will fail!
# Note! STAR/2.5.4 indices are not compatible with newer versions of STAR.
