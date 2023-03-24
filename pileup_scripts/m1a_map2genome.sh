#! /bin/bash
# This file is for the ND5 m1a miseq mapping to genome
export TMPDIR=/lscratch/$SLURM_JOB_ID
set -o pipefail
set -e

function fail {
    echo "$@" >&2
    exit 1
    }

module load samtools/1.11         || fail "could not load samtools module"
module load STAR/2.7.6a          || fail "could not load STAR module"

cd /data/BatistaLab_NGS/m1A/fastq || fail "no such directory"
mkdir -p genomemapped_STAR2.7.6a
GENOME=/fdb/STAR_indices/2.7.6a/GENCODE/Gencode_human/release_35/genes-100

declare -A INPUTS
INPUTS[0]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC1-bw_S1_L001_R1_001.fastq
INPUTS[1]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC1-bw_S1_L001_R2_001.fastq
INPUTS[2]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC2-bw_S2_L001_R1_001.fastq
INPUTS[3]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC2-bw_S2_L001_R2_001.fastq
INPUTS[4]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC3-bw_S3_L001_R1_001.fastq
INPUTS[5]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC3-bw_S3_L001_R2_001.fastq
INPUTS[6]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC4-bw_S4_L001_R1_001.fastq
INPUTS[7]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC4-bw_S4_L001_R2_001.fastq
INPUTS[8]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC5-bw_S5_L001_R1_001.fastq
INPUTS[9]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC5-bw_S5_L001_R2_001.fastq
INPUTS[10]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC6-bw_S6_L001_R1_001.fastq
INPUTS[11]=/data/BatistaLab_NGS/m1A/fastq/ND5-Illumina-BC6-bw_S6_L001_R2_001.fastq
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}

declare -A OUTNAMES
OUTNAMES[0]=ND5_BC1_1
OUTNAMES[1]=ND5_BC1_2
OUTNAMES[2]=ND5_BC2_1
OUTNAMES[3]=ND5_BC2_2
OUTNAMES[4]=ND5_BC3_1
OUTNAMES[5]=ND5_BC3_2
OUTNAMES[6]=ND5_BC4_1
OUTNAMES[7]=ND5_BC4_2
OUTNAMES[8]=ND5_BC5_1
OUTNAMES[9]=ND5_BC5_2
OUTNAMES[10]=ND5_BC6_1
OUTNAMES[11]=ND5_BC6_2
OUTNAME=${OUTNAMES[$SLURM_ARRAY_TASK_ID]}


STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $GENOME \
    --readFilesIn $INPUT \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix genomemapped_STAR2.7.6a/$OUTNAME \

# sbatch --array=4 --partition=ccr --cpus-per-task=4 --mem=20g --gres=lscratch:20 m1a_map2genome.sh
