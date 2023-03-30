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
mkdir -p 18s_repeatmapped_STAR2.7.6a
GENOME=/data/BatistaLab_NGS/annotations/repeat_Index_STAR2.7.6a

declare -A INPUTS
INPUTS[0]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC8-bw_S7_L001_R1_001.fastq
INPUTS[1]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC8-bw_S7_L001_R2_001.fastq
INPUTS[2]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC9-bw_S8_L001_R1_001.fastq
INPUTS[3]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC9-bw_S8_L001_R2_001.fastq
INPUTS[4]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC10-bw_S9_L001_R1_001.fastq
INPUTS[5]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC10-bw_S9_L001_R2_001.fastq
INPUTS[6]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC11-bw_S10_L001_R1_001.fastq
INPUTS[7]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC11-bw_S10_L001_R2_001.fastq
INPUTS[8]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC12-bw_S11_L001_R1_001.fastq
INPUTS[9]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC12-bw_S11_L001_R2_001.fastq
INPUTS[10]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC13-bw_S12_L001_R1_001.fastq
INPUTS[11]=/data/BatistaLab_NGS/m1A/fastq/18s-Illumina-BC13-bw_S12_L001_R2_001.fastq
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}

declare -A OUTNAMES
OUTNAMES[0]=18s_BC8_1
OUTNAMES[1]=18s_BC8_2
OUTNAMES[2]=18s_BC9_1
OUTNAMES[3]=18s_BC9_2
OUTNAMES[4]=18s_BC10_1
OUTNAMES[5]=18s_BC10_2
OUTNAMES[6]=18s_BC11_1
OUTNAMES[7]=18s_BC11_2
OUTNAMES[8]=18s_BC12_1
OUTNAMES[9]=18s_BC12_2
OUTNAMES[10]=18s_BC13_1
OUTNAMES[11]=18s_BC13_2
OUTNAME=${OUTNAMES[$SLURM_ARRAY_TASK_ID]}


STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $GENOME \
    --readFilesIn $INPUT \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix 18s_repeatmapped_STAR2.7.6a/$OUTNAME \

# sbatch --array=0-11 --partition=ccr --cpus-per-task=4 --mem=50g --gres=lscratch:20 18s_m1a_repeatindex.sh
