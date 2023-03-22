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
mkdir -p 28s_repeatmapped_STAR2.7.6a
GENOME=/data/BatistaLab_NGS/annotations/repeat_Index_STAR2.7.6a

declare -A INPUTS
INPUTS[0]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC14-bw_S13_L001_R1_001.fastq
INPUTS[1]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC14-bw_S13_L001_R2_001.fastq
INPUTS[2]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC15-bw_S14_L001_R1_001.fastq
INPUTS[3]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC15-bw_S14_L001_R2_001.fastq
INPUTS[4]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC16-bw_S15_L001_R1_001.fastq
INPUTS[5]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC16-bw_S15_L001_R2_001.fastq
INPUTS[6]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC18-bw_S16_L001_R1_001.fastq
INPUTS[7]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC18-bw_S16_L001_R2_001.fastq
INPUTS[8]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC19-bw_S17_L001_R1_001.fastq
INPUTS[9]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC19-bw_S17_L001_R2_001.fastq
INPUTS[10]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC20-bw_S18_L001_R1_001.fastq
INPUTS[11]=/data/BatistaLab_NGS/m1A/fastq/28s-Illumina-BC20-bw_S18_L001_R2_001.fastq
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}

declare -A OUTNAMES
OUTNAMES[0]=28s_BC14_1
OUTNAMES[1]=28s_BC14_2
OUTNAMES[2]=28s_BC15_1
OUTNAMES[3]=28s_BC15_2
OUTNAMES[4]=28s_BC16_1
OUTNAMES[5]=28s_BC16_2
OUTNAMES[6]=28s_BC18_1
OUTNAMES[7]=28s_BC18_2
OUTNAMES[8]=28s_BC19_1
OUTNAMES[9]=28s_BC19_2
OUTNAMES[10]=28s_BC20_1
OUTNAMES[11]=28s_BC20_2
OUTNAME=${OUTNAMES[$SLURM_ARRAY_TASK_ID]}


STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $GENOME \
    --readFilesIn $INPUT \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix 28s_repeatmapped_STAR2.7.6a/$OUTNAME \

# sbatch --array=0-11 --partition=ccr --cpus-per-task=4 --mem=50g --gres=lscratch:20 28s_m1a_map2repeatindex.sh
