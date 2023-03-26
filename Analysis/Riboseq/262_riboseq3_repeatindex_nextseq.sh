#! /bin/bash
# This file is for the UOK262 mapping to repeat index
export TMPDIR=/lscratch/$SLURM_JOB_ID
set -o pipefail
set -e

function fail {
    echo "$@" >&2
    exit 1
    }

module load samtools/1.11         || fail "could not load samtools module"
module load STAR/2.7.6a          || fail "could not load STAR module"

cd /data/mandlerm/Riboseq3/splitFastq || fail "no such directory"
mkdir -p STAR2.7.6a_repeat
GENOME=/data/BatistaLab_NGS/annotations/repeat_Index_STAR2.7.6a

declare -A INPUTS
INPUTS[0]=/data/mandlerm/Riboseq3/splitFastq/2.collapsed.trimm.fastq
INPUTS[1]=/data/mandlerm/Riboseq3/splitFastq/4.collapsed.trimm.fastq
INPUTS[2]=/data/mandlerm/Riboseq3/splitFastq/5.collapsed.trimm.fastq
INPUTS[3]=/data/mandlerm/Riboseq3/splitFastq/6.collapsed.trimm.fastq
INPUTS[4]=/data/mandlerm/Riboseq3/splitFastq/7.collapsed.trimm.fastq
INPUTS[5]=/data/mandlerm/Riboseq3/splitFastq/12.collapsed.trimm.fastq
INPUTS[6]=/data/mandlerm/Riboseq3/splitFastq/13.collapsed.trimm.fastq
INPUTS[7]=/data/mandlerm/Riboseq3/splitFastq/14.collapsed.trimm.fastq
INPUTS[8]=/data/mandlerm/Riboseq3/splitFastq/15.collapsed.trimm.fastq
INPUTS[9]=/data/mandlerm/Riboseq3/splitFastq/16.collapsed.trimm.fastq
INPUTS[10]=/data/mandlerm/Riboseq3/splitFastq/18.collapsed.trimm.fastq
INPUTS[11]=/data/mandlerm/Riboseq3/splitFastq/19.collapsed.trimm.fastq
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


STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $GENOME \
    --readFilesIn $INPUT \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix STAR2.7.6a/$OUTNAME \

# sbatch --array=0-11 --partition=ccr --cpus-per-task=4 --mem=10g --gres=lscratch:20 262_riboseq3_repeatindex_nextseq.sh
