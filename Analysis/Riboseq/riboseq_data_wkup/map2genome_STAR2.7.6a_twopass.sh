#! /bin/bash
# This file is for the UOK262 mapping to gecode v35 genome with 100 bp overhangs
export TMPDIR=/lscratch/$SLURM_JOB_ID
set -o pipefail
set -e

function fail {
    echo "$@" >&2
    exit 1
    }

module load samtools/1.11         || fail "could not load samtools module"
module load STAR/2.7.6a          || fail "could not load STAR module"

cd /data/mandlerm/Riboseq3/splitFastq/STAR2.7.6a || fail "no such directory"
mkdir -p genomemapped_STAR2.7.6a
GENOME=/fdb/STAR_indices/2.7.6a/GENCODE/Gencode_human/release_35/genes-100

declare -A INPUTS
INPUTS[0]=BC2Unmapped.out.mate1
INPUTS[1]=BC4Unmapped.out.mate1
INPUTS[2]=BC5Unmapped.out.mate1
INPUTS[3]=BC6Unmapped.out.mate1
INPUTS[4]=BC7Unmapped.out.mate1
INPUTS[5]=BC12Unmapped.out.mate1
INPUTS[6]=BC13Unmapped.out.mate1
INPUTS[7]=BC14Unmapped.out.mate1
INPUTS[8]=BC15Unmapped.out.mate1
INPUTS[9]=BC16Unmapped.out.mate1
INPUTS[10]=BC18Unmapped.out.mate1
INPUTS[11]=BC19Unmapped.out.mate1
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
    --readFilesIn /data/mandlerm/Riboseq3/splitFastq/STAR2.7.6a/$INPUT \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix genomemapped_STAR2.7.6a/$OUTNAME \
    --twopassMode Basic\
# sbatch --array=0-11 --partition=ccr --cpus-per-task=4 --mem=50g --gres=lscratch:20 262_riboseq3_map2genome_STAR2.7.6a_nextseq.sh
