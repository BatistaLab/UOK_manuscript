#! /bin/bash
# author: fitzsimmonscm
# This script removes PCR duplicates from RNAseq data

set -e
export TMPDIR=/lscratch/$SLURM_JOB_ID
function fail {
    echo "$@" >&2
    exit 1
    }

declare -A INPUTS
INPUTS[0]=HCAAligned
INPUTS[1]=HCBAligned
INPUTS[2]=HCCAligned
INPUTS[3]=HFAAligned
INPUTS[4]=HFBAligned
INPUTS[5]=HFCAligned
INPUTS[6]=HWAAligned
INPUTS[7]=HWBAligned
INPUTS[8]=HWCAligned
INPUTS[9]=NCAAligned
INPUTS[10]=NCBAligned
INPUTS[11]=NCCAligned
INPUTS[12]=NFAAligned
INPUTS[13]=NFBAligned
INPUTS[14]=NFCAligned
INPUTS[15]=NWAAligned
INPUTS[16]=NWBAligned
INPUTS[17]=NWCAligned
INPUT=${INPUTS[$SLURM_ARRAY_TASK_ID]}

declare -A OUTNAMES
OUTNAMES[0]=HCA
OUTNAMES[1]=HCB
OUTNAMES[2]=HCC
OUTNAMES[3]=HFA
OUTNAMES[4]=HFB
OUTNAMES[5]=HFC
OUTNAMES[6]=HWA
OUTNAMES[7]=HWB
OUTNAMES[8]=HWC
OUTNAMES[9]=NCA
OUTNAMES[10]=NCB
OUTNAMES[11]=NCC
OUTNAMES[12]=NFA
OUTNAMES[13]=NFB
OUTNAMES[14]=NFC
OUTNAMES[15]=NWA
OUTNAMES[16]=NWB
OUTNAMES[17]=NWC
OUTNAME=${OUTNAMES[$SLURM_ARRAY_TASK_ID]}


module load picard || fail "could not load picard"

java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
    I=/data/BatistaLab_NGS/UOK_manuscript/RNAseq/gencodeV35/262_map_gencodeV35/$INPUT.sortedByCoord.out.bam \
    O=/data/BatistaLab_NGS/UOK_manuscript/RNAseq/gencodeV35/picard_remove_dup/$OUTNAME.remove.dup.bam \
    ASSUME_SORTED=true \
    METRICS_FILE=$OUTNAME.mapped.metric.csv \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    VALIDATION_STRINGENCY=LENIENT \
    REMOVE_DUPLICATES=TRUE
