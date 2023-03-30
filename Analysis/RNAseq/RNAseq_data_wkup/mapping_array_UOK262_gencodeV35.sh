#! /bin/bash
# readFilesIn paired end syntax = path to read1.fastq /space/ path to read2.fastq
# This file is for the UOK262 mapping to gencode v35 indices with 100 bp overhang.
export TMPDIR=/lscratch/$SLURM_JOB_ID
set -o pipefail
set -e

function fail {
    echo "$@" >&2
    exit 1
    }

module load samtools/1.11         || fail "could not load samtools module"
module load STAR/2.7.6a         || fail "could not load STAR module"

cd /data/BatistaLab_NGS/UOK_manuscript/RNAseq/gencodeV35 || fail "no such directory"
mkdir -p 262_map_gencodeV35
GENOME=/fdb/STAR_indices/2.7.6a/GENCODE/Gencode_human/release_35/genes-100

declare -A INPUTS
INPUTS[0]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HC/trimm/HCA_S4_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HC/trimm/HCA_S4_comb_R2.trimm.fastq.gz"
INPUTS[1]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HC/trimm/HCB_S5_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HC/trimm/HCB_S5_comb_R2.trimm.fastq.gz"
INPUTS[2]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HC/trimm/HCC_S6_comb_R1.trimm2.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HC/trimm/HCC_S6_comb_R2.trimm2.fastq.gz"
INPUTS[3]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HF/trimm/HFA_S1_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HF/trimm/HFA_S1_comb_R2.trimm.fastq.gz"
INPUTS[4]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HF/trimm/HFB_S2_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HF/trimm/HFB_S2_comb_R2.trimm.fastq.gz"
INPUTS[5]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HF/trimm/HFC_S3_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HF/trimm/HFC_S3_comb_R2.trimm.fastq.gz"
INPUTS[6]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HW/trimm/HWA_S7_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HW/trimm/HWA_S7_comb_R2.trimm.fastq.gz"
INPUTS[7]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HW/trimm/HWB_S8_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HW/trimm/HWB_S8_comb_R2.trimm.fastq.gz"
INPUTS[8]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HW/trimm/HWC_S9_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/HW/trimm/HWC_S9_comb_R2.trimm.fastq.gz"
INPUTS[9]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NC/trimm/NCA_S13_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NC/trimm/NCA_S13_comb_R2.trimm.fastq.gz"
INPUTS[10]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NC/trimm/NCB_S14_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NC/trimm/NCB_S14_comb_R2.trimm.fastq.gz"
INPUTS[11]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NC/trimm/NCC_S15_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NC/trimm/NCC_S15_comb_R2.trimm.fastq.gz"
INPUTS[12]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NF/trimm/NFA_S10_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NF/trimm/NFA_S10_comb_R2.trimm.fastq.gz"
INPUTS[13]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NF/trimm/NFB_S11_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NF/trimm/NFB_S11_comb_R2.trimm.fastq.gz"
INPUTS[14]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NF/trimm/NFC_S12_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NF/trimm/NFC_S12_comb_R2.trimm.fastq.gz"
INPUTS[15]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NW/trimm/NWA_S16_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NW/trimm/NWA_S16_comb_R2.trimm.fastq.gz"
INPUTS[16]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NW/trimm/NWB_S17_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NW/trimm/NWB_S17_comb_R2.trimm.fastq.gz"
INPUTS[17]="/data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NW/trimm/NWC_S18_comb_R1.trimm.fastq.gz /data/BatistaLab_NGS/RNA_seq/UOK262_NGS/NW/trimm/NWC_S18_comb_R2.trimm.fastq.gz"
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


STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $GENOME \
    --readFilesCommand zcat \
    --sjdbOverhang 100 \
    --readFilesIn $INPUT \
    --outSAMtype BAM SortedByCoordinate \
    --outTmpDir=/lscratch/$SLURM_JOB_ID/STARtmp \
    --outReadsUnmapped Fastx \
    --outFileNamePrefix 262_map_gencodeV35/$OUTNAME \
    –-outFilterType BySJout \
    --outFilterMultimapNmax 10 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --outFilterMatchNminOverLread 0.66 \
    --quantMode TranscriptomeSAM \
    --peOverlapNbasesMin 10 \
    --alignEndsProtrude 10 ConcordantPair
