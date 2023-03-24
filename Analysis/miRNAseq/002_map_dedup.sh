#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=15
#SBATCH --gres=lscratch:64
#SBATCH --job-name=MAP_DEDUP
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

cd $SLURM_SUBMIT_DIR

module load STAR/2.7.6a
module load umitools/1.1.1
module load samtools/1.11

directory=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" samples_to_map.txt)
read1=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$2}" samples_to_map.txt)
read2=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$3}" samples_to_map.txt)
sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$4}" samples_to_map.txt)
size=100


#align

STAR \
    --readFilesIn Processed_${sample}.fastq.gz \
    --genomeDir /fdb/STAR_indices/2.7.6a/UCSC/hg38/genes-${size} \
    --runThreadN 14 --genomeLoad LoadAndRemove \
    --limitBAMsortRAM 20000000000 \
    --readFilesCommand zcat \
    --outFileNamePrefix ${sample}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outReadsUnmapped Fastx \
    --outFilterMismatchNmax 1 \
    --outFilterMismatchNoverLmax 1 \
    --outFilterMismatchNoverReadLmax 1 \
    --outFilterMatchNmin 16 \
    --outFilterMatchNminOverLread 0 \
    --outFilterScoreMinOverLread 0 \
    --outFilterMultimapNmax 5000 \
    --winAnchorMultimapNmax 5000 \
    --seedSearchStartLmax 30 \
    --alignTranscriptsPerReadNmax 30000 \
    --alignWindowsPerReadNmax 30000 \
    --alignTranscriptsPerWindowNmax 300 \
    --seedPerReadNmax 3000 \
    --seedPerWindowNmax 300 \
    --seedNoneLociPerWindow 1000 \
    --outFilterMultimapScoreRange 0 \
    --alignIntronMax 1 \
    --alignSJDBoverhangMin 999999999999


#sort and index
samtools sort ${sample}_Aligned.sortedByCoord.out.bam -o ${sample}.sorted.bam
samtools index ${sample}.sorted.bam

#dedup
umi_tools \
    dedup -I ${sample}.sorted.bam \
    -S ${sample}.dedup.bam \
    --log=${sample}.dedup_log.txt \
	--method=unique

#bam to fastq
samtools fastq -@ 12 ${sample}.dedup.bam > ${sample}.dedup.fastq
