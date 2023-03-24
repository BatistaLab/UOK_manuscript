#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=60GB
#SBATCH --cpus-per-task=15
#SBATCH --job-name=mature_miRNA
#SBATCH --error=%x_%A_%a.err
#SBATCH --output=%x_%A_%a.out

cd $SLURM_SUBMIT_DIR

module load STAR

directory=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" samples_to_map.txt)
read1=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$2}" samples_to_map.txt)
read2=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$3}" samples_to_map.txt)
sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$4}" samples_to_map.txt)
size=100

#mature miRNAs

STAR \
    --readFilesIn ${sample}_matureMIRNA.fastq \
	--genomeDir /fdb/STAR_indices/2.7.6a/UCSC/hg38/genes-${size} \
	--runThreadN 14 --genomeLoad NoSharedMemory \
	--sjdbGTFfile /data/maligireddyss/datasets/hsa.gff3 \
	--sjdbGTFfeatureExon miRNA --sjdbGTFtagExonParentTranscript ID --sjdbGTFtagExonParentGene Name \
	--outFileNamePrefix mature_miRNA/${sample}_matureseq_MIR \
	--alignEndsType EndToEnd \
	--outFilterMismatchNmax 1 \
	--outFilterMultimapScoreRange 0 \
	--quantMode TranscriptomeSAM GeneCounts \
	--outReadsUnmapped Fastx \
	--outSAMtype BAM SortedByCoordinate \
	--outFilterMultimapNmax 10 \
	--outSAMunmapped Within \
	--outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 \
	--outFilterMatchNmin 16 \
	--alignSJDBoverhangMin 1000 \
	--alignIntronMax 1 \
	--outWigType wiggle \
	--outWigStrand Stranded \
	--outWigNorm RPM \


 
