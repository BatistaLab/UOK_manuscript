#!/bin/bash

#sbatch --mem=10g --cpus-per-task=8

module load bedtools

bedtools genomecov -ibam BC2Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.91 > BC2.plus.strand.split.bedgraph
bedtools genomecov -ibam BC2Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.91 > BC2.minus.strand.split.bedgraph

bedtools genomecov -ibam BC4Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 1.05 > BC4.plus.split.scale.bedgraph
bedtools genomecov -ibam BC4Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 1.05 > BC4.minus.split.scale.bedgraph

bedtools genomecov -ibam BC5Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.76 > BC5.plus.split.scale.bedgraph
bedtools genomecov -ibam BC5Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.76 > BC5.minus.split.scale.bedgraph

bedtools genomecov -ibam BC6Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.81 > BC6.plus.split.scale.bedgraph
bedtools genomecov -ibam BC6Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.81 > BC6.minus.split.scale.bedgraph

bedtools genomecov -ibam BC7Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.5 > BC7.plus.split.scale.bedgraph
bedtools genomecov -ibam BC7Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.5 > BC7.minus.split.scale.bedgraph

bedtools genomecov -ibam BC12Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.89 > BC12.plus.split.scale.bedgraph
bedtools genomecov -ibam BC12Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.89 > BC12.minus.split.scale.bedgraph

bedtools genomecov -ibam BC13Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 1.01 > BC13.plus.split.scale.bedgraph
bedtools genomecov -ibam BC13Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 1.01 > BC13.minus.split.scale.bedgraph

bedtools genomecov -ibam BC14Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 1.13 > BC14.plus.split.scale.bedgraph
bedtools genomecov -ibam BC14Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 1.13 > BC14.minus.split.scale.bedgraph

bedtools genomecov -ibam BC15Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 1.14 > BC15.plus.split.scale.bedgraph
bedtools genomecov -ibam BC15Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 1.14 > BC15.minus.split.scale.bedgraph

bedtools genomecov -ibam BC16Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.86 > BC16.plus.split.scale.bedgraph
bedtools genomecov -ibam BC16Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.86 > BC16.minus.split.scale.bedgraph

bedtools genomecov -ibam BC18Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 2.34 > BC18.plus.split.scale.bedgraph
bedtools genomecov -ibam BC18Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 2.34 > BC18.minus.split.scale.bedgraph

bedtools genomecov -ibam BC19Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 1.09 > BC19.plus.split.scale.bedgraph
bedtools genomecov -ibam BC19Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 1.09 > BC19.minus.split.scale.bedgraph
