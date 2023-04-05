#!/bin/bash

#sbatch --mem=10g --cpus-per-task=8

module load bedtools

bedtools genomecov -ibam BC02Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.31 > riboBC2.plus.split.scale.bedgraph
bedtools genomecov -ibam BC02Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.31 > riboBC2.minus.split.scale.bedgraph

bedtools genomecov -ibam BC04Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.61 > riboBC4.plus.split.scale.bedgraph
bedtools genomecov -ibam BC04Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.61 > riboBC4.minus.split.scale.bedgraph

bedtools genomecov -ibam BC05Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.24 > riboBC5.plus.split.scale.bedgraph
bedtools genomecov -ibam BC05Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.24 > riboBC5.minus.split.scale.bedgraph

bedtools genomecov -ibam BC06Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.60 > riboBC6.plus.split.scale.bedgraph
bedtools genomecov -ibam BC06Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.60 > riboBC6.minus.split.scale.bedgraph

bedtools genomecov -ibam BC07Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.53 > riboBC7.plus.split.scale.bedgraph
bedtools genomecov -ibam BC07Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.53 > riboBC7.minus.split.scale.bedgraph

bedtools genomecov -ibam BC12Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 1.07 > riboBC12.plus.split.scale.bedgraph
bedtools genomecov -ibam BC12Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 1.07 > riboBC12.minus.split.scale.bedgraph

bedtools genomecov -ibam BC13Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.29 > riboBC13.plus.split.scale.bedgraph
bedtools genomecov -ibam BC13Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.29 > riboBC13.minus.split.scale.bedgraph

bedtools genomecov -ibam BC14Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.33 > riboBC14.plus.split.scale.bedgraph
bedtools genomecov -ibam BC14Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.33 > riboBC14.minus.split.scale.bedgraph

bedtools genomecov -ibam BC15Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.43 > riboBC15.plus.split.scale.bedgraph
bedtools genomecov -ibam BC15Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.43 > riboBC15.minus.split.scale.bedgraph

bedtools genomecov -ibam BC16Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.78 > riboBC16.plus.split.scale.bedgraph
bedtools genomecov -ibam BC16Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.78 > riboBC16.minus.split.scale.bedgraph

bedtools genomecov -ibam BC18Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.74 > riboBC18.plus.split.scale.bedgraph
bedtools genomecov -ibam BC18Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.74 > riboBC18.minus.split.scale.bedgraph

bedtools genomecov -ibam BC19Aligned.sortedByCoord.out.bam -bg -strand + -split -scale 0.84 > riboBC19.plus.split.scale.bedgraph
bedtools genomecov -ibam BC19Aligned.sortedByCoord.out.bam -bg -strand - -split -scale 0.84 > riboBC19.minus.split.scale.bedgraph
