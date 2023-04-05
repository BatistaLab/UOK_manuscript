#!/bin/bash

#sbatch --mem=10g --cpus-per-task=8

./bedGraphToBigWig riboBC2.plus.split.scale.bedgraph hg38.chrom.sizes riboBC2.plus.split.scale.bw
./bedGraphToBigWig riboBC2.minus.split.scale.bedgraph hg38.chrom.sizes riboBC2.minus.split.scale.bw

./bedGraphToBigWig riboBC4.plus.split.scale.bedgraph hg38.chrom.sizes riboBC4.plus.split.scale.bw
./bedGraphToBigWig riboBC4.minus.split.scale.bedgraph hg38.chrom.sizes riboBC4.minus.split.scale.bw

./bedGraphToBigWig riboBC5.plus.split.scale.bedgraph hg38.chrom.sizes riboBC5.plus.split.scale.bw
./bedGraphToBigWig riboBC5.minus.split.scale.bedgraph hg38.chrom.sizes riboBC5.minus.split.scale.bw

./bedGraphToBigWig riboBC6.plus.split.scale.bedgraph hg38.chrom.sizes riboBC6.plus.split.scale.bw
./bedGraphToBigWig riboBC6.minus.split.scale.bedgraph hg38.chrom.sizes riboBC6.minus.split.scale.bw

./bedGraphToBigWig riboBC7.plus.split.scale.bedgraph hg38.chrom.sizes riboBC7.plus.split.scale.bw
./bedGraphToBigWig riboBC7.minus.split.scale.bedgraph hg38.chrom.sizes riboBC7.minus.split.scale.bw

./bedGraphToBigWig riboBC12.plus.split.scale.bedgraph hg38.chrom.sizes riboBC12.plus.split.scale.bw
./bedGraphToBigWig riboBC12.minus.split.scale.bedgraph hg38.chrom.sizes riboBC12.minus.split.scale.bw

./bedGraphToBigWig riboBC13.plus.split.scale.bedgraph hg38.chrom.sizes riboBC13.plus.split.scale.bw
./bedGraphToBigWig riboBC13.minus.split.scale.bedgraph hg38.chrom.sizes riboBC13.minus.split.scale.bw

./bedGraphToBigWig riboBC14.plus.split.scale.bedgraph hg38.chrom.sizes riboBC14.plus.split.scale.bw
./bedGraphToBigWig riboBC14.minus.split.scale.bedgraph hg38.chrom.sizes riboBC14.minus.split.scale.bw

./bedGraphToBigWig riboBC15.plus.split.scale.bedgraph hg38.chrom.sizes riboBC15.plus.split.scale.bw
./bedGraphToBigWig riboBC15.minus.split.scale.bedgraph hg38.chrom.sizes riboBC15.minus.split.scale.bw

./bedGraphToBigWig riboBC16.plus.split.scale.bedgraph hg38.chrom.sizes riboBC16.plus.split.scale.bw
./bedGraphToBigWig riboBC16.minus.split.scale.bedgraph hg38.chrom.sizes riboBC16.minus.split.scale.bw

./bedGraphToBigWig riboBC18.plus.split.scale.bedgraph hg38.chrom.sizes riboBC18.plus.split.scale.bw
./bedGraphToBigWig riboBC18.minus.split.scale.bedgraph hg38.chrom.sizes riboBC18.minus.split.scale.bw

./bedGraphToBigWig riboBC19.plus.split.scale.bedgraph hg38.chrom.sizes riboBC19.plus.split.scale.bw
./bedGraphToBigWig riboBC19.minus.split.scale.bedgraph hg38.chrom.sizes riboBC19.minus.split.scale.bw
