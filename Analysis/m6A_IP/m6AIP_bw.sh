#!/bin/bash

#sbatch --mem=10g --cpus-per-task=8

./bedGraphToBigWig BC2.minus.split.scale.bedgraph hg38.chrom.sizes BC2.minus.split.scale.bw
./bedGraphToBigWig BC2.plus.split.scale.bedgraph hg38.chrom.sizes BC2.plus.split.scale.bw

./bedGraphToBigWig BC4.minus.split.scale.bedgraph hg38.chrom.sizes BC4.minus.split.scale.bw
./bedGraphToBigWig BC4.plus.split.scale.bedgraph hg38.chrom.sizes BC4.plus.split.scale.bw

./bedGraphToBigWig BC5.minus.split.scale.bedgraph hg38.chrom.sizes BC5.minus.split.scale.bw
./bedGraphToBigWig BC5.plus.split.scale.bedgraph hg38.chrom.sizes BC5.plus.split.scale.bw

./bedGraphToBigWig BC6.minus.split.scale.bedgraph hg38.chrom.sizes BC6.minus.split.scale.bw
./bedGraphToBigWig BC6.plus.split.scale.bedgraph hg38.chrom.sizes BC6.plus.split.scale.bw

./bedGraphToBigWig BC7.minus.split.scale.bedgraph hg38.chrom.sizes BC7.minus.split.scale.bw
./bedGraphToBigWig BC7.plus.split.scale.bedgraph hg38.chrom.sizes BC7.plus.split.scale.bw

./bedGraphToBigWig BC12.minus.split.scale.bedgraph hg38.chrom.sizes BC12.minus.split.scale.bw
./bedGraphToBigWig BC12.plus.split.scale.bedgraph hg38.chrom.sizes BC12.plus.split.scale.bw

./bedGraphToBigWig BC13.minus.split.scale.bedgraph hg38.chrom.sizes BC13.minus.split.scale.bw
./bedGraphToBigWig BC13.plus.split.scale.bedgraph hg38.chrom.sizes BC13.plus.split.scale.bw

./bedGraphToBigWig BC14.minus.split.scale.bedgraph hg38.chrom.sizes BC14.minus.split.scale.bw
./bedGraphToBigWig BC14.plus.split.scale.bedgraph hg38.chrom.sizes BC14.plus.split.scale.bw

./bedGraphToBigWig BC15.minus.split.scale.bedgraph hg38.chrom.sizes BC15.minus.split.scale.bw
./bedGraphToBigWig BC15.plus.split.scale.bedgraph hg38.chrom.sizes BC15.plus.split.scale.bw

./bedGraphToBigWig BC16.minus.split.scale.bedgraph hg38.chrom.sizes BC16.minus.split.scale.bw
./bedGraphToBigWig BC16.plus.split.scale.bedgraph hg38.chrom.sizes BC16.plus.split.scale.bw

./bedGraphToBigWig BC18.minus.split.scale.bedgraph hg38.chrom.sizes BC18.minus.split.scale.bw
./bedGraphToBigWig BC18.plus.split.scale.bedgraph hg38.chrom.sizes BC18.plus.split.scale.bw

./bedGraphToBigWig BC19.minus.split.scale.bedgraph hg38.chrom.sizes BC19.minus.split.scale.bw
./bedGraphToBigWig BC19.plus.split.scale.bedgraph hg38.chrom.sizes BC19.plus.split.scale.bw
