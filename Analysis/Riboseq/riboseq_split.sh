#!/bin/bash
set -e

perl /data/fitzsimmonscm/scripts/icSHAPE-master/scripts/splitFastq.pl \
-U riboseq3_1.fastq \
-b 6:6 \
-l CGTGAT:1::ACATCG:2::GCCTAA:3::TGGTCA:4::CACTGT:5::ATTGGC:6::GATCTG:7::TCAAGT:8::CTGAGC:9::AAGCTA:10::GTAGCC:11::TACAAG:12::TTGACT:13::GGAACT:14::TGACAT:15::GGACGG:16::GCGGAC:17::TTTCAC:18::GGCCAC:19::CGAAAC:20 \
-d riboseq31_split
