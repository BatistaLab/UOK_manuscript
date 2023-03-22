#!/bin/bash
set -e

module load python/2.7
python /data/fitzsimmonscm/extract-transcript-regions-master/extract_transcript_regions.py \
-i /data/BatistaLab_NGS/annotations/v32/gencode.v32.annotation.gtf \
-o SLAMbeds_v32 \
--gtf
