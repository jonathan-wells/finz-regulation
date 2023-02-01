#!/usr/bin/env bash

# Extracts DGE matrix from STAR-aligned reads.
for acc in $(cat ../data/expression/farrell2018-bam/bam_stages.txt); do
DigitalExpression \
    -m 120g \
    --INPUT ../data/expression/STARaligned/${acc}.final.bam \
    --OUTPUT ../data/expression/dropseq-dge/${acc}.dge.txt.gz \
    --SUMMARY ../data/expression/dropseq-dge/${acc}.summary.txt \
    --READ_MQ 1 \
    --MIN_NUM_GENES_PER_CELL 625 \
    --MIN_NUM_READS_PER_CELL 1000 \
    --STRAND_STRATEGY SENSE
done
