#!/usr/bin/env bash

# Extract DGE matrices from STAR-aligned reads using dropseq tools Version:2.5.1.
# Different number of cells removed for each sample/stage.

# ZFHIGH
echo "processing ZFHIGH"
DigitalExpression -m 60g \
--INPUT ../data/expression/STARaligned/ZFHIGH-DS5.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZFHIGH-DS5.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZFHIGH-DS5.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 1000 \
--MIN_NUM_READS_PER_CELL 1500 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 60g \
--INPUT ../data/expression/STARaligned/ZFHIGH-DS5b.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZFHIGH-DS5b.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZFHIGH-DS5b.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 1000 \
--MIN_NUM_READS_PER_CELL 1500 \
--STRAND_STRATEGY SENSE

wait

# ZFOBLONG
echo "processing ZFOBLONG"
DigitalExpression -m 60g \
--INPUT ../data/expression/STARaligned/ZFOBLONG-DS5.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZFOBLONG-DS5.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZFOBLONG-DS5.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 625 \
--MIN_NUM_READS_PER_CELL 1500 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 60g \
--INPUT ../data/expression/STARaligned/ZFOBLONG-DS5b.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZFOBLONG-DS5b.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZFOBLONG-DS5b.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 625 \
--MIN_NUM_READS_PER_CELL 1500 \
--STRAND_STRATEGY SENSE

wait

# ZFDOME
echo "processing ZFDOME"
DigitalExpression -m 120g \
--INPUT ../data/expression/STARaligned/ZFDOME-DS5.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZFDOME-DS5.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZFDOME-DS5.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 800 \
--MIN_NUM_READS_PER_CELL 2000 \
--STRAND_STRATEGY SENSE

# ZF30
echo "processing ZF30"
DigitalExpression -m 60g \
--INPUT ../data/expression/STARaligned/ZF30-DS5.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF30-DS5.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF30-DS5.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 625 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 60g \
--INPUT ../data/expression/STARaligned/ZF30-DS5b.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF30-DS5b.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF30-DS5b.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 625 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE

wait

# ZF50
echo "processing ZF50"
DigitalExpression -m 30g \
--INPUT ../data/expression/STARaligned/ZF50-DS2.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF50-DS2.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF50-DS2.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 600 \
--MIN_NUM_READS_PER_CELL 1500 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 30g \
--INPUT ../data/expression/STARaligned/ZF50-DS3.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF50-DS3.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF50-DS3.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 600 \
--MIN_NUM_READS_PER_CELL 1500 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 30g \
--INPUT ../data/expression/STARaligned/ZF50-DS4.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF50-DS4.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF50-DS4.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 600 \
--MIN_NUM_READS_PER_CELL 1500 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 30g \
--INPUT ../data/expression/STARaligned/ZF50-DS4b.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF50-DS4b.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF50-DS4b.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 600 \
--MIN_NUM_READS_PER_CELL 1500 \
--STRAND_STRATEGY SENSE

wait

# ZFS
echo "processing ZFS"
DigitalExpression -m 120g \
--INPUT ../data/expression/STARaligned/ZFS-DS5.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZFS-DS5.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZFS-DS5.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 600 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE

# ZF60
echo "processing ZF60"
DigitalExpression -m 40g \
--INPUT ../data/expression/STARaligned/ZF60-DS2.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF60-DS2.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF60-DS2.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 600 \
--MIN_NUM_READS_PER_CELL 1500 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 40g \
--INPUT ../data/expression/STARaligned/ZF60-DS3.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF60-DS3.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF60-DS3.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 600 \
--MIN_NUM_READS_PER_CELL 1500 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 40g \
--INPUT ../data/expression/STARaligned/ZF60-DS4.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF60-DS4.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF60-DS4.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 600 \
--MIN_NUM_READS_PER_CELL 1500 \
--STRAND_STRATEGY SENSE

wait
 
# ZF75
echo "processing ZF75"
DigitalExpression -m 40g \
--INPUT ../data/expression/STARaligned/ZF75-DS2.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF75-DS2.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF75-DS2.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 600 \
--MIN_NUM_READS_PER_CELL 1400 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 40g \
--INPUT ../data/expression/STARaligned/ZF75-DS3.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF75-DS3.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF75-DS3.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 600 \
--MIN_NUM_READS_PER_CELL 1400 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 40g \
--INPUT ../data/expression/STARaligned/ZF75-DS4.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF75-DS4.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF75-DS4.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 600 \
--MIN_NUM_READS_PER_CELL 1400 \
--STRAND_STRATEGY SENSE

wait

# ZF90
echo "processing ZF90"
DigitalExpression -m 40g \
--INPUT ../data/expression/STARaligned/ZF90-DS2.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF90-DS2.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF90-DS2.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 500 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 40g \
--INPUT ../data/expression/STARaligned/ZF90-DS3.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF90-DS3.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF90-DS3.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 500 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 40g \
--INPUT ../data/expression/STARaligned/ZF90-DS4.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF90-DS4.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF90-DS4.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 500 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE

wait

# ZFB
echo "processing ZFB"
DigitalExpression -m 30g \
--INPUT ../data/expr:ession/STARaligned/ZFB-DS2.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZFB-DS2.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZFB-DS2.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 500 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE & 

DigitalExpression -m 30g \
--INPUT ../data/expression/STARaligned/ZFB-DS2b.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZFB-DS2b.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZFB-DS2b.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 500 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 30g \
--INPUT ../data/expression/STARaligned/ZFB-DS3.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZFB-DS3.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZFB-DS3.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 500 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 30g \
--INPUT ../data/expression/STARaligned/ZFB-DS4.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZFB-DS4.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZFB-DS4.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 500 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE

wait

# ZF3S
echo "processing ZF3S"
DigitalExpression -m 120g \
--INPUT ../data/expression/STARaligned/ZF3S-DS5.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF3S-DS5.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF3S-DS5.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 500 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE &

# ZF6S
echo "processing ZF6S"
DigitalExpression -m 60g \
--INPUT ../data/expression/STARaligned/ZF6S-DS5.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF6S-DS5.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF6S-DS5.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 500 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE &

DigitalExpression -m 60g \
--INPUT ../data/expression/STARaligned/ZF6S-DS5b.final.bam \
--OUTPUT ../data/expression/dropseq-dge/ZF6S-DS5b.dge.txt.gz \
--SUMMARY ../data/expression/dropseq-dge/ZF6S-DS5b.summary.txt \
--TMP_DIR /workdir/jnw72/tmp \
--READ_MQ -1 \
--MIN_NUM_GENES_PER_CELL 500 \
--MIN_NUM_READS_PER_CELL 1000 \
--STRAND_STRATEGY SENSE

echo "All finished"
