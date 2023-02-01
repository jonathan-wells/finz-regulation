#!/usr/bin/env bash

# The following script loosely follows the workflow described in the Drop-seq Core Computational
# Protocol (Drop-seq cookbook), available at https://github.com/broadinstitute/Drop-seq. 

# Build STAR index
/programs/STAR-2.7.10a/STAR \
    --genomeDir "../data/expression/STARgenome" \
    --genomeFastaFiles "../data/expression/STARgenome/GCF_000002035.6_GRCz11_genomic.nonalt.fna" \
    --runMode genomeGenerate \
    --runThreadN 10 \
    --sjdbGTFfile "../data/expression/STARgenome/danRer11.nonalt.genes_tes.gtf" \
    --sjdbOverhang 99 

# Build dict for picard-tools
java -jar /programs/picard-tools-2.9.0/picard.jar CreateSequenceDictionary \
    R="../data/expression/STARgenome/GCF_000002035.6_GRCz11_genomic.nonalt.fna" \
    O="../data/expression/STARgenome/GCF_000002035.6_GRCz11_genomic.nonalt.dict"

for readfile in $(cat ../data/expression/farrell2018-bam/bam_stages.txt); do
    echo Processing $readfile ...
    
    # Extract reads from aligned bam file.
    java -jar /programs/picard-tools-2.9.0/picard.jar SamToFastq \
        INPUT="../data/expression/farrell2018-bam/${readfile}.bam" \
        FASTQ="../data/expression/farrell2018-fastq/${readfile}.fastq"

    # Align reads with STAR
    /programs/STAR-2.7.10a/STAR \
        --genomeDir "../data/expression/STARgenome" \
        --readFilesIn "../data/expression/farrell2018-fastq/${readfile}.fastq" \
        --runThreadN 10 \
        --outSAMtype BAM Unsorted \
        --runMode alignReads \
        --outFilterMultimapNmax 250 \
        --winAnchorMultimapNmax 500 \
        --outMultimapperOrder Random \
        --alignIntronMax 500000 \
        --alignMatesGapMax 500000 \
        --outFileNamePrefix "../data/expression/STARaligned/${readfile}.re"

    # Create unaligned bam file from original. The important parameters are
    # REMOVE_ALIGNMENT_INFORMATION and ATTRIBUTE_TO_CLEAR. Note that XM and XC are not cleared, as
    # these give the molecular barcode and cell barcode, respectively. These are later added back to
    # the STAR re-aligned bamfile.
    java -jar /programs/picard-tools-2.9.0/picard.jar RevertSam \
        I="../data/expression/farrell2018-bam/${readfile}.bam" \
        O="../data/expression/farrell2018-bam/${readfile}.unaligned.bam" \
        SANITIZE=true \
        MAX_DISCARD_FRACTION=0.005 \
        ATTRIBUTE_TO_CLEAR=XT \
        ATTRIBUTE_TO_CLEAR=GE \
        ATTRIBUTE_TO_CLEAR=GS \
        ATTRIBUTE_TO_CLEAR=XF \
        ATTRIBUTE_TO_CLEAR=ZP \
        ATTRIBUTE_TO_CLEAR=XN \
        ATTRIBUTE_TO_CLEAR=AS \
        ATTRIBUTE_TO_CLEAR=OC \
        ATTRIBUTE_TO_CLEAR=OP \
        SORT_ORDER=queryname \
        RESTORE_ORIGINAL_QUALITIES=true \
        REMOVE_DUPLICATE_INFORMATION=true \
        REMOVE_ALIGNMENT_INFORMATION=true 

    # Merge the new STAR alignment bam file with the unaligned bam file.
    java -jar /programs/picard-tools-2.9.0/picard.jar MergeBamAlignment \
        REFERENCE_SEQUENCE="../data/expression/STARgenome/GCF_000002035.6_GRCz11_genomic.nonalt.fna" \
        UNMAPPED_BAM="../data/expression/farrell2018-bam/${readfile}.unaligned.bam" \
        ALIGNED_BAM="../data/expression/STARaligned/${readfile}.reAligned.out.bam"\
        OUTPUT="../data/expression/STARaligned/${readfile}.realigned.barcoded.bam" \
        INCLUDE_SECONDARY_ALIGNMENTS=false \
        PAIRED_RUN=false \
        VALIDATION_STRINGENCY=SILENT
    
    # Cleaning up
    rm "../data/expression/STARaligned/${readfile}.reAligned.out.bam"
    rm "../data/expression/farrell2018-bam/${readfile}.unaligned.bam"
done
