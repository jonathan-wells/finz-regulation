#!/usr/bin/env bash

####################################################################################################
## Takes STAR-aligned BAM files and adds gene/TE annotations prior to generation of digital gene
# expression (DGE) matrix.
####################################################################################################

# Remove TEs overlapping gene transcripts.
rg "\ttranscript\t" ../data/expression/STARgenome/Danio_rerio.GRCz11.101.curated_finz.gtf > transcript.gtf

bedtools intersect \
    -a ../data/expression/STARgenome/danRer11.nonalt.tetranscripts.gtf \
    -b transcript.gtf \
    -wa | sort -u > overlapping_transcript.gtf

bedtools subtract \
    -a ../data/expression/STARgenome/danRer11.nonalt.genes_tes.gtf \
    -b overlapping_transcript.gtf \
    -A \
    -f 1.0 \
    -F 1.0 \
    -s > tmp.gtf 

# Stupid but neccessary hack to temporarily rename TEs (and genes in the process) so that
# TagReadWithGeneFunction doesn't choke on the fact that TEs from the same family are present on
# different chromosomes. Gets changed back in BAM file after tagging is complete.
awk \
    -F '"' \
    -v OFS='"' \
    '{ 
    if (index($2, "ENSDART") != 0 || $2 ~ /^g[0-9]*/) 
        print $1, $2, "; transcript_name \"" $2, $3,  $4, "; gene_name \"" $4, $5, $6, $7, $8,"" 
    else
        print $1, $2, "; transcript_name \"" $2, $3,  $4"_dup"NR, "; gene_name \"" $4"_dup"NR , $5, $6, $7, $8,"" 
    }' tmp.gtf |
    sed 's/;""""$//' > "../data/expression/STARgenome/danRer11.nonalt.genes_tes_notranscriptoverlap.gtf"

# Remove temporary gtf files
rm transcript.gtf
rm overlapping_transcript.gtf
rm tmp.gtf

# Tag reads using Drop-seq tools version:2.5.1
for acc in $(cat ../data/expression/farrell2018-bam/bam_stages.txt); do
    echo processing ${acc} ...
    # Tag reads with gene attributes
    TagReadWithGeneFunction \
        -m 90G \
        --INPUT "../data/expression/STARaligned/${acc}.realigned.barcoded.bam" \
        --OUTPUT tmp.bam \
        --ANNOTATIONS_FILE "../data/expression/STARgenome/danRer11.nonalt.genes_tes_notranscriptoverlap.gtf"

    # Clean up _dup additions to gene names in BAM file.
    samtools view -h tmp.bam |
        sed 's/_dup[0-9]*//g' |
        samtools view -h -b > "../data/expression/STARaligned/${acc}.final.bam"

    rm tmp.bam
done

# Clean up _dup additions to gene names in GTF file.
sed 's/_dup[0-9]*//g' "../data/expression/STARgenome/danRer11.nonalt.genes_tes_notranscriptoverlap.gtf" \
    > tmp
mv tmp "../data/expression/STARgenome/danRer11.nonalt.genes_tes_notranscriptoverlap.gtf"

