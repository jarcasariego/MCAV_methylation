#!/bin/bash

# Script intersects CpG bed tracks with the genome track file, counting the number of CpGs from
# each track in each gene. It also calculates mean summary values for several CpG related stats
# for each gene (i.e. summarizes across CpGs within each gene)

#Path variables and date
saveDate="20210120"
outputDir="/scratch/jeirinlo/jrodr979/DML_analysis/"
CpGpath="/scratch/jeirinlo/jrodr979/DML_analysis/20210119_FeatureCpGDist/bed_files/all.colonies/"
genomeTrack="/scratch/jeirinlo/jrodr979/2021-01-08-Mcav_Genomic-feature-track_corrected/"

# Different CpG bed files
all="${CpGpath}CpG_allCpGs.bed"
cov5="${CpGpath}CpG_5xCov.bed"

# Genome track bed file
geneL=${genomeTrack}Mcav.GFFannotation.gene_sorted.gff
exonL=${genomeTrack}Mcav.GFFannotation.exon_sorted.gff

outputFile=${outputDir}${saveDate}_CpGbyGeneSummary
mkdir -p ${outputFile}

#### Intersect and Count all CpGs for each gene ####
echo "Intersect and Count all CpGs for each gene... "

${pathBedTools}bedtools intersect \
-wa -wb \
-a ${geneL} \
-b ${all} \
> ${outputFile}/${saveDate}_temp_allIntersect.txt

${pathBedTools}bedtools groupby \
-i ${outputFile}/${saveDate}_temp_allIntersect.txt \
-g 1,4,5 \
-c 13 \
-o count \
> ${outputFile}/${saveDate}_allCpG_Count.txt

#### Intersect and Count and Summarize Methylation for all CpGs with 5x coverage for each gene ####
echo "Counting up all CpGs with 5x coverage in genes ... "

${pathBedTools}bedtools intersect \
-wa -wb \
-a ${geneL} \
-b ${cov5} \
> ${outputFile}/${saveDate}_temp_cov5Intersect.txt

${pathBedTools}bedtools groupby \
-i ${outputFile}/${saveDate}_temp_cov5Intersect.txt \
-g 1,4,5 \
-c 13 \
-o count \
> ${outputFile}/${saveDate}_cov5CpG_Count.txt

# Summarize the methylation for CpGs with coverage 
${pathBedTools}bedtools groupby \
-i ${outputFile}/${saveDate}_temp_cov5Intersect.txt \
-g 1,4,5 \
-c 14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34 \
-o mean \
> ${outputFile}/${saveDate}_cov5CpG_Summarize.txt


# We are also going to generate a gene to exon intersect and count to count the
# number of exons in each gene.
echo "Counting all exons in each gene... "

${pathBedTools}bedtools intersect \
-wa -wb \
-a ${geneL} \
-b ${exonL} \
> ${outputFile}/${saveDate}_temp_exon_Intersect.txt

${pathBedTools}bedtools groupby \
-i ${outputFile}/${saveDate}_temp_exon_Intersect.txt \
-g 1,4,5 \
-c 13 \
-o count \
> ${outputFile}/${saveDate}_exonInGene_Count.txt

rm -rf ${outputFile}/*temp*

