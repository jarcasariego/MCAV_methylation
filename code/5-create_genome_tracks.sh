#!/bin/bash

#determine annotation types and display counts

grep -v '^#' Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.coding.gff3 | cut -s -f 3 | sort | uniq -c | sort -rn > all_features.txt

#extract feature types and generate individual gff

grep $'\tmRNA\t' Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.coding.gff3 > Mcav.GFFannotation.mRNA.gff
grep $'\tgene\t' Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.coding.gff3 > Mcav.GFFannotation.gene.gff
grep $'\texon\t' Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.coding.gff3 > Mcav.GFFannotation.exon.gff
grep $'\tCDS\t' Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.coding.gff3 > Mcav.GFFannotation.CDS.gff
grep $'\tfive_prime_UTR\t' Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.coding.gff3 > Mcav.GFFannotation.5-UTR.gff
grep $'\tthree_prime_UTR\t' Mcav_genome/Mcavernosa_annotation/Mcavernosa.maker.coding.gff3 > Mcav.GFFannotation.3-UTR.gff

##### CREATE OTHER GENOME TRACKS

# extract scaffold lenghts

cat Mcav_genome/Mcavernosa_July2018.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Mcav.Chromosome_lenghts.txt

# extract scaffold names  	

cut -f1 Mcav.Chromosome_lenghts.txt > Mcav.Chromosome-Names.txt

#Sort GFF files for downstream use

sortBed -faidx Mcav.Chromosome-Names.txt -i Mcav.GFFannotation.gene.gff > Mcav.GFFannotation.gene_sorted.gff

sortBed -faidx Mcav.Chromosome-Names.txt -i Mcav.GFFannotation.exon.gff > Mcav.GFFannotation.exon_sorted.gff

sortBed -faidx Mcav.Chromosome-Names.txt -i Mcav.GFFannotation.CDS.gff > Mcav.GFFannotation.CDS_sorted.gff

# Intergenic regions (By definition, these are regions that aren't genes. I can use complementBed to find all regions that aren't genes, and subtractBed to remove exons and create this track)

complementBed -i Mcav.GFFannotation.gene_sorted.gff -sorted -g Mcav.Chromosome_lenghts.txt | subtractBed -a - -b Mcav.GFFannotation.exon_sorted.gff > Mcav.GFFannotation.intergenic.gff # track resulting here has an overlap of the first base with the last base of the gene track so I corrected it below

awk '{print $1"\t"$2+1"\t"$3}' Mcav.GFFannotation.intergenic.gff > Mcav.GFFannotation.intergenic_corrected.gff #additionally the start region will be changed from 0 to 1 so need to be corrected also.
sed -i 's/\<1\>/0/g' Mcav.GFFannotation.intergenic_corrected.gff

#Non-coding Sequences (I can use complementBed to create a non-coding sequence track. This track can then be used to create an intron track)

complementBed -i  Mcav.GFFannotation.exon_sorted.gff -g Mcav.Chromosome_lenghts.txt > Mcav.GFFannotation.noncoding.gff

# Introns (The intersections betwen the non-coding sequences and genes are by definition introns)

intersectBed -a Mcav.GFFannotation.noncoding.gff3 -b Mcav.GFFannotation.gene_sorted.gff -sorted > Mcav.GFFannotation.intron.gff3 # track resulting here has an overlap of the first base with the last base of the exon track so I corrected it below

awk '{print $1"\t"$2+1"\t"$3}' Mcav.GFFannotation.intron.gff3 > Mcav.GFFannotation.intron_corrected.gff


# Putative promoter track

samtools faidx Mcav_genome/Mcavernosa_July2018.fasta #index genome

flankBed -i Mcav.GFFannotation.gene_sorted.gff -g Mcav_genome/Mcavernosa_July2018.fasta.fai -l 1000 -r 0 -s | awk '{ gsub("gene","put_promoter",$3); print $0 }'| awk '{if($5-$4 > 3)print $0}'| tr ' ' '\t' > Mcav.GFFannotation.put_promoter.gff
subtractBed -a Mcav.GFFannotation.put_promoter.gff -b Mcav.GFFannotation.gene_sorted.gff > Mcav.GFFannotation.put_promoter.gff # when genes are close to each other promoter region overlaps with the gene. 

##### Create repetitive region tracks

# Run repeat masker and use gff output > Mcav_TE_all.gff

RepeatMasker Mcav_genome/Mcavernosa_July2018.fasta -species "all" -par 8 -gff -excln 1> stdout.txt 2> stderr.txt  	

#### Create Concatenated Gff files

cat \
Mcav.GFFannotation.CDS_sorted.gff \
Mcav.GFFannotation.intron_corrected.gff \
Mcav.GFFannotation.5-UTR.gff \
Mcav.GFFannotation.3-UTR.gff \
Mcav.GFFannotation.put_promoter.gff \
Mcav.GFFannotation.intergenic_corrected.gff \
Mcav_TE_all.gff \
| awk -F"\t" '{if(($1 !~/gff-version/)&&($1 !~/Generated/)&&($1 !~/Project/))print $0}' \
| awk -F"\t" '{print $1"\t"$4"\t"$5"\t"$3}' \
| sortBed -i - \
| uniq \
| windowMaker \
-b - \
-w 2000 \
-i src \
| sortBed -i \
| uniq \
> Mcav.GFFannotation.concat_all2Kbin.corrected_uniq.bed
