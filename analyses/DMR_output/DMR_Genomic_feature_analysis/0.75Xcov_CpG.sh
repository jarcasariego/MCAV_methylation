%%bash

set -ex

#first create variables for each file grouping
declare -a symbC_control="/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CC-20-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CC-22-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CC-22-2_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CC-26-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CC-26-2_5x.tab"

declare -a symbC_heated="/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CH-20-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CH-20-2_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CH-20-3_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CH-22-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CH-22-2_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CH-26-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/CH-26-2_5x.tab"

declare -a symbD_heated="/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BH-20-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BH-20-2_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BH-22-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BH-22-2_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BH-26-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BH-26-2_5x.tab"

declare -a symbD_control="/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BC-20-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BC-20-2_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BC-22-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BC-22-2_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BC-26-1_5x.tab
/Users/jarcasariego/Desktop/2020-11-11-Genomic_feature_analysis/_5x.tab_files/BC-26-2_5x.tab"


#next create another list of variable names

groups=(symbC_control symbD_control symbC_heated symbD_heated)

#next loop through variable names to combine the files within each experimental group
#cat all group files together
#if the context is CG print the chromosome and position
#sort data
#unique and count the frequency of positions
#if the position is overlapping in 3/4 samples keep


for f in "${groups[@]}"
do
x="${f}" 
cat ${!x} | \
awk '{print $1,$2,$3}' | \
sort | \
uniq -c | \
awk '{if($1>2)print $2"\t"$3"\t"$4}' \
> ${f}.3xCpG.bed
done
