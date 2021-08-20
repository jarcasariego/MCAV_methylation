
#### Input ####
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse")
# Take p values from the limma-voom planned comparisons.

# Dataframe of mean methylation for each gene
dnam_gSummary <- readRDS("data/DNAm/CG_unstranded_5_geneLevelSummary_allIntragenicRegions.RData")
dnam_gSummary_5 <- dnam_gSummary[dnam_gSummary$n >= 5,]
# Subset to look at difference at Day 9
dnam_gSummary_5_D9diff <- dnam_gSummary_5$beta_9_2800_mean-dnam_gSummary_5$beta_9_400_mean
dnam_D9_final <- data.frame(gene_id=as.character(dnam_gSummary_5$geneID),diff=dnam_gSummary_5_D9diff)
# Subset to look at difference at Day 80
dnam_gSummary_5_D80diff <- dnam_gSummary_5$beta_80_2800_mean-dnam_gSummary_5$beta_80_400_mean
dnam_D80_final <- data.frame(gene_id=as.character(dnam_gSummary_5$geneID),diff=dnam_gSummary_5_D80diff)
# Save (already in github repo)
write.csv(dnam_D9_final,file = "data/Analysis/gomwu/GOMWU_CpG_D9.csv",row.names = FALSE)
write.csv(dnam_D80_final,file = "data/Analysis/gomwu/GOMWU_CpG_D80.csv",row.names = FALSE)
#### Performing the GO Mann-Whitney U (GOMWU) Test ####

## NOTE: the GO MWU functions assume all data and scripts are together in the current working directory.
# This may require some manual moving around. The simplest solution is creating a temperory folder that combines:
# the data and src gomwu folders and setting your wd directory to that path.

# Set working directory to folder with github repo
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/temp_goMWU/")
# Source GO MWU functions
source("gomwu.functions.R")

### Default settings and Go slim databases
## GO slim enrichment database
goDatabase="go.obo"
# download from http://www.geneontology.org/GO.downloads.ontology.shtml
## Enrichment categories
goDivision=c("MF","BP","CC") # Check all three divisions
## Annotation file
goAnnotations=paste0("GOMWU_Goterms.tab")

## Run log Bayes Factor for CpGs 
cpg_D9 <- list.files(pattern="CpG_D9.csv")
cpg_D80 <- list.files(pattern="CpG_D80.csv")
input <- c(cpg_D9,cpg_D80)

## Run log folder gene expression change  
ge_D9 <- list.files(pattern="LogFoldChangeGE_D9.tab")
ge_D80 <- list.files(pattern="LogFoldChangeGE_D80.tab")
input <- c(ge_D9,ge_D80)


## Running GO MWU function
# Loops through the two timepoints (runs them separately)
for(i in 1:length(input)){
  for(j in goDivision){
    gomwuStats(input[i], goDatabase, goAnnotations,j,
               perlPath="perl",
               largest=0.1,
               smallest=2,
               clusterCutHeight=0.25)
  }
}

## Options to create cladogram of significantly enriched genes 
# gomwuPlot(input[2],goAnnotations,goDivision[2],
#           absValue=2,  # genes with the measure value exceeding this will be counted as "good genes".
#           level1=0.1, # FDR threshold for plotting.
#           #Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
#           level2=0.05, # FDR cutoff to print in regular (not italic) font.
#           level3=0.01, # FDR cutoff to print in large bold font.
#           txtsize=1.2, # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
#           treeHeight=0.5, # height of the hierarchical clustering tree
#           #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
# )

### Read in all GO MWU outputs and compile into single table ###
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/data/Analysis/gomwu/outputs/")
# DNA methylation
D9_CpG_BP <- read.delim("MWU_BP_GOMWU_CpG_D9.csv",header = TRUE,sep=" ")
D9_CpG_BP <- data.frame(Seq="Methylation",Day="9",Category="Biological Process",D9_CpG_BP)
D80_CpG_BP <- read.delim("MWU_BP_GOMWU_CpG_D80.csv",header = TRUE,sep=" ")
D80_CpG_BP <- data.frame(Seq="Methylation",Day="80",Category="Biological Process",D80_CpG_BP)
D9_CpG_MF <- read.delim("MWU_MF_GOMWU_CpG_D9.csv",header = TRUE,sep=" ")
D9_CpG_MF <- data.frame(Seq="Methylation",Day="9",Category="Molecular Function",D9_CpG_MF)
D80_CpG_MF <- read.delim("MWU_MF_GOMWU_CpG_D80.csv",header = TRUE,sep=" ")
D80_CpG_MF <- data.frame(Seq="Methylation",Day="80",Category="Molecular Function",D80_CpG_MF)
D9_CpG_CC <- read.delim("MWU_CC_GOMWU_CpG_D9.csv",header = TRUE,sep=" ")
D9_CpG_CC <- data.frame(Seq="Methylation",Day="9",Category="Cellular Component",D9_CpG_CC)
D80_CpG_CC <- read.delim("MWU_CC_GOMWU_CpG_D80.csv",header = TRUE,sep=" ")
D80_CpG_CC <- data.frame(Seq="Methylation",Day="80",Category="Cellular Component",D80_CpG_CC)
# Gene Expression
D9_GE_BP <- read.delim("MWU_BP_GOMWU_LogFoldChangeGE_D9.tab",header = TRUE,sep=" ")
D9_GE_BP <- data.frame(Seq="Gene Expression",Day="9",Category="Biological Process",D9_GE_BP)
D80_GE_BP <- read.delim("MWU_BP_GOMWU_LogFoldChangeGE_D80.tab",header = TRUE,sep=" ")
D80_GE_BP <- data.frame(Seq="Gene Expression",Day="80",Category="Biological Process",D80_GE_BP)
D9_GE_MF <- read.delim("MWU_MF_GOMWU_LogFoldChangeGE_D9.tab",header = TRUE,sep=" ")
D9_GE_MF <- data.frame(Seq="Gene Expression",Day="9",Category="Molecular Function",D9_GE_MF)
D80_GE_MF <- read.delim("MWU_MF_GOMWU_LogFoldChangeGE_D80.tab",header = TRUE,sep=" ")
D80_GE_MF <- data.frame(Seq="Gene Expression",Day="80",Category="Molecular Function",D80_GE_MF)
D9_GE_CC <- read.delim("MWU_CC_GOMWU_LogFoldChangeGE_D9.tab",header = TRUE,sep=" ")
D9_GE_CC <- data.frame(Seq="Gene Expression",Day="9",Category="Cellular Component",D9_GE_CC)
D80_GE_CC <- read.delim("MWU_CC_GOMWU_LogFoldChangeGE_D80.tab",header = TRUE,sep=" ")
D80_GE_CC <- data.frame(Seq="Gene Expression",Day="80",Category="Cellular Component",D80_GE_CC)

#Combining everything (manually)
final_df <- rbind(D9_CpG_BP,D80_CpG_BP,D9_CpG_MF,D80_CpG_MF,D9_CpG_CC,
                  D80_CpG_CC,D9_GE_BP,D80_GE_BP,D9_GE_MF,D80_GE_MF,
                  D9_GE_CC,D80_GE_CC)
final_df <- final_df[final_df$p.adj <= 0.05,]
write.csv(final_df,"GOMWU_finalSummary.csv")


