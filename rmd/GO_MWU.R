
#Read methylation data

meth.data <- list.files(path = "../analyses/Gene_Methylation/bed/", pattern = ".bed$", full.names=TRUE) %>%
  set_names(.) %>% 
  map_dfr(read.csv,.id="Sample.ID", header=FALSE, sep="\t", na.string="NA", stringsAsFactors = FALSE) %>% 
  dplyr::select(-c(V3,V7:V14)) %>%
  group_by(Sample.ID)
colnames(meth.data) <- c("Sample.ID", "scaffold", "position","per.meth","meth","unmeth","gene")
meth.data$gene <- gsub(";.*","",meth.data$gene) #remove extra characters
meth.data$gene <- gsub("ID=","",meth.data$gene) #remove extra characters
meth.data$Sample.ID <- gsub("../analyses/Gene_Methylation/bed/","",meth.data$Sample.ID) #remove extra characters
meth.data$Sample.ID <- gsub("_.*","",meth.data$Sample.ID) #remove extra characters 
MD.All <- merge(meth.data, sample.info, by="Sample.ID")
MD.All <- MD.All %>%
  mutate(sym = recode(trt1, b = "D", c = "C"),
         treatment = paste(sym, trt2, sep = ""),
         beta = per.meth/100)

meth.means <- aggregate(per.meth ~ treatment*gene, data=MD.All, FUN=mean) %>%
  spread(., key = treatment, value = c(per.meth))

#all genes are covered with at least 5 CpG. NO NEED TO FILTER

dnam.gene_summ_5 <- meth.means 

# Subset to look at diff between contrasts 
##Ch_Cc
dnam.gene_summ_5_ChCcdiff <- dnam.gene_summ_5$Ch-dnam.gene_summ_5$Cc
dnam_SymbC_final <- data.frame(gene_id=as.character(dnam.gene_summ_5$gene),diff=dnam.gene_summ_5_ChCcdiff)
##Dc_Cc
dnam.gene_summ_5_DcCcdiff <- dnam.gene_summ_5$Dc-dnam.gene_summ_5$Cc
dnam_Ctrol_final <- data.frame(gene_id=as.character(dnam.gene_summ_5$gene),diff=dnam.gene_summ_5_DcCcdiff)
##Dh_Dc
dnam.gene_summ_5_DhDcdiff <- dnam.gene_summ_5$Dh-dnam.gene_summ_5$Dc
dnam_SymbD_final <- data.frame(gene_id=as.character(dnam.gene_summ_5$gene),diff=dnam.gene_summ_5_DhDcdiff)
##Dh_Ch
dnam.gene_summ_5_DhChdiff <- dnam.gene_summ_5$Dh-dnam.gene_summ_5$Ch
dnam_Heated_final <- data.frame(gene_id=as.character(dnam.gene_summ_5$gene),diff=dnam.gene_summ_5_DhChdiff)

write.csv(dnam_Ctrol_final,file = "../analyses/gomwu/GOMWU_CpG_DcCc.csv",row.names = FALSE)
write.csv(dnam_Heated_final,file = "../analyses/gomwu/GOMWU_CpG_DhCh.csv",row.names = FALSE)
write.csv(dnam_SymbC_final,file = "../analyses/gomwu/GOMWU_CpG_ChCc.csv",row.names = FALSE)
write.csv(dnam_SymbD_final,file = "../analyses/gomwu/GOMWU_CpG_DhDc.csv",row.names = FALSE)

# Subset and calculate log2fold for each contrasts 
##Ch_Cc
dnam.gene_summ_5_ChCclf2 <- log2(dnam.gene_summ_5$Ch/dnam.gene_summ_5$Cc)
dnam_SymbC_final_lf2 <- data.frame(gene_id=as.character(dnam.gene_summ_5$gene),log2fold=dnam.gene_summ_5_ChCclf2)
##Dc_Cc
dnam.gene_summ_5_DcCclf2 <- log2(dnam.gene_summ_5$Dc/dnam.gene_summ_5$Cc)
dnam_Ctrol_final_lf2 <- data.frame(gene_id=as.character(dnam.gene_summ_5$gene),log2fold=dnam.gene_summ_5_DcCclf2)
##Dh_Dc
dnam.gene_summ_5_DhDclf2 <- log2(dnam.gene_summ_5$Dh/dnam.gene_summ_5$Dc)
dnam_SymbD_final_lf2 <- data.frame(gene_id=as.character(dnam.gene_summ_5$gene),log2fold=dnam.gene_summ_5_DhDclf2)
##Dh_Ch
dnam.gene_summ_5_DhChlf2 <- log2(dnam.gene_summ_5$Dh/dnam.gene_summ_5$Ch)
dnam_Heated_final_lf2 <- data.frame(gene_id=as.character(dnam.gene_summ_5$gene),log2fold=dnam.gene_summ_5_DhChlf2)

write.csv(dnam_Ctrol_final_lf2,file = "../analyses/gomwu/GOMWU_CpG_lf2_DcCc.csv",row.names = FALSE)
write.csv(dnam_Heated_final_lf2,file = "../analyses/gomwu/GOMWU_CpG_lf2_DhCh.csv",row.names = FALSE)
write.csv(dnam_SymbC_final_lf2,file = "../analyses/gomwu/GOMWU_CpG_lf2_ChCc.csv",row.names = FALSE)
write.csv(dnam_SymbD_final_lf2,file = "../analyses/gomwu/GOMWU_CpG_lf2_DhDc.csv",row.names = FALSE)

## NOTE: the GO MWU functions assume all data and scripts are together in the current working directory.
# the data and src gomwu folders and setting your wd directory to that path.

library(ape)

# Set working directory to folder with github repo
setwd("../analyses/gomwu/")
# Source GO MWU functions
source("gomwu.functions.R")

goDatabase="go.obo"
# download from http://www.geneontology.org/GO.downloads.ontology.shtml
## Enrichment categories
goDivision=c("MF","BP","CC") # Check all three divisions
## Annotation file
goAnnotations=paste0("Mcavernosa_gene2go.tab")

## Run log Bayes Factor for CpGs diff 
cpg_Ctrol <- "GOMWU_CpG_DcCc.csv"
cpg_Heated <- "GOMWU_CpG_DhCh.csv"
cpg_SymbC <- "GOMWU_CpG_ChCc.csv"
cpg_SymbD <- "GOMWU_CpG_DhDc.csv"
input <- c(cpg_Ctrol,cpg_Heated,cpg_SymbC,cpg_SymbD)

## Run log Bayes Factor for CpGs log2fold 
cpg_lf2_Ctrol <- "GOMWU_CpG_lf2_DcCc.csv"
cpg_lf2_Heated <- "GOMWU_CpG_lf2_DhCh.csv"
cpg_lf2_SymbC <- "GOMWU_CpG_lf2_ChCc.csv"
cpg_lf2_SymbD <- "GOMWU_CpG_lf2_DhDc.csv"
input2 <- c(cpg_lf2_Ctrol,cpg_lf2_Heated,cpg_lf2_SymbC,cpg_lf2_SymbD)

## GO for diff
for(i in 1:length(input)){
  for(j in goDivision){
    gomwuStats(input[i], goDatabase, goAnnotations,j,
               perlPath="C:/Strawberry/perl/bin/perl.exe",
               largest=0.1,
               smallest=2,
               clusterCutHeight=0.25)
  }
}

# DNA methylation
DcCc_CpG_BP <- read.delim("MWU_BP_GOMWU_CpG_DcCc.csv",header = TRUE,sep=" ")
DcCc_CpG_BP <- data.frame(Seq="Methylation",Contrast="Dc.versus.Cc",Category="Biological Process",DcCc_CpG_BP)
DhCh_CpG_BP <- read.delim("MWU_BP_GOMWU_CpG_DhCh.csv",header = TRUE,sep=" ")
DhCh_CpG_BP <- data.frame(Seq="Methylation",Contrast="Dh.versus.Ch",Category="Biological Process",DhCh_CpG_BP)
ChCc_CpG_BP <- read.delim("MWU_BP_GOMWU_CpG_ChCc.csv",header = TRUE,sep=" ")
ChCc_CpG_BP <- data.frame(Seq="Methylation",Contrast="Ch.versus.Cc",Category="Biological Process",ChCc_CpG_BP)
DhDc_CpG_BP <- read.delim("MWU_BP_GOMWU_CpG_DhDc.csv",header = TRUE,sep=" ")
DhDc_CpG_BP <- data.frame(Seq="Methylation",Contrast="Dh.versus.Dc",Category="Biological Process",DhDc_CpG_BP)

DcCc_CpG_MF <- read.delim("MWU_MF_GOMWU_CpG_DcCc.csv",header = TRUE,sep=" ")
DcCc_CpG_MF <- data.frame(Seq="Methylation",Contrast="Dc.versus.Cc",Category="Biological Process",DcCc_CpG_MF)
DhCh_CpG_MF <- read.delim("MWU_MF_GOMWU_CpG_DhCh.csv",header = TRUE,sep=" ")
DhCh_CpG_MF <- data.frame(Seq="Methylation",Contrast="Dh.versus.Ch",Category="Biological Process",DhCh_CpG_MF)
ChCc_CpG_MF <- read.delim("MWU_MF_GOMWU_CpG_ChCc.csv",header = TRUE,sep=" ")
ChCc_CpG_MF <- data.frame(Seq="Methylation",Contrast="Ch.versus.Cc",Category="Biological Process",ChCc_CpG_MF)
DhDc_CpG_MF <- read.delim("MWU_MF_GOMWU_CpG_DhDc.csv",header = TRUE,sep=" ")
DhDc_CpG_MF <- data.frame(Seq="Methylation",Contrast="Dh.versus.Dc",Category="Biological Process",DhDc_CpG_MF)

DcCc_CpG_CC <- read.delim("MWU_CC_GOMWU_CpG_DcCc.csv",header = TRUE,sep=" ")
DcCc_CpG_CC <- data.frame(Seq="Methylation",Contrast="Dc.versus.Cc",Category="Biological Process",DcCc_CpG_CC)
DhCh_CpG_CC <- read.delim("MWU_CC_GOMWU_CpG_DhCh.csv",header = TRUE,sep=" ")
DhCh_CpG_CC <- data.frame(Seq="Methylation",Contrast="Dh.versus.Ch",Category="Biological Process",DhCh_CpG_CC)
ChCc_CpG_CC <- read.delim("MWU_CC_GOMWU_CpG_ChCc.csv",header = TRUE,sep=" ")
ChCc_CpG_CC <- data.frame(Seq="Methylation",Contrast="Ch.versus.Cc",Category="Biological Process",ChCc_CpG_CC)
DhDc_CpG_CC <- read.delim("MWU_CC_GOMWU_CpG_DhDc.csv",header = TRUE,sep=" ")
DhDc_CpG_CC <- data.frame(Seq="Methylation",Contrast="Dh.versus.Dc",Category="Biological Process",DhDc_CpG_CC)


#GO for L2F
for(i in 1:length(input2)){
  for(j in goDivision){
    gomwuStats(input2[i], goDatabase, goAnnotations,j,
               perlPath="C:/Strawberry/perl/bin/perl.exe",
               largest=0.1,
               smallest=2,
               clusterCutHeight=0.25)
  }
}

# DNA methylation
DcCc_CpG_lf2_BP <- read.delim("MWU_BP_GOMWU_CpG_lf2_DcCc.csv",header = TRUE,sep=" ")
DcCc_CpG_lf2_BP <- data.frame(Seq="Methylation",Contrast="Dc.versus.Cc",Category="Biological Process",DcCc_CpG_lf2_BP)
DhCh_CpG_lf2_BP <- read.delim("MWU_BP_GOMWU_CpG_lf2_DhCh.csv",header = TRUE,sep=" ")
DhCh_CpG_lf2_BP <- data.frame(Seq="Methylation",Contrast="Dh.versus.Ch",Category="Biological Process",DhCh_CpG_lf2_BP)
ChCc_CpG_lf2_BP <- read.delim("MWU_BP_GOMWU_CpG_lf2_ChCc.csv",header = TRUE,sep=" ")
ChCc_CpG_lf2_BP <- data.frame(Seq="Methylation",Contrast="Ch.versus.Cc",Category="Biological Process",ChCc_CpG_lf2_BP)
DhDc_CpG_lf2_BP <- read.delim("MWU_BP_GOMWU_CpG_lf2_DhDc.csv",header = TRUE,sep=" ")
DhDc_CpG_lf2_BP <- data.frame(Seq="Methylation",Contrast="Dh.versus.Dc",Category="Biological Process",DhDc_CpG_lf2_BP)

DcCc_CpG_lf2_MF <- read.delim("MWU_MF_GOMWU_CpG_lf2_DcCc.csv",header = TRUE,sep=" ")
DcCc_CpG_lf2_MF <- data.frame(Seq="Methylation",Contrast="Dc.versus.Cc",Category="Biological Process",DcCc_CpG_lf2_MF)
DhCh_CpG_lf2_MF <- read.delim("MWU_MF_GOMWU_CpG_lf2_DhCh.csv",header = TRUE,sep=" ")
DhCh_CpG_lf2_MF <- data.frame(Seq="Methylation",Contrast="Dh.versus.Ch",Category="Biological Process",DhCh_CpG_lf2_MF)
ChCc_CpG_lf2_MF <- read.delim("MWU_MF_GOMWU_CpG_lf2_ChCc.csv",header = TRUE,sep=" ")
ChCc_CpG_lf2_MF <- data.frame(Seq="Methylation",Contrast="Ch.versus.Cc",Category="Biological Process",ChCc_CpG_lf2_MF)
DhDc_CpG_lf2_MF <- read.delim("MWU_MF_GOMWU_CpG_lf2_DhDc.csv",header = TRUE,sep=" ")
DhDc_CpG_lf2_MF <- data.frame(Seq="Methylation",Contrast="Dh.versus.Dc",Category="Biological Process",DhDc_CpG_lf2_MF)

DcCc_CpG_lf2_CC <- read.delim("MWU_CC_GOMWU_CpG_lf2_DcCc.csv",header = TRUE,sep=" ")
DcCc_CpG_lf2_CC <- data.frame(Seq="Methylation",Contrast="Dc.versus.Cc",Category="Biological Process",DcCc_CpG_lf2_CC)
DhCh_CpG_lf2_CC <- read.delim("MWU_CC_GOMWU_CpG_lf2_DhCh.csv",header = TRUE,sep=" ")
DhCh_CpG_lf2_CC <- data.frame(Seq="Methylation",Contrast="Dh.versus.Ch",Category="Biological Process",DhCh_CpG_lf2_CC)
ChCc_CpG_lf2_CC <- read.delim("MWU_CC_GOMWU_CpG_lf2_ChCc.csv",header = TRUE,sep=" ")
ChCc_CpG_lf2_CC <- data.frame(Seq="Methylation",Contrast="Ch.versus.Cc",Category="Biological Process",ChCc_CpG_lf2_CC)
DhDc_CpG_lf2_CC <- read.delim("MWU_CC_GOMWU_CpG_lf2_DhDc.csv",header = TRUE,sep=" ")
DhDc_CpG_lf2_CC <- data.frame(Seq="Methylation",Contrast="Dh.versus.Dc",Category="Biological Process",DhDc_CpG_lf2_CC)


setwd("../../output/")

final_df <- rbind(DcCc_CpG_BP,DhCh_CpG_BP,ChCc_CpG_BP,DhDc_CpG_BP,DcCc_CpG_MF,DhCh_CpG_MF,
                  ChCc_CpG_MF,DhDc_CpG_MF,DcCc_CpG_CC,DhCh_CpG_CC,ChCc_CpG_CC,DhDc_CpG_CC)
final_df <- final_df[final_df$p.adj <= 0.05,]
write.csv(final_df,"GOMWU_finalSummary.csv")

final_lf2_df <- rbind(DcCc_CpG_lf2_BP,
                      DhCh_CpG_lf2_BP,
                      ChCc_CpG_lf2_BP,
                      DhDc_CpG_lf2_BP,
                      DcCc_CpG_lf2_MF,
                      DhCh_CpG_lf2_MF,
                      ChCc_CpG_lf2_MF,
                      DhDc_CpG_lf2_MF,
                      DcCc_CpG_lf2_CC,
                      DhCh_CpG_lf2_CC,
                      ChCc_CpG_lf2_CC,
                      DhDc_CpG_lf2_CC)
final_lf2_df <- final_lf2_df[final_lf2_df$p.adj <= 0.05,]



