
library(tidyr)
library(data.table)
library(gplots)
library(ggplot2) #install.packages('ggplot2')
library(dplyr)
library(broom)
library(RColorBrewer)
library(egg) #install.packages('egg')
library(purrr)
library(nlme)
library(vegan)

meta_data <- read.csv("../data/Treatment_metadata.csv", stringsAsFactors = FALSE)
meta_data <- meta_data[,c(3:11)]
colnames(meta_data) <- c("species", "colony", "core", "sample", "trt1", "trt2", "TreatComb", "Tank","Sample.ID")

expr_rep <- read.csv("../analyses/Repeat_analysis_MCAV/expression/expr_vs_treat_rep.csv")
colnames(expr_rep)[1] <- "repeat_type" 
expr_rep_stacked <- expr_rep %>% 
  group_by(repeat_type) %>%
  slice(1) %>%
  gather(key = "Sample.ID", value = "perc_expression", -repeat_type) 
expr_rep_stacked$Sample.ID <- gsub("\\.","-",expr_rep_stacked$Sample.ID)
expr_rep_stacked <- merge(expr_rep_stacked, meta_data, by = "Sample.ID")
#create a group variable for the combination of treatments
expr_rep_stacked <- expr_rep_stacked %>%
  mutate(sym = recode(trt1, b = "D", c = "C"),
         group = paste(sym, trt2, sep = ""))
expr_rep_stacked <- expr_rep_stacked[-7]

#arcsin sqrt transformation function
asinTransform <- function(p) { asin(sqrt(p))}

#arcsin transform data 
expr_rep_stacked_asin <- expr_rep_stacked

expr_rep_2way_aov <- expr_rep_stacked_asin %>% group_by(repeat_type) %>%
  do(expression_aov_models = aov(perc_expression ~sym*trt2, data = expr_rep_stacked_asin))
#summarize ANOVA data
expr_rep_2way_aov_modelsumm  <- expr_rep_2way_aov %>% ungroup %>% 
  pull(expression_aov_models) %>% 
  map_dfr(tidy, .id = 'grp')
#spread out pvalues
expr_rep_2way_aov_modelsumm_wide <- data.frame(tidyr::pivot_wider(expr_rep_2way_aov_modelsumm, names_from = term, values_from = c("df", "sumsq", "meansq","statistic","p.value"),names_sep = "_" ))
expr_rep_2way_aov_modelsumm_wide <- cbind(expr_rep_2way_aov[,1], expr_rep_2way_aov_modelsumm_wide[,-c(1)])
all_aov_sig_sym <- expr_rep_2way_aov_modelsumm_wide[which(expr_rep_2way_aov_modelsumm_wide$p.value_sym <= 0.1), 1]
all_aov_sig_temp <- expr_rep_2way_aov_modelsumm_wide[which(expr_rep_2way_aov_modelsumm_wide$p.value_trt2 <= 0.1), 1]
all_aov_sig_Inter <- expr_rep_2way_aov_modelsumm_wide[which(expr_rep_2way_aov_modelsumm_wide$p.value_sym.trt2 <= 0.1), 1]

expr_rep_1way_aov <- expr_rep_stacked_asin %>% group_by(repeat_type) %>%
  do(expression_aov_models = aov(perc_expression ~group, data = . ))
#summarize ANOVA data
expr_rep_1way_aov_modelsumm  <- expr_rep_1way_aov %>% ungroup %>% 
  pull(expression_aov_models) %>% 
  map_dfr(tidy, .id = 'grp')
#spread out pvalues
expr_rep_1way_aov_modelsumm_wide <- data.frame(tidyr::pivot_wider(expr_rep_1way_aov_modelsumm, names_from = term, values_from = c("df", "sumsq", "meansq","statistic","p.value"),names_sep = "_" ))
expr_rep_1way_aov_modelsumm_wide <- cbind(expr_rep_1way_aov[,1], expr_rep_1way_aov_modelsumm_wide[,-c(1)])
aov_sig_group <- expr_rep_1way_aov_modelsumm_wide[which(expr_rep_1way_aov_modelsumm_wide$p.value_gro <= 0.05), 1]

