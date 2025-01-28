############################################################
###
###   MB rxns to host RNA correlation - stratified by AgeGroup
###
############################################################

library(RColorBrewer)
library(DESeq2)
library(gplots)
library(clusterProfiler)

source("customFunctions.R")

setwd("combineDataLayers/")


#######################################
## load matching metadata
df_brainMetaData <- as.data.frame( readRDS(file = "../mouseMetagenome/df_metaData.rds" ))
rownames(df_brainMetaData) <- df_brainMetaData$SampleName


#######################################
## Load metabolic modeling derived reaction abundances (microbiome side - metagenomics based)
df_brainMMrxn <- readRDS("df_rxnAbundancesMM.rds")
# relabel samples from SeqID to SampleName
df_brainMMrxn <- df_brainMMrxn[,df_brainMetaData$SequencingID]
colnames(df_brainMMrxn) <- df_brainMetaData$SampleName
# data is already normalized
colSums(df_brainMMrxn)


######################
## load host side transcriptomics data
obj_DifAbundBrainAllAges <- readRDS(file ="../mouseTranscriptome/obj_DifAbundBrainAllAges.rds")
df_brainRNAnorm <- as.data.frame( getVarianceStabilizedData(obj_DifAbundBrainAllAges) )
df_brainAgeDependance <- as.data.frame( results(obj_DifAbundBrainAllAges, name = "age", alpha = 0.05) )
rm(obj_DifAbundBrainAllAges)
# same order as metadata-table
df_brainRNAnorm <- df_brainRNAnorm[,df_brainMetaData$SampleName]


#########################################################
## check all three tables for matching sample identifiers
all.equal(colnames(df_brainRNAnorm), colnames(df_brainMMrxn), df_brainMetaData$SampleName)
# TRUE


############################
### For each Age-Group:
for(int_age in unique(df_brainMetaData$Age)) {
  print(int_age) 
  # remove rxn features with 0-variance (= equal values across all samples)
  if(!isEmpty(caret::nearZeroVar(x = t(df_brainMMrxn), saveMetrics = F))) {
    df_tmpBrainMMrxnFiltered <- df_brainMMrxn[-caret::nearZeroVar(x = t(df_brainMMrxn[,df_brainMetaData[df_brainMetaData$Age == int_age,"SampleName"]]), saveMetrics = F),
                                              df_brainMetaData[df_brainMetaData$Age == int_age,"SampleName"]]
  }
  # nothing to discard
  else df_tmpBrainMMrxnFiltered <- df_brainMMrxn[,df_brainMetaData[df_brainMetaData$Age == int_age,"SampleName"]]
  # remove transcript features with 0-variance (= equal values across all samples)
  df_tmpBrainRNAnormFiltered <- df_brainRNAnorm[-caret::nearZeroVar(x = t(df_brainRNAnorm[,df_brainMetaData[df_brainMetaData$Age == int_age,"SampleName"]]), saveMetrics = F),
                                                df_brainMetaData[df_brainMetaData$Age == int_age,"SampleName"]]
  # same colnames? (matching samples in both tables)
  if(!all.equal(colnames(df_tmpBrainMMrxnFiltered), colnames(df_tmpBrainRNAnormFiltered))) break;
  
  if(length(unique(df_brainMetaData[df_brainMetaData$Age == int_age,"Batch"])) < 2) {
    fct_batch <- NULL
  } else {
    fct_batch <- df_brainMetaData[df_brainMetaData$Age == int_age,"Batch"]
  }
  
  # run correlation all 2 all without partial cor.
  df_tmpCorRxn2BrainRNA <- df_pcorAll2All(df_table1 = df_tmpBrainMMrxnFiltered,
                                          df_table2 = df_tmpBrainRNAnormFiltered,
                                          var_pcorZ = fct_batch,
                                          str_labels = c("Rxn","Gene"),
                                          int_cores = 15) #method = "spearman"
  saveRDS(df_tmpCorRxn2BrainRNA, paste0("df_", int_age, "monthCorRxn2BrainRNA.rds"))
  
  # how many do we actually have?
  print( table(df_tmpCorRxn2BrainRNA$p.adj <= 0.05) )
  #
  print( table( df_tmpCorRxn2BrainRNA$p.adj <= 0.1 &
                abs(df_tmpCorRxn2BrainRNA$Spearman_rho) >= 0.55 ))
  
  rm(df_tmpCorRxn2BrainRNA, df_tmpBrainMMrxnFiltered, df_tmpBrainRNAnormFiltered)
  gc()
}
# [1] 24
# [1] "Starting correlation of 34089 Genes to 2129 Rxns within 10 samples! Wed Feb 28 05:44:32 2024"
# [1] "Finished all correlations! Preparing results. Wed Feb 28 06:06:46 2024"
# 
# FALSE     TRUE 
# 72564107    11374 
# 
# FALSE     TRUE 
# 72563987    11494 
# [1] 2
# [1] "Starting correlation of 33557 Genes to 2129 Rxns within 10 samples! Wed Feb 28 06:19:12 2024"
# [1] "Finished all correlations! Preparing results. Wed Feb 28 06:41:18 2024"
# 
# FALSE     TRUE 
# 71426487    16366 
# 
# FALSE     TRUE 
# 71426370    16483 
# [1] 15
# [1] "Starting correlation of 33776 Genes to 2129 Rxns within 10 samples! Wed Feb 28 06:53:12 2024"
# [1] "Finished all correlations! Preparing results. Wed Feb 28 07:15:16 2024"
# 
# FALSE     TRUE 
# 71898943    10161 
# 
# FALSE     TRUE 
# 71898791    10313 
# [1] 9
# [1] "Partial correlation enabled!"
# [1] "Starting partial correlation of 33404 Genes to 2129 Rxns within 10 samples! Wed Feb 28 07:26:56 2024"
# [1] "Finished all correlations! Preparing results. Wed Feb 28 13:41:04 2024"
# 
# FALSE     TRUE 
# 70365274   751839 
# 
# FALSE     TRUE 
# 69994633  1122480 
###
# [1] 30
# [1] "Partial correlation enabled!"
# [1] "Starting partial correlation of 34374 Genes to 2129 Rxns within 12 samples! Wed Feb 28 13:54:35 2024"
# [1] "Finished all correlations! Preparing results. Wed Feb 28 20:37:45 2024"
# 
# FALSE     TRUE 
# 72891797   290448 
# 
# FALSE     TRUE 
# 72719835   462410 

dbl_percentCorRxn2Meta <- list("2" = 16483/(71426370+16483)*100, 
                               "15" = 10313/(71898791+10313)*100, 
                               "24" = 11494/(72563987+11494)*100,
                               # pcor
                               "9" = 1122480/(69994633+1122480)*100, 
                               "30" = 462410/(72719835+462410)*100 )
pdf(file = "brain_hostMBcorrelationAgeStratified.pdf", width = 185/72, height = 185/72, pointsize = 10)
par(mar = c(3,5,3,1)+.1)
layout(matrix(c(1,1,1,2,2), nrow = 1))
barplot(height = as.numeric(dbl_percentCorRxn2Meta)[1:3],
        # names.arg = names(dbl_percentCorRxn2Meta),
        las = 2, ylim = c(0,0.026),
        main = "cor Host-MB", #"Correlations of Microbiome Reactions\n to Host Transcripts",
        col = brewer.pal("Pastel1", n = 5)[c(1,3,4)])
title(ylab = "Strong correlations [%]", line=4)
axis(side = 1, at = seq(from = 0.65, by = 1.225, length.out = 3),
     labels = names(dbl_percentCorRxn2Meta)[1:3], tick = F, line = F)
lines(x = c(0,1)+.75, y = rep(0.024,2)); text(x= 0.5+0.75, y= 0.025, labels = "***", cex = .8)
lines(x = c(1,2)+1, y = rep(0.019,2)); text(x= 1.5+1, y= 0.02, labels = "***", cex = .8)
# pcors
par(mar = c(3,2,3,3)+.1)
barplot(height = as.numeric(dbl_percentCorRxn2Meta)[4:5],
        # names.arg = names(dbl_percentCorRxn2Meta),
        las = 2, ylim = c(0,2.6),
        main = "pcor seq. batch", #"Partial Correlations of Microbiome Reactions\n to Host Transcripts corrected for Seq. Batch",
        col = brewer.pal("Pastel1", n = 5)[c(2,5)])
# title(ylab = "Strong correlations [%]", line=4)
axis(side = 1, at = seq(from = 0.65, by = 1.225, length.out = 2),
     labels = names(dbl_percentCorRxn2Meta)[4:5], tick = F, line = F)
lines(x = c(0,1)+.75, y = rep(1.9,2)); text(x= 0.5+0.75, y= 2, labels = "***", cex = .8)
dev.off()
# export at 440x400 barplot_strongCorHost2MetaLiverByAge


###########
# Pearson's chi-squared test:
# Null hypothesis: Joint distribution of the cell counts in the contingency table is the product of row and column sums (marginals)
####
# 2 to 15:
tbl_microbiomeCorrelation <- as.table( rbind(c(16483, 71426370), 
                                             c(10313, 71898791) ))
dimnames(tbl_microbiomeCorrelation) <- list(Age = c("2months", "15months"), 
                                            MicrobiomeCorrelation = c("strongMBcor", "weakMBcor"))
tbl_microbiomeCorrelation
#             MicrobiomeCorrelation
# Age        strongMBcor weakMBcor
# 2months          11248  56769182
# 15months          9300  57064932
####
# 15 to 24:
tbl_microbiomeCorrelation <- as.table( rbind(c(10313, 71898791),
                                             c(11494, 72563987) ))
dimnames(tbl_microbiomeCorrelation) <- list(Age = c("15months", "24months"), 
                                            MicrobiomeCorrelation = c("strongMBcor", "weakMBcor"))
####
# 9 to 30:
tbl_microbiomeCorrelation <- as.table( rbind(c(1122480, 69994633),
                                             c(462410, 72719835) ))
dimnames(tbl_microbiomeCorrelation) <- list(Age = c("9months", "30months"), 
                                            MicrobiomeCorrelation = c("strongMBcor", "weakMBcor"))

####
# Pearson's Chi-squared test with Yates' continuity correction
obj_chiSqMBcor <- chisq.test(x = tbl_microbiomeCorrelation)
obj_chiSqMBcor
# data:  2 months to 15 months:
# X-squared = 1460.9, df = 1, p-value < 2.2e-16
# data:  15 to 24 months
# X-squared = 53.44, df = 1, p-value = 2.667e-13
# data:  9 to 30 months
# X-squared = 297448, df = 1, p-value < 2.2e-16
obj_chiSqMBcor$expected     # expected counts under the null
obj_chiSqMBcor$residuals    # Pearson residuals  (positive numbers are overrepresented, negative numbers are underrepresented compared to null hpyothesis)

p2star(x = p.adjust(p = c(2.2e-16, 2.667e-13, 2.2e-16), method = "bonferroni") ) # "BH" 
# [1] "***" "***" "***"

