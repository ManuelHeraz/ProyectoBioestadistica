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
df_liverMetaData <- as.data.frame( readRDS(file = "df_liverMetaData.rds" ))
rownames(df_liverMetaData) <- df_liverMetaData$SampleName


#######################################
## Load metabolic modeling derived reaction abundances (microbiome side - metagenomics based)
df_liverMMrxn <- readRDS("df_rxnAbundancesMM.rds")
# relabel samples from SeqID to SampleName
df_liverMMrxn <- df_liverMMrxn[,df_liverMetaData$SequencingID]
colnames(df_liverMMrxn) <- df_liverMetaData$SampleName
# data is already normalized
colSums(df_liverMMrxn)


######################
## load host side transcriptomics data
obj_DifAbundLiverAllAges <- readRDS(file = "../mouseTranscriptome/obj_DifAbundLiverAllAges.rds")
df_liverRNAnorm <- as.data.frame( getVarianceStabilizedData(obj_DifAbundLiverAllAges) )
df_liverAgeDependance <- as.data.frame( results(obj_DifAbundLiverAllAges, name = "age", alpha = 0.05) )
rm(obj_DifAbundLiverAllAges)
# same order as metadata-table
df_liverRNAnorm <- df_liverRNAnorm[,df_liverMetaData$SampleName]


#########################################################
## check all three tables for matching sample identifiers
all.equal(colnames(df_liverRNAnorm), colnames(df_liverMMrxn), df_liverMetaData$SampleName)


############################
### For each Age-Group:
for(int_age in unique(df_liverMetaData$Age)) {
  print(int_age) 
  # remove rxn features with 0-variance (= equal values across all samples)
  if(!isEmpty(caret::nearZeroVar(x = t(df_liverMMrxn), saveMetrics = F))) {
    df_tmpLiverMMrxnFiltered <- df_liverMMrxn[-caret::nearZeroVar(x = t(df_liverMMrxn[,df_liverMetaData[df_liverMetaData$Age == int_age,"SampleName"]]), saveMetrics = F),
                                              df_liverMetaData[df_liverMetaData$Age == int_age,"SampleName"]]
  }
  # nothing to discard
  else df_tmpLiverMMrxnFiltered <- df_liverMMrxn[,df_liverMetaData[df_liverMetaData$Age == int_age,"SampleName"]]
  # remove transcript features with 0-variance (= equal values across all samples)
  df_tmpLiverRNAnormFiltered <- df_liverRNAnorm[-caret::nearZeroVar(x = t(df_liverRNAnorm[,df_liverMetaData[df_liverMetaData$Age == int_age,"SampleName"]]), saveMetrics = F),
                                                df_liverMetaData[df_liverMetaData$Age == int_age,"SampleName"]]
  # same colnames? (matching samples in both tables)
  if(!all.equal(colnames(df_tmpLiverMMrxnFiltered), colnames(df_tmpLiverRNAnormFiltered))) break;
  
  if(length(unique(df_liverMetaData[df_liverMetaData$Age == int_age,"Batch"])) < 2) {
    fct_batch <- NULL
  } else {
    fct_batch <- df_liverMetaData[df_liverMetaData$Age == int_age,"Batch"]
  }
  
  # run correlation all 2 all without partial cor.
  df_tmpCorRxn2LiverRNA <- df_pcorAll2All(df_table1 = df_tmpLiverMMrxnFiltered,
                                          df_table2 = df_tmpLiverRNAnormFiltered,
                                          var_pcorZ = fct_batch,
                                          str_labels = c("Rxn","Gene"),
                                          int_cores = 15) #method = "spearman"
  saveRDS(df_tmpCorRxn2LiverRNA, paste0("df_", int_age, "monthCorRxn2LiverRNA.rds"))
  
  # how many do we actually have?
  print( table(df_tmpCorRxn2LiverRNA$p.adj <= 0.05) )
  #
  print( table( df_tmpCorRxn2LiverRNA$p.adj <= 0.1 &
                abs(df_tmpCorRxn2LiverRNA$Spearman_rho) >= 0.55 ))
  
  rm(df_tmpCorRxn2LiverRNA, df_tmpLiverMMrxnFiltered, df_tmpLiverRNAnormFiltered)
  gc()
}
###
# [1] 15
# [1] "Starting correlation of 26808 Genes to 2129 Rxns within 10 samples! Tue Feb 27 17:38:18 2024"
# [1] "Finished all correlations! Preparing results. Tue Feb 27 17:56:09 2024"
# FALSE     TRUE 
# 57065044     9188 
# 
# FALSE     TRUE 
# 57064932     9300 
###
# [1] 2
# [1] "Starting correlation of 26670 Genes to 2129 Rxns within 10 samples! Sat Jan 14 02:32:45 2023"
# [1] "Finished all correlations! Preparing results. Sat Jan 14 02:50:41 2023"
# FALSE     TRUE 
# 56769292    11138 
# 
# FALSE     TRUE 
# 56769182    11248 
###
# [1] 24
# [1] "Starting correlation of 28767 Genes to 2129 Rxns within 10 samples! Sat Jan 14 02:53:25 2023"
# [1] "Finished all correlations! Preparing results. Sat Jan 14 03:12:45 2023"
# 
#    FALSE     TRUE 
# 61235858     9085 
# 
#    FALSE     TRUE 
# 61235667     9276
# [1] 9
# [1] "Partial correlation enabled!"
# [1] "Starting partial correlation of 27038 Genes to 2129 Rxns within 10 samples! Tue Feb 27 19:10:45 2024"
# [1] "Finished all correlations! Preparing results. Wed Feb 28 00:13:16 2024"
# 
# FALSE     TRUE 
# 56984084   579818 
# 
# FALSE     TRUE 
# 56685971   877931 
# [1] 30
# [1] "Partial correlation enabled!"
# [1] "Starting partial correlation of 27965 Genes to 2129 Rxns within 12 samples! Wed Feb 28 00:20:51 2024"
# [1] "Finished all correlations! Preparing results. Wed Feb 28 05:35:57 2024"
# 
# FALSE     TRUE 
# 59152527   384958 
# 
# FALSE     TRUE 
# 58934622   602863 

dbl_percentCorRxn2Meta <- list(# cor
                               "2" = 11248/(56769182+11248)*100 , 
                               "15" = 9300/(57064932+9300)*100, 
                               "24" = 9276/(61235667+9276)*100,
                               # pcor
                               "9" = 877931/(56685971+877931)*100,
                               "30" = 602863/(58934622+602863)*100
                               )
pdf(file = "liver_hostMBcorrelationAgeStratified.pdf", width = 185/72, height = 185/72, pointsize = 10)
par(mar = c(3,5,3,1)+.1)
layout(matrix(c(1,1,1,2,2), nrow = 1))
barplot(height = as.numeric(dbl_percentCorRxn2Meta)[1:3],
        # names.arg = names(dbl_percentCorRxn2Meta),
        las = 2, ylim = c(0,0.025),
        main = "cor Host-MB", #"Correlations of Microbiome Reactions\n to Host Transcripts",
        col = brewer.pal("Pastel1", n = 5)[c(1,3,4)])
title(ylab = "Strong correlations [%]", line=4)
axis(side = 1, at = seq(from = 0.65, by = 1.225, length.out = 3),
     labels = names(dbl_percentCorRxn2Meta)[1:3], tick = F, line = F)
lines(x = c(0,1)+.75, y = rep(0.022,2)); text(x= 0.5+0.75, y= 0.023, labels = "***", cex = .8)
lines(x = c(1,2)+1, y = rep(0.02,2)); text(x= 1.5+1, y= 0.021, labels = "***", cex = .8)
# pcors
par(mar = c(3,2,3,3)+.1)
barplot(height = as.numeric(dbl_percentCorRxn2Meta)[4:5],
        # names.arg = names(dbl_percentCorRxn2Meta),
        las = 2, ylim = c(0,2.5),
        main = "pcor seq. batch", #"Partial Correlations of Microbiome Reactions\n to Host Transcripts corrected for Seq. Batch",
        col = brewer.pal("Pastel1", n = 5)[c(2,5)])
# title(ylab = "Strong correlations [%]", line=4)
axis(side = 1, at = seq(from = 0.65, by = 1.225, length.out = 2),
     labels = names(dbl_percentCorRxn2Meta)[4:5], tick = F, line = F)
lines(x = c(0,1)+.75, y = rep(1.7,2)); text(x= 0.5+0.75, y= 1.8, labels = "***", cex = .8)
dev.off()
# export at 440x400 barplot_strongCorHost2MetaLiverByAge

###########
# Pearson's chi-squared test:
# Null hypothesis: Joint distribution of the cell counts in the contingency table is the product of row and column sums (marginals)
####
# 2 to 15:
tbl_microbiomeCorrelation <- as.table( rbind(c(11248, 56769182), 
                                             c(9300, 57064932) ))
dimnames(tbl_microbiomeCorrelation) <- list(Age = c("2months", "15months"), 
                                            MicrobiomeCorrelation = c("strongMBcor", "weakMBcor"))
tbl_microbiomeCorrelation
#             MicrobiomeCorrelation
# Age        strongMBcor weakMBcor
# 2months          11248  56769182
# 15months          9300  57064932
####
# 15 to 24:
tbl_microbiomeCorrelation <- as.table( rbind(c(9300, 57064932),
                                             c(9276, 61235667) ))
dimnames(tbl_microbiomeCorrelation) <- list(Age = c("15months", "24months"), 
                                            MicrobiomeCorrelation = c("strongMBcor", "weakMBcor"))
####
# 9 to 30:
tbl_microbiomeCorrelation <- as.table( rbind(c(877931, 56685971),
                                             c(602863, 58934622) ))
dimnames(tbl_microbiomeCorrelation) <- list(Age = c("9months", "30months"), 
                                            MicrobiomeCorrelation = c("strongMBcor", "weakMBcor"))

####
# Pearson's Chi-squared test with Yates' continuity correction
obj_chiSqMBcor <- chisq.test(x = tbl_microbiomeCorrelation)
obj_chiSqMBcor
# data:  2 months to 15 months:
# X-squared = 194.71, df = 1, p-value < 2.2e-16
# data:  15 to 24 months
# X-squared = 24.766, df = 1, p-value = 6.473e-07
# data:  9 to 30 months
# X-squared = 61584, df = 1, p-value < 2.2e-16
obj_chiSqMBcor$expected     # expected counts under the null
obj_chiSqMBcor$residuals    # Pearson residuals  (positive numbers are overrepresented, negative numbers are underrepresented compared to null hpyothesis)

p2star(x = p.adjust(p = c(2.2e-16, 6.473e-07, 2.2e-16), method = "bonferroni") ) # "BH" 
# [1] "***" "***" "***"


