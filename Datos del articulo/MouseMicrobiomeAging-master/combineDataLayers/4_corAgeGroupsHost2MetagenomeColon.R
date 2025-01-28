################################################################
###
###   MB rxns to host RNA correlation - stratified by AgeGroup
###
################################################################

library(RColorBrewer)
library(DESeq2)
library(gplots)
library(clusterProfiler)

source("customFunctions.R")

setwd("combineDataLayers/")


# #######################################
# ## load matching metadata
# df_colonMetaData <- as.data.frame( readRDS(file = "../mouseMetagenome/df_metaData.rds"))
# rownames(df_colonMetaData) <- df_colonMetaData$SampleName
# 
# 
# #######################################
# ## Load metabolic modeling derived reaction abundances (microbiome side - metagenomics based)
# df_colonMMrxn <- readRDS("df_rxnAbundancesMM.rds")
# # relabel samples from SeqID to SampleName
# df_colonMMrxn <- df_colonMMrxn[,df_colonMetaData$SequencingID]
# colnames(df_colonMMrxn) <- df_colonMetaData$SampleName
# # data is already normalized
# colSums(df_colonMMrxn)
# 
# 
# ######################
# ## load host side transcriptomics data
# # Variance stabilized RNA-Seq reads informed with Age as integer
# obj_DifAbundColonAllAges <- readRDS(file = "../mouseTranscriptome/obj_DifAbundColonAllAges.rds")
# df_colonRNAnorm <- as.data.frame( getVarianceStabilizedData(obj_DifAbundColonAllAges) )
# df_colonAgeDependance <- as.data.frame( results(obj_DifAbundColonAllAges, name = "age", alpha = 0.05) )
# rm(obj_DifAbundColonAllAges)
# # same order as metadata-table
# df_colonRNAnorm <- df_colonRNAnorm[,df_colonMetaData$SampleName]


#########################################################
# ## make the 2 data layers connectable
# # Match rna with rxn data tables
# df_colonMMrxn <- df_colonMMrxn[,df_colonMetaData[df_colonMetaData$SampleName %in% colnames(df_colonRNAnorm),"SampleName"]]
# df_colonMetaData <- df_colonMetaData[colnames(df_colonMMrxn),]
# # match column names between tables
# df_colonRNAnorm <- df_colonRNAnorm[,df_colonMetaData$SampleName]
# check all three tables for matching sample identifiers
all.equal(colnames(df_colonRNAnorm), colnames(df_colonMMrxn), rownames(df_colonMetaData))
# TRUE


############################
### For each Age-Group:
for(int_age in unique(df_colonMetaData$Age)) {
 print(int_age) 
  # remove rxn features with 0-variance (= equal values across all samples)
  if(!isEmpty(caret::nearZeroVar(x = t(df_colonMMrxn), saveMetrics = F))) {
    df_tmpColonMMrxnFiltered <- df_colonMMrxn[-caret::nearZeroVar(x = t(df_colonMMrxn[,df_colonMetaData[df_colonMetaData$Age == int_age,"SampleName"]]), saveMetrics = F),
                                              df_colonMetaData[df_colonMetaData$Age == int_age,"SampleName"]]
  }
  # nothing to discard
  else df_tmpColonMMrxnFiltered <- df_colonMMrxn[,df_colonMetaData[df_colonMetaData$Age == int_age,"SampleName"]]
  # remove transcript features with 0-variance (= equal values across all samples)
  df_tmpColonRNAnormFiltered <- df_colonRNAnorm[-caret::nearZeroVar(x = t(df_colonRNAnorm[,df_colonMetaData[df_colonMetaData$Age == int_age,"SampleName"]]), saveMetrics = F),
                                                df_colonMetaData[df_colonMetaData$Age == int_age,"SampleName"]]
  # same colnames? (matching samples in both tables)
  if(!all.equal(colnames(df_tmpColonMMrxnFiltered), colnames(df_tmpColonRNAnormFiltered))) break;
  
  # run correlation all 2 all without partial cor.
  df_tmpCorRxn2ColonRNA <- df_pcorAll2All(df_table1 = df_tmpColonMMrxnFiltered,
                                          df_table2 = df_tmpColonRNAnormFiltered,
                                          var_pcorZ = NULL,
                                          str_labels = c("Rxn","Gene"),
                                          int_cores = 15) #method = "spearman"
  saveRDS(df_tmpCorRxn2ColonRNA, paste0("df_", int_age, "monthCorRxn2ColonRNA.rds"))
  
  # how many do we actually have?
  print( table(df_tmpCorRxn2ColonRNA$p.adj <= 0.05) )
  #
  print( table( df_tmpCorRxn2ColonRNA$p.adj <= 0.1 &
                abs(df_tmpCorRxn2ColonRNA$Spearman_rho) >= 0.55 ))
  
  rm(df_tmpCorRxn2ColonRNA, df_tmpColonMMrxnFiltered, df_tmpColonRNAnormFiltered)
  gc()
}
###
# [1] 2
# [1] "Starting correlation of 31577 Genes to 2129 Rxns within 10 samples! Wed May 17 13:25:55 2023"
# [1] "Finished all correlations! Preparing results. Wed May 17 13:46:16 2023"
# 
# FALSE     TRUE 
# 67207592    19841 
# 
# FALSE     TRUE 
# 67207314    20119 
# 
# [1] 24
# [1] "Starting correlation of 31562 Genes to 2129 Rxns within 10 samples! Wed May 17 13:50:54 2023"
# [1] "Finished all correlations! Preparing results. Wed May 17 14:12:52 2023"
# 
# FALSE     TRUE 
# 67184890    10608 
# 
# FALSE     TRUE 
# 67184774    10724 
#
# [1] 30
# [1] "Starting correlation of 32221 Genes to 2129 Rxns within 12 samples! Wed May 17 14:18:03 2023"
# [1] "Finished all correlations! Preparing results. Wed May 17 14:40:21 2023"
# 
# FALSE     TRUE 
# 68588453    10056 
# 
# FALSE     TRUE 
# 68588262    10247 
# 
# [1] 9
# [1] "Starting correlation of 31677 Genes to 2129 Rxns within 10 samples! Wed May 17 14:45:30 2023"
# [1] "Finished all correlations! Preparing results. Wed May 17 15:07:00 2023"
# 
# FALSE     TRUE 
# 67426401    13932 
# 
# FALSE     TRUE 
# 67426115    14218 
# 
# [1] 15
# [1] "Starting correlation of 31508 Genes to 2129 Rxns within 10 samples! Wed May 17 15:13:08 2023"
# [1] "Finished all correlations! Preparing results. Wed May 17 15:34:25 2023"
# 
# FALSE     TRUE 
# 67067675  12857 
# 
# FALSE     TRUE 
# 67067531  13001 


###########
# Pearson's chi-squared test:
# Null hypothesis: Joint distribution of the cell counts in the contingency table is the product of row and column sums (marginals)
####
# 2 to 9:
tbl_microbiomeCorrelation <- as.table( rbind(c(19841, 67207592), 
                                             c(13932, 67426401) ))
dimnames(tbl_microbiomeCorrelation) <- list(Age = c("2months", "9months"), 
                                            MicrobiomeCorrelation = c("strongMBcor", "weakMBcor"))
tbl_microbiomeCorrelation
#         MicrobiomeCorrelation
# Age       strongMBcor weakMBcor
# 2months         19841  67207592
# 9months         13932  67426401
####
# 9 to 15:
tbl_microbiomeCorrelation <- as.table( rbind(c(13932, 67426401),
                                             c(12857, 67067675) ))
dimnames(tbl_microbiomeCorrelation) <- list(Age = c("9months", "15months"), 
                                            MicrobiomeCorrelation = c("strongMBcor", "weakMBcor"))
####
# 15 to 24:
tbl_microbiomeCorrelation <- as.table( rbind(c(12857, 67067675),
                                             c(10608, 67184890) ))
dimnames(tbl_microbiomeCorrelation) <- list(Age = c("15months", "24months"), 
                                            MicrobiomeCorrelation = c("strongMBcor", "weakMBcor"))
####
# 24 to 30:
tbl_microbiomeCorrelation <- as.table( rbind(c(10608, 67184890),
                                             c(10056, 68588453) ))
dimnames(tbl_microbiomeCorrelation) <- list(Age = c("24months", "30months"), 
                                            MicrobiomeCorrelation = c("strongMBcor", "weakMBcor"))
####
# Pearson's Chi-squared test with Yates' continuity correction
obj_chiSqMBcor <- chisq.test(x = tbl_microbiomeCorrelation)
obj_chiSqMBcor
# data:  2 months to 9 months:
# X-squared = 1052.5, df = 1, p-value < 2.2e-16
# data:  9 to 15 months
# X-squared = 37.512, df = 1, p-value = 9.085e-10
# data:  15 to 24 months
# X-squared = 219.27, df = 1, p-value < 2.2e-16
# data:  24 to 30 months
# X-squared = 28.291, df = 1, p-value = 1.044e-07

obj_chiSqMBcor$expected     # expected counts under the null
obj_chiSqMBcor$residuals    # Pearson residuals  (positive numbers are overrepresented, negative numbers are underrepresented compared to null hpyothesis)

p2star( p.adjust(p = c(2.2e-16, 9.085e-10, 2.2e-16, 1.044e-07), method = "bonferroni") ) # "BH" 
# [1] "***" "***"   "***" "***"

####
# now plot results:
dbl_percentCorRxn2Meta <- list("2" = 19841/(67207592+19841)*100 , "9" = 13932/(67426401+13932)*100,
                               "15" = 12857/(67067675+12857)*100, "24" = 10608/(67184890+10608)*100,
                               "30" = 10056/(68588453+10056)*100 )

# width = 220/72, height = 350/72, pointsize = 12
pdf(file = "colon_hostMBcorrelationAgeStratified.pdf", width = 185/72, height = 185/72, pointsize = 10)
par(mar = c(2,4,0,0)+.1)
barplot(height = as.numeric(dbl_percentCorRxn2Meta), 
        # names.arg = names(dbl_percentCorRxn2Meta),
        las = 2, ylim = c(0,0.034),
        # main = "Correlations of Microbiome Reactions\n to Host Transcripts",
        col = brewer.pal("Pastel1", n = 5),
        cex.axis = .8,
        yaxt = "n")
axis(side = 2, at = c(0,1,2,3)/100, las = 2, cex.axis = .8)
lines(x = c(0,1)+.75, y = rep(0.032,2)); text(x= 0.5+0.75, y= 0.033, labels = "***", cex = .8)
lines(x = c(1,2)+1, y = rep(0.023,2)); text(x= 1.5+1, y= 0.024, labels = "***", cex = .8)
lines(x = c(2,3)+1.25, y = rep(0.022,2)); text(x= 2.5+1.25, y= 0.023, labels = "***", cex = .8)
lines(x = c(3,4)+1.5, y = rep(0.021,2)); text(x= 3.5+1.5, y= 0.022, labels = "***", cex = .8)
title(ylab = "Strong correlations [%]", line=3)
axis(side = 1, at = seq(from = 0.65, by = 1.225, length.out = 5), 
     labels = names(dbl_percentCorRxn2Meta), tick = F, line = -1, cex.axis = .8)
title(xlab = "Age [months]", line = 1)
dev.off()
# export at 440x400 barplot_strongCorHost2MetaColonByAge

