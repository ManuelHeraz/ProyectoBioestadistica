############################################################
###
###   Combine data layers with correlations or multiOmics
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
df_colonMetaData <- as.data.frame( readRDS(file = "../mouseMetagenome/df_metaData.rds"))
rownames(df_colonMetaData) <- df_colonMetaData$SampleName


#######################################
## Load metabolic modeling derived reaction abundances (microbiome side - metagenomics based)
df_colonMMrxn <- readRDS("df_rxnAbundancesMM.rds")
# relabel samples from SeqID to SampleName
df_colonMMrxn <- df_colonMMrxn[,df_colonMetaData$SequencingID]
colnames(df_colonMMrxn) <- df_colonMetaData$SampleName
# data is already normalized
colSums(df_colonMMrxn)


######################
## load host side transcriptomics data
# Variance stabilized RNA-Seq reads informed with Age as integer
obj_DifAbundColonAllAges <- readRDS(file = "../mouseTranscriptome/obj_DifAbundColonAllAges.rds")
df_colonRNAnorm <- as.data.frame( getVarianceStabilizedData(obj_DifAbundColonAllAges) )
df_colonAgeDependance <- as.data.frame( results(obj_DifAbundColonAllAges, name = "age", alpha = 0.05) )
rm(obj_DifAbundColonAllAges)
# same order as metadata-table
df_colonRNAnorm <- df_colonRNAnorm[,df_colonMetaData$SampleName]


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


#####################
### remove features with near 0-variance (= near equal values across all samples)
## MB side
caret::nearZeroVar(x = t(df_colonMMrxn), saveMetrics = F)
# nothing to drop
df_colonMMrxnFiltered <- df_colonMMrxn
# df_colonMMrxnFiltered <- df_colonMMrxn[-caret::nearZeroVar(x = t(df_colonMMrxn), saveMetrics = F),]

## Host side
# caret::nearZeroVar(x = t(df_colonRNAnorm), saveMetrics = F)
# 25222 features to drop
df_colonRNAnormFiltered <- df_colonRNAnorm[-caret::nearZeroVar(x = t(df_colonRNAnorm), saveMetrics = F),]

# Move the next step to high perfomance computing cluster - 32 cores - 4 hours
saveRDS(df_MMrxnFiltered, "./df_colonMMrxnFiltered.rds")
saveRDS(df_colonRNAnormFiltered, "./df_colonRNAnormFiltered.rds")
saveRDS(df_colonMetaData, "./df_colonMetaData.rds")
# run script "pcorColonToMB.R"

rm(df_colonRNAnorm, df_colonMMrxn, df_colonMetaData)

#################
## partial correlation where age is regressed out
# df_pcorRxn2ColonRNA <- df_pcorAll2All(df_table1 = df_colonMMrxnFiltered[, df_colonMetaData$SampleID],
#                                       df_table2 = df_colonRNAnormFiltered[, df_colonMetaData$SampleID], 
#                                       var_pcorZ = df_colonMetaData$Age, 
#                                       str_labels = c("Rxn","Gene"), 
#                                       int_cores = 32)
# [1] "Partial correlation enabled!"
# [1] "Starting partial correlation of 30249 Genes to 2129 Rxns within 52 samples! Tue May 16 14:12:04 2023"
# [1] "Finished all correlations! Preparing results. Tue May 16 15:38:07 2023"
# saveRDS(df_pcorRxn2ColonRNA, "df_pcorRxn2ColonRNA.rds")

###############
# load back from cluster
df_pcorRxn2ColonRNA <- readRDS("df_pcorRxn2ColonRNA.rds")
# order by p-values
df_pcorRxn2ColonRNA <- df_pcorRxn2ColonRNA[order(df_pcorRxn2ColonRNA$p.val, abs(df_pcorRxn2ColonRNA$Spearman_rho), decreasing = c(F, T)),]

# p-value histogram
par(mar = c(5,4,4,2)+.1)
hist(df_pcorRxn2ColonRNA$p.val, breaks = 40)

# how many do we actually have?
table(df_pcorRxn2ColonRNA$p.adj <= 0.05)
# FALSE 
# 64400121 

# with strong correlation?
table(df_pcorRxn2ColonRNA$p.adj <= 0.1 &
      abs(df_pcorRxn2ColonRNA$Spearman_rho) >= 0.55)
# FALSE     TRUE 
# 64387389    12732 

# # how many reactions are correlated well enough?
# length(unique(df_pcorRxn2ColonRNA[abs(df_pcorRxn2ColonRNA$Spearman_rho) >= 0.4,"Rxn"])) / length(unique(df_pcorRxn2ColonRNA[,"Rxn"]))
# # 67.9 %
# 
# # correlated genes
# length(unique(df_pcorRxn2ColonRNA[abs(df_pcorRxn2ColonRNA$Spearman_rho) >= 0.4,"Gene"])) / length(unique(df_pcorRxn2ColonRNA[,"Gene"]))
# # 5.6 %

################
# Select only significant results
df_pcorRxn2ColonRNATop <- df_pcorRxn2ColonRNA[df_pcorRxn2ColonRNA$p.adj <= 0.1 &
                                              abs(df_pcorRxn2ColonRNA$Spearman_rho) >= 0.55,]  # Results in matrix of 1606 x 2815
# order results table by how often a gene or reaction is popping up 
df_pcorRxn2ColonRNATop <- df_pcorRxn2ColonRNATop[order(table(df_pcorRxn2ColonRNATop$Gene)[df_pcorRxn2ColonRNATop$Gene],
                                                       table(df_pcorRxn2ColonRNATop$Rxn)[df_pcorRxn2ColonRNATop$Rxn],
                                                       decreasing = T),]

# now we can delete the full table to safe on RAM
rm(df_pcorRxn2ColonRNA); gc()


