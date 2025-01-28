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
df_brainMetaData <- as.data.frame( readRDS(file = "../mouseTranscriptome/df_transcriptMetaData.rds"))
df_brainMetaData <- df_brainMetaData[df_brainMetaData$Tissue == "Brain",]
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
# Variance stabilized RNA-Seq reads informed with Age as integer
obj_DifAbundBrainAllAges <- readRDS(file = "../mouseTranscriptome/obj_DifAbundBrainAllAges.rds")
df_brainRNAnorm <- as.data.frame( getVarianceStabilizedData(obj_DifAbundBrainAllAges) )
df_brainAgeDependance <- as.data.frame( results(obj_DifAbundBrainAllAges, name = "age", alpha = 0.05) )
rm(obj_DifAbundBrainAllAges)
# same order as metadata-table
df_brainRNAnorm <- df_brainRNAnorm[,df_brainMetaData$SampleName]


#########################################################
# ## make the 2 data layers connectable
# # Match rna with rxn data tables
# df_brainMMrxn <- df_brainMMrxn[,df_brainMetaData[df_brainMetaData$SampleName %in% colnames(df_brainRNAnorm),"SampleName"]]
# df_brainMetaData <- df_brainMetaData[colnames(df_brainMMrxn),]
# # match column names between tables
# df_brainRNAnorm <- df_brainRNAnorm[,df_brainMetaData$SampleName]
# check all three tables for matching sample identifiers
all.equal(colnames(df_brainRNAnorm), colnames(df_brainMMrxn), rownames(df_brainMetaData))
# TRUE


#####################
### remove features with near 0-variance (= near equal values across all samples)
## MB side
caret::nearZeroVar(x = t(df_brainMMrxn), saveMetrics = F)
# nothing to drop
df_brainMMrxnFiltered <- df_brainMMrxn
# df_brainMMrxnFiltered <- df_brainMMrxn[-caret::nearZeroVar(x = t(df_brainMMrxn), saveMetrics = F),]

## Host side
# caret::nearZeroVar(x = t(df_brainRNAnorm), saveMetrics = F)
# 25222 features to drop
df_brainRNAnormFiltered <- df_brainRNAnorm[-caret::nearZeroVar(x = t(df_brainRNAnorm), saveMetrics = F),]

# Move the next step to high perfomance computing cluster - 32 cores - 4 hours
saveRDS(df_brainMMrxnFiltered, "./df_brainMMrxnFiltered.rds")
saveRDS(df_brainRNAnormFiltered, "./df_brainRNAnormFiltered.rds")
saveRDS(df_brainMetaData, "./df_brainMetaData.rds")
# run script "HostMBcorrelations/pcorBrainToMB.R"

rm(df_brainRNAnorm, df_brainMMrxn, df_brainMetaData)

#################
## partial correlation where age is regressed out
# df_pcorRxn2BrainRNA <- df_pcorAll2All(df_table1 = df_brainMMrxnFiltered[, df_brainMetaData$SampleID],
#                                       df_table2 = df_brainRNAnormFiltered[, df_brainMetaData$SampleID], 
#                                       var_pcorZ = df_brainMetaData$Age, 
#                                       str_labels = c("Rxn","Gene"), 
#                                       int_cores = 32)
# [1] "Partial correlation enabled!"
# [1] "Starting partial correlation of 32294 Genes to 2129 Rxns within 52 samples! Tue May 16 14:51:22 2023"
# [1] "Finished all correlations! Preparing results. Tue May 16 16:34:23 2023"
# saveRDS(df_pcorRxn2BrainRNA, "df_pcorRxn2BrainRNA.rds")

###############
# load back from cluster
df_pcorRxn2BrainRNA <- readRDS("df_pcorRxn2BrainRNA.rds")

# order by p-values
df_pcorRxn2BrainRNA <- df_pcorRxn2BrainRNA[order(df_pcorRxn2BrainRNA$p.val, abs(df_pcorRxn2BrainRNA$Spearman_rho), decreasing = c(F, T)),]

# p-value histogram
hist(df_pcorRxn2BrainRNA$p.val, breaks = 40)

# how many do we actually have?
table(df_pcorRxn2BrainRNA$p.adj <= 0.05)
# FALSE 
# 68753926

# with strong correlation?
table(df_pcorRxn2BrainRNA$p.adj <= 0.1 &
      abs(df_pcorRxn2BrainRNA$Spearman_rho) >= 0.55)
# FALSE     TRUE 
# 68751427     2499

# # how many reactions are correlated well enough?
# length(unique(df_pcorRxn2BrainRNA[df_pcorRxn2BrainRNA$p.adj <= 0.05,"Rxn"])) / length(unique(df_pcorRxn2BrainRNA[,"Rxn"]))
# # 15.8 %
# 
# # correlated genes
# length(unique(df_pcorRxn2BrainRNA[df_pcorRxn2BrainRNA$p.adj <= 0.05,"Gene"])) / length(unique(df_pcorRxn2BrainRNA[,"Gene"]))
# # .8 %


################
# Select only significant results
df_pcorRxn2BrainRNATop <- df_pcorRxn2BrainRNA[df_pcorRxn2BrainRNA$p.adj <= 0.1 &
                                              abs(df_pcorRxn2BrainRNA$Spearman_rho) >= 0.55,] # Results in matrix of 1236 x 926

# order results table by how often a gene or reaction is popping up 
df_pcorRxn2BrainRNATop <- df_pcorRxn2BrainRNATop[order(table(df_pcorRxn2BrainRNATop$Gene)[df_pcorRxn2BrainRNATop$Gene],
                                                       table(df_pcorRxn2BrainRNATop$Rxn)[df_pcorRxn2BrainRNATop$Rxn],
                                                       decreasing = T),]

# how many reactions and host transcripts correlated
print(c(length(unique(df_pcorRxn2BrainRNATop$Rxn)),
        length(unique(df_pcorRxn2BrainRNATop$Gene))))

rm(df_pcorRxn2BrainRNA); gc()

