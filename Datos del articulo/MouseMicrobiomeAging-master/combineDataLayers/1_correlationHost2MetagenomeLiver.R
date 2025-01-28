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
df_liverMetaData <- as.data.frame( readRDS(file = "../mouseTranscriptome/df_transcriptMetaData.rds"))
df_liverMetaData <- df_liverMetaData[df_liverMetaData$Tissue == "Liver",]
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
# Variance stabilized RNA-Seq reads informed with Age as integer
obj_DifAbundLiverAllAges <- readRDS(file = "../mouseTranscriptome/obj_DifAbundLiverAllAges.rds")
df_liverRNAnorm <- as.data.frame( getVarianceStabilizedData(obj_DifAbundLiverAllAges) )
df_liverAgeDependance <- as.data.frame( results(obj_DifAbundLiverAllAges, name = "age", alpha = 0.05) )
rm(obj_DifAbundLiverAllAges)
# same order as metadata-table
df_liverRNAnorm <- df_liverRNAnorm[,df_liverMetaData$SampleName]


#########################################################
# ## make the 2 data layers connectable
# # Match rna with rxn data tables
# df_liverMMrxn <- df_liverMMrxn[,df_liverMetaData[df_liverMetaData$SampleName %in% colnames(df_liverRNAnorm),"SampleName"]]
# df_liverMetaData <- df_liverMetaData[colnames(df_liverMMrxn),]
# # match column names between tables
# df_liverRNAnorm <- df_liverRNAnorm[,df_liverMetaData$SampleName]
# check all three tables for matching sample identifiers
all.equal(colnames(df_liverRNAnorm), colnames(df_liverMMrxn), rownames(df_liverMetaData))
# TRUE


#####################
### remove features with near 0-variance (= near equal values across all samples)
## MB side
caret::nearZeroVar(x = t(df_liverMMrxn), saveMetrics = F)
# nothing to drop
df_liverMMrxnFiltered <- df_liverMMrxn
# df_liverMMrxnFiltered <- df_liverMMrxn[-caret::nearZeroVar(x = t(df_liverMMrxn), saveMetrics = F),]

## Host side
# caret::nearZeroVar(x = t(df_liverRNAnorm), saveMetrics = F)
# many features to drop
df_liverRNAnormFiltered <- df_liverRNAnorm[-caret::nearZeroVar(x = t(df_liverRNAnorm), saveMetrics = F),]

# Move the next step to high perfomance computing cluster - 32 cores - 4 hours
saveRDS(df_liverMMrxnFiltered, "./df_liverMMrxnFiltered.rds")
saveRDS(df_liverRNAnormFiltered, "./df_liverRNAnormFiltered.rds")
saveRDS(df_liverMetaData, "./df_liverMetaData.rds")
# run script "pcorLiverToMB.R"

rm(df_liverRNAnorm, df_liverMMrxn, df_liverMetaData)

#################
## partial correlation where age is regressed out
# df_pcorRxn2LiverRNA <- df_pcorAll2All(df_table1 = df_liverMMrxnFiltered[, df_liverMetaData$SampleID],
#                                       df_table2 = df_liverRNAnormFiltered[, df_liverMetaData$SampleID], 
#                                       var_pcorZ = df_liverMetaData$Age, 
#                                       str_labels = c("Rxn","Gene"), 
#                                       int_cores = 32)
# [1] "Partial correlation enabled!"
# [1] "Starting partial correlation of 26061 Genes to 2129 Rxns within 52 samples! Tue May 16 14:49:48 2023"
# [1] "Finished all correlations! Preparing results. Tue May 16 16:13:58 2023"


# saveRDS(df_pcorRxn2LiverRNA, "df_pcorRxn2LiverRNA.rds")

###############
# load back from cluster
df_pcorRxn2LiverRNA <- readRDS("df_pcorRxn2LiverRNA.rds")

# order by p-values
df_pcorRxn2LiverRNA <- df_pcorRxn2LiverRNA[order(df_pcorRxn2LiverRNA$p.val, abs(df_pcorRxn2LiverRNA$Spearman_rho), decreasing = c(F, T)),]

# p-value histogram
par(mar = c(5,4,4,2)+.1)
hist(df_pcorRxn2LiverRNA$p.val, breaks = 40)

# how many do we actually have?
table(df_pcorRxn2LiverRNA$p.adj <= 0.05)
#    FALSE     TRUE 
# 55482243     1626  

# with strong correlation?
table(df_pcorRxn2LiverRNA$p.adj <= 0.1 &
        abs(df_pcorRxn2LiverRNA$Spearman_rho) >= 0.55)
#    FALSE     TRUE 
# 55480444     3425 
# 
# # how many reactions are correlated well enough?
# length(unique(df_pcorRxn2LiverRNA[df_pcorRxn2LiverRNA$p.adj <= 0.05,"Rxn"])) / length(unique(df_pcorRxn2LiverRNA[,"Rxn"]))
# # 31.1 %
# 
# # correlated genes
# length(unique(df_pcorRxn2LiverRNA[df_pcorRxn2LiverRNA$p.adj <= 0.05,"Gene"])) / length(unique(df_pcorRxn2LiverRNA[,"Gene"]))
# # 1.4 %


################
# Select only significant results
df_pcorRxn2LiverRNATop <- df_pcorRxn2LiverRNA[df_pcorRxn2LiverRNA$p.adj <= 0.1 &
                                              abs(df_pcorRxn2LiverRNA$Spearman_rho) >= 0.55,] # Results in matrix of 1359 x 1277

# order results table by how often a gene or reaction is popping up 
df_pcorRxn2LiverRNATop <- df_pcorRxn2LiverRNATop[order(table(df_pcorRxn2LiverRNATop$Gene)[df_pcorRxn2LiverRNATop$Gene],
                                                       table(df_pcorRxn2LiverRNATop$Rxn)[df_pcorRxn2LiverRNATop$Rxn],
                                                       decreasing = T),]

# how many reactions and host transcripts correlated
print(c(length(unique(df_pcorRxn2LiverRNATop$Rxn)),
        length(unique(df_pcorRxn2LiverRNATop$Gene))))

# now we can delete the full table to safe on RAM
rm(df_pcorRxn2LiverRNA); gc()

