###########################################################
###
###   DESeq of RNA by Age, stratified by tissue
###
###                   L.B. 01-2024
###
###########################################################

################
##
library(DESeq2)
library(clusterProfiler)
source("customFunctions.R")
setwd("mouseTranscriptome/")


########
# Load metadata and counn data
df_transcriptMetaData <- readRDS("df_transcriptMetaData.rds")
mtx_CountsColonAllAges <- readRDS("mtx_CountsColonAllAges.rds")
mtx_CountsLiverAllAges <- readRDS("mtx_CountsLiverAllAges.rds")
mtx_CountsBrainAllAges <- readRDS("mtx_CountsBrainAllAges.rds")


####################
## Colon:
# select sample-metadata for this organ only
df_transcriptMetaDataColon <- df_transcriptMetaData[df_transcriptMetaData$Tissue == "Colon",]
rownames(df_transcriptMetaDataColon) <- df_transcriptMetaDataColon$SampleName
# scaling of the age variable
df_transcriptMetaDataColon$AgeScaled <- scale(df_transcriptMetaDataColon$Age, center = T, scale = T)

# visual QC
void_MDSplot(df_sampleMatrix = mtx_CountsColonAllAges[,df_transcriptMetaDataColon$id], 
             df_metaData = df_transcriptMetaDataColon, bol_labels = T, var_colorBy = "Age", 
             var_shapeBy = NULL, bol_continousColor = F)

# DESeq-object from existing count matrix
obj_DESeqColonAllAges <- DESeqDataSetFromMatrix(countData = mtx_CountsColonAllAges[,df_transcriptMetaDataColon$SampleName], 
                                                colData = df_transcriptMetaDataColon,
                                                design = ~ AgeScaled) 
obj_DifAbundColonAllAges <- DESeq(obj_DESeqColonAllAges)
# get result as data.frame
df_DifAbundColonAllAges <- as.data.frame( results(obj_DifAbundColonAllAges, independentFiltering = T, alpha = 0.05, name = "AgeScaled") )
# add gene-symbol as column
df_DifAbundColonAllAges$GeneSymbol <- df_mm10Ensembl2Genesymbol[rownames(df_DifAbundColonAllAges),"GeneSymbol"]
# order by p-values
df_DifAbundColonAllAges <- df_DifAbundColonAllAges[order(df_DifAbundColonAllAges$pvalue, decreasing = F),]
# drop NA-entries (usually genes where testing could not be done due to all-zero values)
df_DifAbundColonAllAges <- df_DifAbundColonAllAges[!is.na(df_DifAbundColonAllAges$log2FoldChange),]

# export significant hits
df_tmp <- df_DifAbundColonAllAges[df_DifAbundColonAllAges$padj <= 0.05,]
df_tmp <- df_tmp[!is.na(df_tmp$baseMean),]
write.table(x = df_tmp,
            file = paste0("DESeqColonAllAges_", Sys.Date(), ".tsv"),
            quote = F, sep = "\t", row.names = T, col.names = T)


########
## Liver:
df_transcriptMetaDataLiver <- df_transcriptMetaData[df_transcriptMetaData$Tissue == "Liver",]
rownames(df_transcriptMetaDataLiver) <- df_transcriptMetaDataLiver$SampleName

df_transcriptMetaDataLiver$AgeScaled <- scale(df_transcriptMetaDataLiver$Age, center = T, scale = T)

# visual QC
void_MDSplot(df_sampleMatrix = mtx_CountsLiverAllAges[,df_transcriptMetaDataLiver$id], 
             df_metaData = df_transcriptMetaDataLiver, bol_labels = F, var_colorBy = "Age", 
             var_shapeBy = "Batch", bol_continousColor = F)

# DESeq-object from existing count matrix
obj_DESeqLiverAllAges <- DESeqDataSetFromMatrix(countData = mtx_CountsLiverAllAges[,df_transcriptMetaDataLiver$SampleName], 
                                                colData = df_transcriptMetaDataLiver,
                                                design = ~ AgeScaled + Batch)
obj_DifAbundLiverAllAges <- DESeq(obj_DESeqLiverAllAges)
#
df_DifAbundLiverAllAges <- as.data.frame( results(obj_DifAbundLiverAllAges, independentFiltering = T, alpha = 0.05, name = "AgeScaled") )
df_DifAbundLiverAllAges$GeneSymbol <- df_mm10Ensembl2Genesymbol[rownames(df_DifAbundLiverAllAges),"GeneSymbol"]
df_DifAbundLiverAllAges <- df_DifAbundLiverAllAges[order(df_DifAbundLiverAllAges$pvalue, decreasing = F),]
df_DifAbundLiverAllAges <- df_DifAbundLiverAllAges[!is.na(df_DifAbundLiverAllAges$log2FoldChange),]
#
df_tmp <- df_DifAbundLiverAllAges[df_DifAbundLiverAllAges$padj <= 0.05,]
df_tmp <- df_tmp[!is.na(df_tmp$baseMean),]
write.table(x = df_tmp,
            file = paste0("DESeqLiverAllAges_", Sys.Date(), ".tsv"),
            quote = F, sep = "\t", row.names = T, col.names = T)


########
## Brain:
df_transcriptMetaDataBrain <- df_transcriptMetaData[df_transcriptMetaData$Tissue == "Brain",]
rownames(df_transcriptMetaDataBrain) <- df_transcriptMetaDataBrain$SampleName

df_transcriptMetaDataBrain$AgeScaled <- scale(df_transcriptMetaDataBrain$Age, center = T, scale = T)

# visual QC
void_MDSplot(df_sampleMatrix = mtx_CountsBrainAllAges[,df_transcriptMetaDataBrain$id], 
             df_metaData = df_transcriptMetaDataBrain, bol_labels = T, var_colorBy = "Age", 
             var_shapeBy = "Batch", bol_continousColor = F)

# DESeq-object from existing count matrix
obj_DESeqBrainAllAges <- DESeqDataSetFromMatrix(countData = mtx_CountsBrainAllAges[,df_transcriptMetaDataBrain$SampleName], 
                                                colData = df_transcriptMetaDataBrain,
                                                design = ~ AgeScaled + Batch)
obj_DifAbundBrainAllAges <- DESeq(obj_DESeqBrainAllAges)
#
df_DifAbundBrainAllAges <- as.data.frame( results(obj_DifAbundBrainAllAges, independentFiltering = T, alpha = 0.05, name = "AgeScaled") )
df_DifAbundBrainAllAges$GeneSymbol <- df_mm10Ensembl2Genesymbol[rownames(df_DifAbundBrainAllAges),"GeneSymbol"]
df_DifAbundBrainAllAges <- df_DifAbundBrainAllAges[order(df_DifAbundBrainAllAges$pvalue, decreasing = F),]
df_DifAbundBrainAllAges <- df_DifAbundBrainAllAges[!is.na(df_DifAbundBrainAllAges$log2FoldChange),]
#
df_tmp <- df_DifAbundBrainAllAges[df_DifAbundBrainAllAges$padj <= 0.05,]
df_tmp <- df_tmp[!is.na(df_tmp$baseMean),]
write.table(x = df_tmp,
            file = paste0("DESeqBrainAllAges_", Sys.Date(), ".tsv"),
            quote = F, sep = "\t", row.names = T, col.names = T)


#########################
## shared among all three organs
# going up with age
write.table(x = na.omit( Reduce(intersect, list(rownames(df_DifAbundColonAllAges)[df_DifAbundColonAllAges$padj <= 0.05 &
                                                                                  df_DifAbundColonAllAges$log2FoldChange > 0],
                                                rownames(df_DifAbundLiverAllAges)[df_DifAbundLiverAllAges$padj <= 0.05 &
                                                                                  df_DifAbundLiverAllAges$log2FoldChange > 0],
                                                rownames(df_DifAbundBrainAllAges)[df_DifAbundBrainAllAges$padj <= 0.05 &
                                                                                  df_DifAbundBrainAllAges$log2FoldChange > 0] ))),
            file = paste0("transcriptsUpInAgeColonLiverBrain_",Sys.Date(),".txt"),
            quote = F, sep = "\t", row.names = F, col.names = F)
# going down with age
write.table(x = na.omit( Reduce(intersect, list(rownames(df_DifAbundColonAllAges)[df_DifAbundColonAllAges$padj <= 0.05 &
                                                                                  df_DifAbundColonAllAges$log2FoldChange < 0],
                                                rownames(df_DifAbundLiverAllAges)[df_DifAbundLiverAllAges$padj <= 0.05 &
                                                                                  df_DifAbundLiverAllAges$log2FoldChange < 0],
                                                rownames(df_DifAbundBrainAllAges)[df_DifAbundBrainAllAges$padj <= 0.05 &
                                                                                  df_DifAbundBrainAllAges$log2FoldChange < 0] ))),
            file = paste0("transcriptsDownInAgeColonLiverBrain_",Sys.Date(),".txt"),
            quote = F, sep = "\t", row.names = F, col.names = F)
# GO shared organ background
write.table(x = na.omit( unique(c(rownames(df_DifAbundColonAllAges),
                                  rownames(df_DifAbundLiverAllAges),
                                  rownames(df_DifAbundBrainAllAges) ))),
            file = paste0("transcriptsBackgroundColonLiverBrain_",Sys.Date(),".txt"),
            quote = F, sep = "\t", row.names = F, col.names = F)

###############
# GO for shared lists
## Positive Correlations Genes
obj_Gene2GOPositive <- enricher(gene = na.omit( Reduce(intersect, list(df_DifAbundColonAllAges[df_DifAbundColonAllAges$padj <= 0.05 &
                                                                                               df_DifAbundColonAllAges$log2FoldChange > 0,"GeneSymbol"],
                                                                       df_DifAbundLiverAllAges[df_DifAbundLiverAllAges$padj <= 0.05 &
                                                                                               df_DifAbundLiverAllAges$log2FoldChange > 0,"GeneSymbol"],
                                                                       df_DifAbundBrainAllAges[df_DifAbundBrainAllAges$padj <= 0.05 &
                                                                                               df_DifAbundBrainAllAges$log2FoldChange > 0,"GeneSymbol"] ))),
                                minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                                TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                                TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                                universe = na.omit( unique(c(df_DifAbundColonAllAges$GeneSymbol,
                                                             df_DifAbundLiverAllAges$GeneSymbol,
                                                             df_DifAbundBrainAllAges$GeneSymbol)))
                                )
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GOPositive <- obj_Gene2GOPositive@result[order(obj_Gene2GOPositive@result$qvalue,
                                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GOPositive@result$ID], 
                                                       decreasing = c(F,F)),]
df_Gene2GOPositive <- df_Gene2GOPositive[df_Gene2GOPositive$p.adjust <= 0.05 &
                                         df_Gene2GOPositive$Count >= 3,]
# Incompatible GO-Id versions:
# df_Gene2GOPositive["GO:1902187","ID"] <- "GO:0044790" # same term but ID changed at some point O_o

write.table(x = df_Gene2GOPositive[,c(1:6,8)], 
            file = paste0("allOrgansUpByAgeSharedGOBioProcess_", Sys.Date(),".tsv"), 
            quote = F, sep = "\t", row.names = F, col.names = T)

######
## Negative Correlations Genes
obj_Gene2GONegative <- enricher(gene = na.omit( Reduce(intersect, list(df_DifAbundColonAllAges[df_DifAbundColonAllAges$padj <= 0.05 &
                                                                                               df_DifAbundColonAllAges$log2FoldChange < 0,"GeneSymbol"],
                                                                       df_DifAbundLiverAllAges[df_DifAbundLiverAllAges$padj <= 0.05 &
                                                                                               df_DifAbundLiverAllAges$log2FoldChange < 0,"GeneSymbol"],
                                                                       df_DifAbundBrainAllAges[df_DifAbundBrainAllAges$padj <= 0.05 &
                                                                                               df_DifAbundBrainAllAges$log2FoldChange < 0,"GeneSymbol"] ))),
                                minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                                TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                                TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                                universe = na.omit( unique(c(df_DifAbundColonAllAges$GeneSymbol,
                                                             df_DifAbundLiverAllAges$GeneSymbol,
                                                             df_DifAbundBrainAllAges$GeneSymbol)))
)
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GONegative <- obj_Gene2GONegative@result[order(obj_Gene2GONegative@result$qvalue,
                                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GONegative@result$ID], 
                                                       decreasing = c(F,F)),]
df_Gene2GONegative <- df_Gene2GONegative[df_Gene2GONegative$pvalue <= 0.05 &   # NOTE: CUTOFF IS WITH RAW PVALUE NOT FDR HERE!!!
                                         df_Gene2GONegative$Count >= 3,]
write.table(x = df_Gene2GONegative[,c(1:6,8)], 
            file = paste0("allOrgansDownByAgeSharedGOBioProcess_", Sys.Date(),".tsv"), 
            quote = F, sep = "\t", row.names = F, col.names = T)


##################################
## Plot a shared heatmap 
library(RColorBrewer)

str_ensemblIDs <- unique( df_mm10Ensembl2Genesymbol[df_mm10Ensembl2Genesymbol$GeneSymbol %in% 
                                              unlist(strsplit(x = df_Gene2GOPositive[,"geneID"], split = "/", fixed = T)),"EnsemblID"])
str_ensemblIDs <- c(str_ensemblIDs, 
                    unique( df_mm10Ensembl2Genesymbol[df_mm10Ensembl2Genesymbol$GeneSymbol %in% 
                                              unlist(strsplit(x = df_Gene2GONegative[,"geneID"], split = "/", fixed = T)),"EnsemblID"]))
par(mfrow = c(1,3))

mtx_tmp <- getVarianceStabilizedData(obj_DifAbundColonAllAges)[str_ensemblIDs,]
mtx_tmp <- mtx_tmp[order(df_DifAbundColonAllAges[str_ensemblIDs,"log2FoldChange"], decreasing = T),
                   order(df_transcriptMetaDataColon[colnames(mtx_tmp),"Age"],
                         hclust(dist(t(mtx_tmp)))$order)]
gplots::heatmap.2(x = mtx_tmp,
                  Colv = NULL,
                  Rowv = NULL,
                  col = c(brewer.pal("Blues", n = 9)[9:3], "#ffffff", brewer.pal("Reds", n = 9)[3:9]),
                  dendrogram = "n", trace = "n", scale = "row", labCol = "", key = F,
                  labRow = "", #df_mm10Ensembl2Genesymbol[rownames(mtx_tmp),"GeneSymbol"],
                  ColSideColors = brewer.pal("Pastel1", n = 5)[df_transcriptMetaDataColon[colnames(mtx_tmp),"AgeGroup"]]
)
# export at 400x800 heatmap_sharedTranscriptsColon (pdf: 4x6)
# liver
mtx_tmp <- getVarianceStabilizedData(obj_DifAbundLiverAllAges)[str_ensemblIDs,]
mtx_tmp <- mtx_tmp[order(df_DifAbundColonAllAges[str_ensemblIDs,"log2FoldChange"], decreasing = T), # sorting by colon FC is intentional here
                   order(df_transcriptMetaDataLiver[colnames(mtx_tmp),"Age"],
                         hclust(dist(t(mtx_tmp)))$order)]
gplots::heatmap.2(x = mtx_tmp,
                  Colv = NULL,
                  Rowv = NULL,
                  col = c(brewer.pal("Blues", n = 9)[9:3], "#ffffff", brewer.pal("Reds", n = 9)[3:9]),
                  dendrogram = "n", trace = "n", scale = "row", labCol = "", key = F,
                  labRow = "", #df_mm10Ensembl2Genesymbol[rownames(mtx_tmp),"GeneSymbol"],
                  ColSideColors = brewer.pal("Pastel1", n = 5)[df_transcriptMetaDataLiver[colnames(mtx_tmp),"AgeGroup"]]
)
# export at 400x800 heatmap_sharedTranscriptsLiver
# brain
mtx_tmp <- getVarianceStabilizedData(obj_DifAbundBrainAllAges)[str_ensemblIDs,]
mtx_tmp <- mtx_tmp[order(df_DifAbundColonAllAges[str_ensemblIDs,"log2FoldChange"], decreasing = T), # sorting by colon FC is intentional here
                   order(df_transcriptMetaDataBrain[colnames(mtx_tmp),"Age"],
                         hclust(dist(t(mtx_tmp)))$order)]
gplots::heatmap.2(x = mtx_tmp,
                  Colv = NULL,
                  Rowv = NULL,
                  col = c(brewer.pal("Blues", n = 9)[9:3], "#ffffff", brewer.pal("Reds", n = 9)[3:9]),
                  dendrogram = "n", trace = "n", scale = "row", labCol = "", key = F,
                  labRow = "",#df_mm10Ensembl2Genesymbol[rownames(mtx_tmp),"GeneSymbol"], 
                  ColSideColors = brewer.pal("Pastel1", n = 5)[df_transcriptMetaDataBrain[colnames(mtx_tmp),"AgeGroup"]]
)
# export at 400x800 heatmap_sharedTranscriptsBrain

# legend
gplots::heatmap.2(x = mtx_tmp,
                  Colv = NULL,
                  Rowv = NULL,
                  col = c(brewer.pal("Blues", n = 9)[9:3], "#ffffff", brewer.pal("Reds", n = 9)[3:9]),
                  dendrogram = "n", trace = "n", scale = "row", labCol = "", key = T, density.info = "n"
)
# 840x510 heatmapLegend

