############################################################
###
###   Gene Ontology enrichment for Colon genes
###
############################################################

library(RColorBrewer)
library(clusterProfiler)

source("customFunctions.R")
setwd("mouseTranscriptome/")


################################
## Load matching metadata
df_colonMetaData


#################################
## Gene expression abundance data
df_colonRNAnormFiltered


#################################
## DESeq association of gene expression to age
df_DifAbundColonAllAges

### GO enrichment
## Positive Correlations Genes
obj_Gene2GOPositive <- enricher(gene = unique(na.omit(df_DifAbundColonAllAges[df_DifAbundColonAllAges$padj <= 0.05 &
                                                                              df_DifAbundColonAllAges$log2FoldChange > 0,"GeneSymbol"])),
                                minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                                TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                                TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                                universe = unique(na.omit(df_mm10Ensembl2Genesymbol[rownames(df_colonRNAnormFiltered),"GeneSymbol"])))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GOPositive <- obj_Gene2GOPositive@result[order(obj_Gene2GOPositive@result$qvalue,
                                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GOPositive@result$ID], 
                                                       decreasing = c(F,F)),]
df_Gene2GOPositive <- df_Gene2GOPositive[df_Gene2GOPositive$p.adjust <= 0.05 &
                                           df_Gene2GOPositive$Count >= 3,]
# Do we need to filter out a bunch of GO-terms?
df_Gene2GOPositiveColon <- df_trimGO(df_GOresultTable = df_Gene2GOPositive[df_Gene2GOPositive$p.adjust <= 1e-4,])
# saveRDS(df_Gene2GOPositiveColon,"df_Gene2GOPositiveColon.rds")
# df_Gene2GOPositiveColon <- readRDS("df_Gene2GOPositiveColon.rds")
table(df_Gene2GOPositiveColon$softTrim)
# FALSE  TRUE 
#    28    11
# df_Gene2GOPositiveColon <- df_Gene2GOPositiveColon[df_Gene2GOPositiveColon$softTrim != T,]
df_Gene2GOPositiveColon$ShortDescription <- str_simplifyDescriptions( df_Gene2GOPositiveColon$Description, bol_stripRomanNums = F )
unique(df_Gene2GOPositiveColon[df_Gene2GOPositiveColon$softTrim != T,"ShortDescription"] )
# 26

##########
## Negative Correlations Genes
obj_Gene2GONegative <- enricher(gene = unique(na.omit(df_DifAbundColonAllAges[df_DifAbundColonAllAges$padj <= 0.05 &
                                                                                df_DifAbundColonAllAges$log2FoldChange < 0,"GeneSymbol"])),
                                minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                                TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                                TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                                universe = unique(na.omit(df_mm10Ensembl2Genesymbol[rownames(df_colonRNAnormFiltered),"GeneSymbol"])))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GONegative <- obj_Gene2GONegative@result[order(obj_Gene2GONegative@result$qvalue,
                                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GONegative@result$ID], 
                                                       decreasing = c(F,F)),]
df_Gene2GONegative <- df_Gene2GONegative[df_Gene2GONegative$p.adjust <= 0.05 &
                                           df_Gene2GONegative$Count >= 3,]
# Do we need to filter out a bunch of GO-terms?
df_Gene2GONegativeColon <- df_trimGO(df_GOresultTable = df_Gene2GONegative[df_Gene2GONegative$p.adjust <= 1e-4,])
# saveRDS(df_Gene2GONegativeColon,"df_Gene2GONegativeColon.rds")
# df_Gene2GONegativeColon <- readRDS("df_Gene2GONegativeColon.rds")
table(df_Gene2GONegativeColon$softTrim)
# FALSE  
# 2
df_Gene2GONegativeColon$ShortDescription <- str_simplifyDescriptions( df_Gene2GONegativeColon$Description, bol_stripRomanNums = F )
unique(df_Gene2GONegativeColon[df_Gene2GONegativeColon$softTrim != T,"ShortDescription"] )
# 2

#########
## Make a heatmap of mean rxn activity from each GO-group 
df_ColonGO <- rbind(df_Gene2GOPositiveColon[df_Gene2GOPositiveColon$softTrim != T,],
                    df_Gene2GONegativeColon[df_Gene2GONegativeColon$softTrim != T,])
# # filter more stringent
# df_ColonGO <- df_ColonGO[df_ColonGO$p.adjust <= 1e-4,]

# combine duplicated descriptions on feature level
for(str_GO in df_ColonGO$ShortDescription) {
  df_ColonGO[df_ColonGO$ShortDescription == str_GO,"geneID"] <- paste(collapse = "/", 
                            unique(unlist(strsplit(x = df_ColonGO[df_ColonGO$ShortDescription == str_GO,"geneID"], "/", fixed = T))))
}
# drop duplicated GOs
df_ColonGO <- df_ColonGO[!duplicated(df_ColonGO$ShortDescription),]

# mean of feature abundance over each GO term
mtx_tmp2 <- matrix(data = NA, nrow = nrow(df_ColonGO), ncol = 5)
# cluster by shared features
rownames(mtx_tmp2) <- df_ColonGO[hclust(as.dist(mtx_mydist(df_ColonGO)))$order,"ShortDescription"]
# calc. mean of features
for(str_GOId in df_ColonGO$ID) {
  lst_genes <- unlist(strsplit(x = df_ColonGO[str_GOId,"geneID"], "/", fixed = T))
  lst_genes <- rownames(df_DifAbundColonAllAges)[df_DifAbundColonAllAges$GeneSymbol %in% lst_genes]
  
  # abundance of genes involved in this GO-term
  mtx_tmp <- df_colonRNAnormFiltered[na.omit(rownames(df_DifAbundColonAllAges)[df_DifAbundColonAllAges$padj <= 0.05]),
                                     df_colonMetaData$SampleName]
  mtx_tmp <- as.matrix(mtx_tmp[lst_genes, 
                               order(df_colonMetaData[colnames(mtx_tmp),"AgeGroup"])])
  # normalise matrix so each reaction is similarily abundant (percentage of abundance)
  mtx_tmp <- mtx_tmp / rowSums(mtx_tmp) * 1000
  
  mtx_tmp2[df_ColonGO[str_GOId,"ShortDescription"],] <- c( mean( mtx_tmp[, df_colonMetaData[colnames(mtx_tmp),"AgeGroup"] == 1 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_colonMetaData[colnames(mtx_tmp),"AgeGroup"] == 2 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_colonMetaData[colnames(mtx_tmp),"AgeGroup"] == 3 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_colonMetaData[colnames(mtx_tmp),"AgeGroup"] == 4 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_colonMetaData[colnames(mtx_tmp),"AgeGroup"] == 5 ], na.rm = T )
  )
}
View(mtx_tmp2)

# # export enrichment to supplements
# df_tmp <- merge(x = df_ColonGO, y = mtx_tmp2,
#                 by.x = "Description", by.y = "row.names")
# df_tmp <- df_tmp[,c(2,1,3:5,8,10:14)]
# colnames(df_tmp) <- c("TermID","SubSystem","FeatureRatio","BackgroundRatio","pvalue","FeatureIDs","Mean2months","Mean9months","Mean15months","Mean24months","Mean30months")
# write.table(x = df_tmp, file = "enrichmentAgeAssocInternRxnFlux.csv", 
#             quote = T, sep = ",", row.names = F, col.names = T)
# rm(df_tmp)

str_rowLabels <- rownames(mtx_tmp2)
str_rowLabels <- Hmisc::capitalize(gsub(pattern = "\n", replacement = " ", str_rowLabels))
str_rowLabels[20] <- "MHC class Ib antigen presentation"
str_rowLabels[22] <- "MHC class I antigen presentation via ER"
str_rowLabels[27] <- "MHC class II antigen presentation"

# rownames(mtx_tmp2) <- str_rowLabels
  # # drop duplicates by value
  # duplicated(mtx_tmp2)
  # mtx_tmp2 <- mtx_tmp2[!duplicated(mtx_tmp2),]
 
# # simplify rownames a bit more:
# rownames(mtx_tmp2)[1] <- "Folate synth."
# rownames(mtx_tmp2)[2] <- "Sphingosine metab."
# rownames(mtx_tmp2)[4] <- "Reductive acetyl-CoA pwy."
# rownames(mtx_tmp2)[9] <- "Adeninyl adenosylcobamide synth."
# rownames(mtx_tmp2)[13] <- "Phosphatidyl-serine & -ethanolamine synth."

# create a heatmap plot
pdf(file = "heatmap_colonGeneExpMeanGOAbundanceByAge.pdf", width = 360/72, height = 480/72, pointsize = 12)
gplots::heatmap.2(x = mtx_tmp2,
                  dendrogram = "none",
                  Colv = NA,
                  Rowv = NA,
                  scale = "row",
                  revC = T,
                  ColSideColors = brewer.pal("Pastel1", n = 5)[1:5],
                  col = colorRampPalette(brewer.pal("Purples", n=9))(50),
                  colsep = F, rowsep = F, 
                  key = T, trace = "n", density.info = "n", 
                  cexRow = 1.1,
                  labRow = str_rowLabels,#rownames(mtx_tmp2),
                  adjCol = c(1,.5),
                  labCol = c(2,9,15,24,30),
                  margins = c(5,22)
)
dev.off()

# rm(df_tmpMetadataFBA, obj_Rxn2SubsysPositive, df_Rxn2SubsysPositive, obj_Rxn2SubsysNegative, df_Rxn2SubsysNegative, str_rowLabels, mtx_tmp, mtx_tmp2, df_Rxn2Subsys)


##########################################
## Make a shared heatmap
# no soft trimming in this step
df_ColonGO <- rbind(df_Gene2GOPositiveColon,
                    df_Gene2GONegativeColon)
# combine duplicated descriptions on feature level
for(str_GO in df_ColonGO$ShortDescription) {
  df_ColonGO[df_ColonGO$ShortDescription == str_GO,"geneID"] <- paste(collapse = "/", 
                                                                      unique(unlist(strsplit(x = df_ColonGO[df_ColonGO$ShortDescription == str_GO,"geneID"], "/", fixed = T))))
}
# drop duplicated GOs
df_ColonGO <- df_ColonGO[!duplicated(df_ColonGO$ShortDescription),]
# mean of feature abundance over each GO term
mtx_tmp2 <- matrix(data = NA, nrow = nrow(df_ColonGO), ncol = 5)
# cluster by shared features
rownames(mtx_tmp2) <- df_ColonGO[hclust(as.dist(mtx_mydist(df_ColonGO)))$order,"ShortDescription"]
# calc. mean of features
for(str_GOId in df_ColonGO$ID) {
  lst_genes <- unlist(strsplit(x = df_ColonGO[str_GOId,"geneID"], "/", fixed = T))
  lst_genes <- rownames(df_DifAbundColonAllAges)[df_DifAbundColonAllAges$GeneSymbol %in% lst_genes]
  
  # abundance of genes involved in this GO-term
  mtx_tmp <- df_colonRNAnormFiltered[na.omit(rownames(df_DifAbundColonAllAges)[df_DifAbundColonAllAges$padj <= 0.05]),
                                     df_colonMetaData$SampleName]
  mtx_tmp <- as.matrix(mtx_tmp[lst_genes, 
                               order(df_colonMetaData[colnames(mtx_tmp),"AgeGroup"])])
  # normalise matrix so each reaction is similarily abundant (percentage of abundance)
  mtx_tmp <- mtx_tmp / rowSums(mtx_tmp) * 1000
  
  mtx_tmp2[df_ColonGO[str_GOId,"ShortDescription"],] <- c( mean( mtx_tmp[, df_colonMetaData[colnames(mtx_tmp),"AgeGroup"] == 1 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_colonMetaData[colnames(mtx_tmp),"AgeGroup"] == 2 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_colonMetaData[colnames(mtx_tmp),"AgeGroup"] == 3 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_colonMetaData[colnames(mtx_tmp),"AgeGroup"] == 4 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_colonMetaData[colnames(mtx_tmp),"AgeGroup"] == 5 ], na.rm = T )
  )
}
mtx_colonGOmeans <- mtx_tmp2; rm(mtx_tmp2)
#####
# Liver
df_LiverGO <- rbind(df_Gene2GOPositiveLiver,
                    df_Gene2GONegativeLiver)
# combine duplicated descriptions on feature level
for(str_GO in df_LiverGO$ShortDescription) {
  df_LiverGO[df_LiverGO$ShortDescription == str_GO,"geneID"] <- paste(collapse = "/", 
                                                                      unique(unlist(strsplit(x = df_LiverGO[df_LiverGO$ShortDescription == str_GO,"geneID"], "/", fixed = T))))
}
# drop duplicated GOs
df_LiverGO <- df_LiverGO[!duplicated(df_LiverGO$ShortDescription),]
# mean of feature abundance over each GO term
mtx_tmp2 <- matrix(data = NA, nrow = nrow(df_LiverGO), ncol = 5)
# cluster by shared features
rownames(mtx_tmp2) <- df_LiverGO[hclust(as.dist(mtx_mydist(df_LiverGO)))$order,"ShortDescription"]
# calc. mean of features
for(str_GOId in df_LiverGO$ID) {
  lst_genes <- unlist(strsplit(x = df_LiverGO[str_GOId,"geneID"], "/", fixed = T))
  lst_genes <- rownames(df_DifAbundLiverAllAges)[df_DifAbundLiverAllAges$GeneSymbol %in% lst_genes]
  
  # abundance of genes involved in this GO-term
  mtx_tmp <- df_liverRNAnormFiltered[na.omit(rownames(df_DifAbundLiverAllAges)[df_DifAbundLiverAllAges$padj <= 0.05]),
                                     df_liverMetaData$SampleName]
  mtx_tmp <- as.matrix(mtx_tmp[lst_genes, 
                               order(df_liverMetaData[colnames(mtx_tmp),"AgeGroup"])])
  # normalise matrix so each reaction is similarily abundant (percentage of abundance)
  mtx_tmp <- mtx_tmp / rowSums(mtx_tmp) * 1000
  
  mtx_tmp2[df_LiverGO[str_GOId,"ShortDescription"],] <- c( mean( mtx_tmp[, df_liverMetaData[colnames(mtx_tmp),"AgeGroup"] == 1 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_liverMetaData[colnames(mtx_tmp),"AgeGroup"] == 2 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_liverMetaData[colnames(mtx_tmp),"AgeGroup"] == 3 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_liverMetaData[colnames(mtx_tmp),"AgeGroup"] == 4 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_liverMetaData[colnames(mtx_tmp),"AgeGroup"] == 5 ], na.rm = T )
  )
}
mtx_liverGOmeans <- mtx_tmp2; rm(mtx_tmp2)
#####
# Brain
df_BrainGO <- rbind(df_Gene2GOPositiveBrain,
                    df_Gene2GONegativeBrain)
# combine duplicated descriptions on feature level
for(str_GO in df_BrainGO$ShortDescription) {
  df_BrainGO[df_BrainGO$ShortDescription == str_GO,"geneID"] <- paste(collapse = "/", 
                                                                      unique(unlist(strsplit(x = df_BrainGO[df_BrainGO$ShortDescription == str_GO,"geneID"], "/", fixed = T))))
}
# drop duplicated GOs
df_BrainGO <- df_BrainGO[!duplicated(df_BrainGO$ShortDescription),]
# mean of feature abundance over each GO term
mtx_tmp2 <- matrix(data = NA, nrow = nrow(df_BrainGO), ncol = 5)
# cluster by shared features
rownames(mtx_tmp2) <- df_BrainGO[hclust(as.dist(mtx_mydist(df_BrainGO)))$order,"ShortDescription"]
# calc. mean of features
for(str_GOId in df_BrainGO$ID) {
  lst_genes <- unlist(strsplit(x = df_BrainGO[str_GOId,"geneID"], "/", fixed = T))
  lst_genes <- rownames(df_DifAbundBrainAllAges)[df_DifAbundBrainAllAges$GeneSymbol %in% lst_genes]
  
  # abundance of genes involved in this GO-term
  mtx_tmp <- df_brainRNAnormFiltered[na.omit(rownames(df_DifAbundBrainAllAges)[df_DifAbundBrainAllAges$padj <= 0.05]),
                                     df_brainMetaData$SampleName]
  mtx_tmp <- as.matrix(mtx_tmp[lst_genes, 
                               order(df_brainMetaData[colnames(mtx_tmp),"AgeGroup"])])
  # normalise matrix so each reaction is similarily abundant (percentage of abundance)
  mtx_tmp <- mtx_tmp / rowSums(mtx_tmp) * 1000
  
  mtx_tmp2[df_BrainGO[str_GOId,"ShortDescription"],] <- c( mean( mtx_tmp[, df_brainMetaData[colnames(mtx_tmp),"AgeGroup"] == 1 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_brainMetaData[colnames(mtx_tmp),"AgeGroup"] == 2 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_brainMetaData[colnames(mtx_tmp),"AgeGroup"] == 3 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_brainMetaData[colnames(mtx_tmp),"AgeGroup"] == 4 ], na.rm = T ),
                                                           mean( mtx_tmp[, df_brainMetaData[colnames(mtx_tmp),"AgeGroup"] == 5 ], na.rm = T )
  )
}
mtx_brainGOmeans <- mtx_tmp2; rm(mtx_tmp2)
#####
# GOs in at least two:
lst_sharedGOs <- names(which(table(c(rownames(mtx_colonGOmeans), rownames(mtx_liverGOmeans), rownames(mtx_brainGOmeans))) >= 2))
# some redundancy trimming on only those terms that are shared across all three tissues
df_tmp <- df_LiverGO[df_LiverGO$ShortDescription %in% lst_sharedGOs,c("ID","Count","geneID")]
df_tmp <- df_trimGO(df_GOresultTable = df_tmp)
# remove from list
lst_sharedGOs <- lst_sharedGOs[!lst_sharedGOs %in% df_LiverGO[df_tmp[df_tmp$softTrim == T,"ID"],"ShortDescription"]]

# add missing GO terms to the other tables
mtx_tmp2 <- rbind(mtx_colonGOmeans,
                  matrix(data = NA, ncol = 5, nrow = sum(!lst_sharedGOs %in% rownames(mtx_colonGOmeans))) )
rownames(mtx_tmp2) <- c(rownames(mtx_colonGOmeans),
                        lst_sharedGOs[!lst_sharedGOs %in% rownames(mtx_colonGOmeans)])
mtx_colonGOmeans <- mtx_tmp2; rm(mtx_tmp2)
# liver
mtx_tmp2 <- rbind(mtx_liverGOmeans,
                  matrix(data = NA, ncol = 5, nrow = sum(!lst_sharedGOs %in% rownames(mtx_liverGOmeans))) )
rownames(mtx_tmp2) <- c(rownames(mtx_liverGOmeans),
                        lst_sharedGOs[!lst_sharedGOs %in% rownames(mtx_liverGOmeans)])
mtx_liverGOmeans <- mtx_tmp2; rm(mtx_tmp2)
# brain
mtx_tmp2 <- rbind(mtx_brainGOmeans,
                  matrix(data = NA, ncol = 5, nrow = sum(!lst_sharedGOs %in% rownames(mtx_brainGOmeans))) )
rownames(mtx_tmp2) <- c(rownames(mtx_brainGOmeans),
                        lst_sharedGOs[!lst_sharedGOs %in% rownames(mtx_brainGOmeans)])
mtx_brainGOmeans <- mtx_tmp2; rm(mtx_tmp2)

# combined matrix
mtx_sharedGOs <- cbind(t( apply(X = mtx_colonGOmeans[lst_sharedGOs,], MARGIN = 1, FUN = scale) ),
                       t( apply(X = mtx_liverGOmeans[lst_sharedGOs,], MARGIN = 1, FUN = scale) ),
                       t( apply(X = mtx_brainGOmeans[lst_sharedGOs,], MARGIN = 1, FUN = scale) ) )
# rearange by values
mtx_sharedGOs <- mtx_sharedGOs[hclust(dist(mtx_sharedGOs))$order,]
# colon NAs to bottom, then liver NAs, then brain NAs
mtx_sharedGOs <- mtx_sharedGOs[order(is.na(mtx_sharedGOs[,1]), is.na(mtx_sharedGOs[,6]), is.na(mtx_sharedGOs[,11])),]

# # export enrichment to supplements
# df_tmp <- merge(x = df_ColonGO, y = mtx_tmp2,
#                 by.x = "Description", by.y = "row.names")
# df_tmp <- df_tmp[,c(2,1,3:5,8,10:14)]
# colnames(df_tmp) <- c("TermID","SubSystem","FeatureRatio","BackgroundRatio","pvalue","FeatureIDs","Mean2months","Mean9months","Mean15months","Mean24months","Mean30months")
# write.table(x = df_tmp, file = "enrichmentAgeAssocInternRxnFlux.csv", 
#             quote = T, sep = ",", row.names = F, col.names = T)
# rm(df_tmp)

str_rowLabels <- rownames(mtx_sharedGOs)
str_rowLabels <- Hmisc::capitalize(gsub(pattern = "\n", replacement = " ", str_rowLabels))
str_rowLabels[str_rowLabels == "Antigen proc. and presentation of exogenous peptide antigen via MHC class II"] <- "MHC class II antigen presentation"
str_rowLabels

# create a heatmap plot
pdf(file = "heatmap_sharedGeneExpMeanGOAbundanceByAge.pdf", width = 550/72, height = 620/72, pointsize = 12)
gplots::heatmap.2(x = mtx_sharedGOs,
                  dendrogram = "none",
                  Colv = NA,
                  Rowv = NA,
                  scale = "row",
                  revC = F, na.color = "#d1d1d1",
                  ColSideColors = rep(brewer.pal("Pastel1", n = 5)[1:5], times = 3),
                  # col = colorRampPalette(brewer.pal("Purples", n=9))(50),
                  col = c(colorRampPalette(brewer.pal("Blues", n = 9)[7:2])(20), 
                          colorRampPalette(brewer.pal("Reds", n = 9)[2:7])(20)),
                  colsep = F, rowsep = F, 
                  key = T, trace = "n", density.info = "n", 
                  cexRow = 1.1,
                  labRow = str_rowLabels,#rownames(mtx_tmp2),
                  adjCol = c(1,.5),
                  # labCol = c(2,9,15,24,30),
                  margins = c(5,22)
)
dev.off()


#####################
## unique heatmaps of remaining GO terms
table(df_ColonGO[df_ColonGO$ShortDescription %in% rownames(mtx_colonGOmeans)[(!rownames(mtx_colonGOmeans) %in% lst_sharedGOs)],"softTrim"])
# FALSE  TRUE 
#     7     1
table(df_LiverGO[df_LiverGO$ShortDescription %in% rownames(mtx_liverGOmeans)[(!rownames(mtx_liverGOmeans) %in% lst_sharedGOs)],"softTrim"])
# FALSE  TRUE 
#    27     5 
table(df_BrainGO[df_BrainGO$ShortDescription %in% rownames(mtx_brainGOmeans)[(!rownames(mtx_brainGOmeans) %in% lst_sharedGOs)],"softTrim"])
# FALSE  TRUE 
#    15     2

####
# Colon
df_tmp <- df_ColonGO[df_ColonGO$ShortDescription %in% rownames(mtx_colonGOmeans)[(!rownames(mtx_colonGOmeans) %in% lst_sharedGOs)],]
df_tmp <- df_tmp[df_tmp$softTrim != T,]
#
mtx_tmp2 <- mtx_colonGOmeans[df_tmp$ShortDescription,]
#
str_rowLabels <- rownames(mtx_tmp2)
str_rowLabels <- Hmisc::capitalize(gsub(pattern = "\n", replacement = " ", str_rowLabels))
str_rowLabels[4] <- "MHC class Ib antigen presentation"
str_rowLabels[3] <- "MHC class I antigen presentation via ER"
# create a heatmap plot
pdf(file = "heatmap_colonGeneExpMeanGOAbundanceByAge.pdf", width = 360/72, height = 300/72, pointsize = 12)
gplots::heatmap.2(x = mtx_tmp2,
                  dendrogram = "none",
                  Colv = NA,
                  Rowv = NA,
                  scale = "row",
                  revC = T,
                  ColSideColors = brewer.pal("Pastel1", n = 5)[1:5],
                  col = c(colorRampPalette(brewer.pal("Blues", n = 9)[7:2])(20), 
                          colorRampPalette(brewer.pal("Reds", n = 9)[2:7])(20)),
                  colsep = F, rowsep = F, 
                  key = T, trace = "n", density.info = "n", 
                  cexRow = 1.1,
                  labRow = str_rowLabels,#rownames(mtx_tmp2),
                  adjCol = c(1,.5),
                  labCol = c(2,9,15,24,30),
                  margins = c(12,22)
)
dev.off()

#####
# Liver
df_tmp <- df_LiverGO[df_LiverGO$ShortDescription %in% rownames(mtx_liverGOmeans)[(!rownames(mtx_liverGOmeans) %in% lst_sharedGOs)],]
df_tmp <- df_tmp[df_tmp$softTrim != T,]
#
mtx_tmp2 <- mtx_liverGOmeans[df_tmp$ShortDescription,]
#
str_rowLabels <- rownames(mtx_tmp2)
str_rowLabels <- Hmisc::capitalize(gsub(pattern = "\n", replacement = " ", str_rowLabels))
str_rowLabels[25] <- "Mitochondrial complex I assembly"
str_rowLabels[26] <- "Proton driven ATP synth."
str_rowLabels[16] <- "Transmembrane tyrosine kinase sign. pwy."
str_rowLabels[1] <- "Immunoglobulin prod. in immune resp."
# create a heatmap plot
pdf(file = "heatmap_liverGeneExpMeanGOAbundanceByAge.pdf", width = 360/72, height = 540/72, pointsize = 12)
gplots::heatmap.2(x = mtx_tmp2,
                  dendrogram = "none",
                  Colv = NA,
                  Rowv = NA,
                  scale = "row",
                  revC = T,
                  ColSideColors = brewer.pal("Pastel1", n = 5)[1:5],
                  col = c(colorRampPalette(brewer.pal("Blues", n = 9)[7:2])(20), 
                          colorRampPalette(brewer.pal("Reds", n = 9)[2:7])(20)),
                  colsep = F, rowsep = F, 
                  key = T, trace = "n", density.info = "n", 
                  cexRow = 1.1,
                  labRow = str_rowLabels,#rownames(mtx_tmp2),
                  adjCol = c(1,.5),
                  labCol = c(2,9,15,24,30),
                  margins = c(5,22)
)
dev.off()

#####
# Brain
df_tmp <- df_BrainGO[df_BrainGO$ShortDescription %in% rownames(mtx_brainGOmeans)[(!rownames(mtx_brainGOmeans) %in% lst_sharedGOs)],]
df_tmp <- df_tmp[df_tmp$softTrim != T,]
#
mtx_tmp2 <- mtx_brainGOmeans[df_tmp$ShortDescription,]
#
str_rowLabels <- rownames(mtx_tmp2)
str_rowLabels <- Hmisc::capitalize(gsub(pattern = "\n", replacement = " ", str_rowLabels))
str_rowLabels[12] <- "Homophilic cell adhesion"  
# create a heatmap plot
pdf(file = "heatmap_brainGeneExpMeanGOAbundanceByAge.pdf", width = 360/72, height = 320/72, pointsize = 12)
gplots::heatmap.2(x = mtx_tmp2,
                  dendrogram = "none",
                  Colv = NA,
                  Rowv = NA,
                  scale = "row",
                  revC = T,
                  ColSideColors = brewer.pal("Pastel1", n = 5)[1:5],
                  col = c(colorRampPalette(brewer.pal("Blues", n = 9)[7:2])(20), 
                          colorRampPalette(brewer.pal("Reds", n = 9)[2:7])(20)),
                  colsep = F, rowsep = F, 
                  key = T, trace = "n", density.info = "n", 
                  cexRow = 1.1,
                  labRow = str_rowLabels,#rownames(mtx_tmp2),
                  adjCol = c(1,.5),
                  labCol = c(2,9,15,24,30),
                  margins = c(5,22)
)
dev.off()
