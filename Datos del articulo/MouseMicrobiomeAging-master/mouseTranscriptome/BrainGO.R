############################################################
###
###   Gene Ontology enrichment for Brain genes
###
############################################################

library(RColorBrewer)
library(clusterProfiler)

source("customFunctions.R")
setwd("mouseTranscriptome/")


################################
## Load matching metadata
df_brainMetaData


#################################
## Gene expression abundance data
df_brainRNAnormFiltered


#################################
## DESeq association of gene expression to age
df_DifAbundBrainAllAges

### GO enrichment
## Positive Correlations Genes
obj_Gene2GOPositive <- enricher(gene = unique(na.omit(df_DifAbundBrainAllAges[df_DifAbundBrainAllAges$padj <= 0.05 &
                                                                                df_DifAbundBrainAllAges$log2FoldChange > 0,"GeneSymbol"])),
                                minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                                TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                                TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                                universe = unique(na.omit(df_mm10Ensembl2Genesymbol[rownames(df_brainRNAnormFiltered),"GeneSymbol"])))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GOPositive <- obj_Gene2GOPositive@result[order(obj_Gene2GOPositive@result$qvalue,
                                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GOPositive@result$ID], 
                                                       decreasing = c(F,F)),]
df_Gene2GOPositive <- df_Gene2GOPositive[df_Gene2GOPositive$p.adjust <= 0.05 &
                                         df_Gene2GOPositive$Count >= 3,]
# Do we need to filter out a bunch of GO-terms?
df_Gene2GOPositiveBrain <- df_trimGO(df_GOresultTable = df_Gene2GOPositive[df_Gene2GOPositive$p.adjust <= 1e-6,])
# saveRDS(df_Gene2GOPositiveBrain,"df_Gene2GOPositiveBrain.rds")
# df_Gene2GOPositiveBrain <- readRDS("df_Gene2GOPositiveBrain.rds")
table(df_Gene2GOPositiveBrain$softTrim)
# FALSE  TRUE 
#    30    12  
# df_Gene2GOPositiveBrain <- df_Gene2GOPositiveBrain[df_Gene2GOPositiveBrain$softTrim != T,]
df_Gene2GOPositiveBrain$ShortDescription <- str_simplifyDescriptions( df_Gene2GOPositiveBrain$Description, bol_stripRomanNums = F )
unique(df_Gene2GOPositiveBrain[df_Gene2GOPositiveBrain$softTrim != T,"ShortDescription"] )
# 29

##########
## Negative Correlations Genes
obj_Gene2GONegative <- enricher(gene = unique(na.omit(df_DifAbundBrainAllAges[df_DifAbundBrainAllAges$padj <= 0.05 &
                                                                                df_DifAbundBrainAllAges$log2FoldChange < 0,"GeneSymbol"])),
                                minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                                TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                                TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                                universe = unique(na.omit(df_mm10Ensembl2Genesymbol[rownames(df_brainRNAnormFiltered),"GeneSymbol"])))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GONegative <- obj_Gene2GONegative@result[order(obj_Gene2GONegative@result$qvalue,
                                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GONegative@result$ID], 
                                                       decreasing = c(F,F)),]
df_Gene2GONegative <- df_Gene2GONegative[df_Gene2GONegative$p.adjust <= 0.05 &
                                           df_Gene2GONegative$Count >= 3,]
# Do we need to filter out a bunch of GO-terms?
df_Gene2GONegativeBrain <- df_trimGO(df_GOresultTable = df_Gene2GONegative[df_Gene2GONegative$p.adjust <= 1e-6,])
# saveRDS(df_Gene2GONegativeBrain,"df_Gene2GONegativeBrain.rds")
# df_Gene2GONegativeBrain <- readRDS("df_Gene2GONegativeBrain.rds")
table(df_Gene2GONegativeBrain$softTrim)
# FALSE  TRUE 
#   10     3
df_Gene2GONegativeBrain$ShortDescription <- str_simplifyDescriptions( df_Gene2GONegativeBrain$Description, bol_stripRomanNums = F )
unique(df_Gene2GONegativeBrain[df_Gene2GONegativeBrain$softTrim != T,"ShortDescription"] )
# 10

#########
## Make a heatmap of mean rxn activity from each GO-group 
df_BrainGO <- rbind(df_Gene2GOPositiveBrain[df_Gene2GOPositiveBrain$softTrim != T,],
                    df_Gene2GONegativeBrain[df_Gene2GONegativeBrain$softTrim != T,])
# # filter more stringent
# df_BrainGO <- df_BrainGO[df_BrainGO$p.adjust <= 1e-6,]

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
View(mtx_tmp2)

# # export enrichment to supplements
# df_tmp <- merge(x = df_BrainGO, y = mtx_tmp2,
#                 by.x = "Description", by.y = "row.names")
# df_tmp <- df_tmp[,c(2,1,3:5,8,10:14)]
# colnames(df_tmp) <- c("TermID","SubSystem","FeatureRatio","BackgroundRatio","pvalue","FeatureIDs","Mean2months","Mean9months","Mean15months","Mean24months","Mean30months")
# write.table(x = df_tmp, file = "enrichmentAgeAssocInternRxnFlux.csv", 
#             quote = T, sep = ",", row.names = F, col.names = T)
# rm(df_tmp)

str_rowLabels <- rownames(mtx_tmp2)
str_rowLabels <- Hmisc::capitalize(gsub(pattern = "\n", replacement = " ", str_rowLabels))
# str_rowLabels[38] <- "MHC class I antigen presentation via ER"
# str_rowLabels[39] <- "MHC class Ib antigen presentation"
str_rowLabels[3] <- "Homophilic cell adhesion"  
str_rowLabels[38] <- "MHC class II antigen presentation"
# str_rowLabels[50] <- "MHC class II antigen assembly"
# str_rowLabels[51] <- "Immunoglobulin prod. in immune resp."

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

# # scale row-wise
# mtx_tmp2 <- t(apply(X = mtx_tmp2, MARGIN = 1, FUN = scales::rescale, to = c(-1,1)))

# create a heatmap plot
pdf(file = "heatmap_brainGeneExpMeanGOAbundanceByAge.pdf", width = 360/72, height = 660/72, pointsize = 12)
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
                  #breaks = 51,
                  margins = c(5,22)
)
dev.off()

# rm(df_tmpMetadataFBA, obj_Rxn2SubsysPositive, df_Rxn2SubsysPositive, obj_Rxn2SubsysNegative, df_Rxn2SubsysNegative, str_rowLabels, mtx_tmp, mtx_tmp2, df_Rxn2Subsys)



