##############################################
##
## MB-correlated genes shared among organs?
##
##############################################

library(RColorBrewer)
library(clusterProfiler)
library(UpSetR)

source("customFunctions.R")
setwd("combineDataLayers/")


##########################
## Display overlap between mb-host and age associated genes via upset plot
# Dataset
int_upset <- c(
  # Colon_Age = sum( df_colonAgeDependance$padj <= 0.05, na.rm = T ),
  # Colon_MB = length(unique(df_pcorRxn2ColonRNATop$Gene)),
  "Colon_Age&Colon_MB" =  length(intersect(unique(df_pcorRxn2ColonRNATop$Gene),
                                           na.omit(rownames(df_colonAgeDependance)[df_colonAgeDependance$padj <= 0.05]))),
  # Brain_Age = sum( df_brainAgeDependance$padj <= 0.05, na.rm = T ),
  # Brain_MB = length(unique(df_pcorRxn2BrainRNATop$Gene)),
  "Brain_Age&Brain_MB" =  length(intersect(unique(df_pcorRxn2BrainRNATop$Gene),
                                           na.omit(rownames(df_brainAgeDependance)[df_brainAgeDependance$padj <= 0.05]))),
  # Liver_Age = sum( df_liverAgeDependance$padj <= 0.05, na.rm = T ),
  # Liver_MB = length(unique(df_pcorRxn2LiverRNATop$Gene))
  "Liver_Age&Liver_MB" =  length(intersect(unique(df_pcorRxn2LiverRNATop$Gene),
                                           na.omit(rownames(df_liverAgeDependance)[df_liverAgeDependance$padj <= 0.05]))),
  
  "Colon_MB&Brain_MB" = length(intersect(unique(df_pcorRxn2ColonRNATop$Gene),
                                         unique(df_pcorRxn2BrainRNATop$Gene))),
  "Colon_MB&Liver_MB" = length(intersect(unique(df_pcorRxn2ColonRNATop$Gene),
                                         unique(df_pcorRxn2LiverRNATop$Gene))),
  "Liver_MB&Brain_MB" = length(intersect(unique(df_pcorRxn2LiverRNATop$Gene),
                                         unique(df_pcorRxn2BrainRNATop$Gene))),
  "Colon_MB&Brain_MB&Liver_MB" = length(Reduce(intersect, list(unique(df_pcorRxn2ColonRNATop$Gene),
                                                               unique(df_pcorRxn2LiverRNATop$Gene),
                                                               unique(df_pcorRxn2BrainRNATop$Gene)))),
  
  "Colon_Age&Brain_Age&Liver_Age" = length(Reduce(intersect, list(unique( na.omit( rownames(df_colonAgeDependance)[df_colonAgeDependance$padj <= 0.05] )),
                                                                  unique( na.omit( rownames(df_brainAgeDependance)[df_brainAgeDependance$padj <= 0.05] )),
                                                                  unique( na.omit( rownames(df_liverAgeDependance)[df_liverAgeDependance$padj <= 0.05] )) ))),
  "Colon_Age&Brain_Age" = length(intersect( unique( na.omit( rownames(df_colonAgeDependance)[df_colonAgeDependance$padj <= 0.05] )),
                                            unique( na.omit( rownames(df_brainAgeDependance)[df_brainAgeDependance$padj <= 0.05] )) )),
  "Liver_Age&Brain_Age" = length(intersect( unique( na.omit( rownames(df_liverAgeDependance)[df_liverAgeDependance$padj <= 0.05] )),
                                            unique( na.omit( rownames(df_brainAgeDependance)[df_brainAgeDependance$padj <= 0.05] )) )),
  "Colon_Age&Liver_Age" = length(intersect( unique( na.omit( rownames(df_colonAgeDependance)[df_colonAgeDependance$padj <= 0.05] )),
                                            unique( na.omit( rownames(df_liverAgeDependance)[df_liverAgeDependance$padj <= 0.05] )) ))
  
)
df_upsetMetadata <- data.frame(sets = c('Colon_Age', 'Colon_MB', 'Brain_Age', 'Brain_MB', 'Liver_Age', 'Liver_MB'),
                               setSize = c(sum( df_colonAgeDependance$padj <= 0.05, na.rm = T ), 
                                           length(unique(df_pcorRxn2ColonRNATop$Gene)),
                                           sum( df_brainAgeDependance$padj <= 0.05, na.rm = T ),
                                           length(unique(df_pcorRxn2BrainRNATop$Gene)),
                                           sum( df_liverAgeDependance$padj <= 0.05, na.rm = T ),
                                           length(unique(df_pcorRxn2LiverRNATop$Gene))), stringsAsFactors = F )
df_upsetMetadata$setSize <- as.numeric(df_upsetMetadata$setSize)
# Not really an upset plot - since each column is a unique overlap pair (you can't sum up the colSums)
pdf(file = "upsettingPlot_overlapHostGenes.pdf", width = 700/72, height = 500/72, pointsize = 12)
upset(fromExpression(int_upset),
      nsets = 6, 
      # sets = c('Colon_Age', 'Colon_MB', 'Brain_Age', 'Brain_MB', 'Liver_Age', 'Liver_MB'),
      sets = rev( c('Colon_Age', 'Brain_Age', 'Liver_Age', 'Colon_MB', 'Brain_MB', 'Liver_MB') ),
      keep.order = T,
      order.by = c("degree", "freq"),
      decreasing = c(T,T), 
      
      set.metadata = list(data = df_upsetMetadata,
                          plots = list(
                            list(type = "hist", column = "setSize", assign = 20)
                          )
      )
)
dev.off()


##########################
### Colon - Overlap of age and MB-assoc. genes
## Hypergeometric test for enrichment of overlap in age-associated to MB-associated colon transcripts
phyper(# Shared genes (overlap): 588-1
  q = length(intersect(unique(df_pcorRxn2ColonRNATop$Gene),
                       na.omit(rownames(df_colonAgeDependance)[df_colonAgeDependance$padj <= 0.05]))) - 1, 
  # MB correlated genes: 2815
  m = length(unique(df_pcorRxn2ColonRNATop$Gene)), 
  # total (background) - MB Correlated genes: 30249 - 2815 = 27434
  n = nrow(df_colonRNAnormFiltered)-length(unique(df_pcorRxn2ColonRNATop$Gene)), 
  # Age assoc. genes: 4715
  k = sum( df_colonAgeDependance$padj <= 0.05, na.rm = T ),
  lower.tail = F)
# probablity of seing an equal or bigger overlap:
# [1] 1.647321e-15
# There is a significant enrichment of MB-associated Genes in Age-associated Genes in Colon

# export results table
df_tmp <- merge(x = df_colonAgeDependance[df_colonAgeDependance$padj <= 0.05,c("log2FoldChange","GeneSymbol")],
                y = df_pcorRxn2ColonRNATop[,c("Gene","Rxn","Spearman_rho")],
                by.x = 0, by.y = "Gene")
colnames(df_tmp) <- c("Host_EnsemblID", "Age_log2FC", "Host_GeneSymbol", "MB_ReactionID", "MB_SpearmanRho")
df_tmp <- df_tmp[,c("Host_EnsemblID", "Host_GeneSymbol", "Age_log2FC", "MB_ReactionID", "MB_SpearmanRho")]
export(df_tmp, filename = "mtx_ColonSharedAgeMBAssocGenes.csv")

######
### Ontology Enrichment for significant correlations only!
## Colon age pos. log2FC (up in old) and MB Correlations Genes 
obj_Gene2GO <- enricher(gene = unique(na.omit(df_mm10Ensembl2Genesymbol[ intersect(unique(df_pcorRxn2ColonRNATop$Gene),
                                                                                   na.omit(rownames(df_colonAgeDependance)[df_colonAgeDependance$padj <= 0.05 &
                                                                                                                           df_colonAgeDependance$log2FoldChange > 0]) ),
                                                                         "GeneSymbol"])),
                        minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                        TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                        TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                        universe = unique(na.omit(df_mm10Ensembl2Genesymbol[rownames(df_colonRNAnormFiltered),"GeneSymbol"])))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GO <- obj_Gene2GO@result[order(obj_Gene2GO@result$qvalue,
                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GO@result$ID], 
                                       decreasing = c(F,F)),]
df_Gene2GOPositive <- df_Gene2GO[df_Gene2GO$p.adjust <= 0.05 &
                                 df_Gene2GO$Count >= 3,]
df_Gene2GOPositive <- df_trimGO(df_GOresultTable = df_Gene2GOPositive)

#####
## Colon age neg. log2FC (up in old) and MB Correlations Genes 
obj_Gene2GO <- enricher(gene = unique(na.omit(df_mm10Ensembl2Genesymbol[ intersect(unique(df_pcorRxn2ColonRNATop$Gene),
                                                                                   na.omit(rownames(df_colonAgeDependance)[df_colonAgeDependance$padj <= 0.05 &
                                                                                                                           df_colonAgeDependance$log2FoldChange < 0]) ),
                                                                         "GeneSymbol"])),
                        minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                        TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                        TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                        universe = unique(na.omit(df_mm10Ensembl2Genesymbol[rownames(df_colonRNAnormFiltered),"GeneSymbol"])))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GO <- obj_Gene2GO@result[order(obj_Gene2GO@result$qvalue,
                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GO@result$ID], 
                                       decreasing = c(F,F)),]
df_Gene2GONegative <- df_Gene2GO[df_Gene2GO$p.adjust <= 0.05 &
                                 df_Gene2GO$Count >= 3,]
df_Gene2GONegative <- df_trimGO(df_GOresultTable = df_Gene2GONegative)
# combine
df_Gene2GO <- rbind(df_Gene2GOPositive[!df_Gene2GOPositive$softTrim,],
                    df_Gene2GONegative[!df_Gene2GONegative$softTrim,])
df_Gene2GOcolon <- df_Gene2GO
# df_Gene2GO <- df_Gene2GOcolon
# add MB correlation
df_Gene2GO$MedianMBcor <- apply(X = df_Gene2GO, MARGIN = 1, FUN = function(obj_row) {
  lst_genes <-unique( unlist(strsplit(split = "/", fixed = T, x = obj_row[["geneID"]])) )
  lst_genes <- df_mm10Ensembl2Genesymbol[df_mm10Ensembl2Genesymbol$GeneSymbol %in% lst_genes, "EnsemblID"]
  # collect median MB-correlation for each gene
  df_tmp <- df_pcorRxn2ColonRNATop[df_pcorRxn2ColonRNATop$Gene %in% lst_genes,]
  df_tmp <- tapply(X = df_tmp$Spearman_rho, INDEX = df_tmp$Gene, FUN = median)
  # summarize median correlations for all genes
  return( median(df_tmp) )
})
## add Age foldChange
# now get median log2FC across all genes related to one GO-term
df_Gene2GO$MedianAgeLog2FC <- apply(X = df_Gene2GO, MARGIN = 1, FUN = function(obj_row) {
  lst_genes <-unique( unlist(strsplit(split = "/", fixed = T, x = obj_row[["geneID"]])) )
  lst_genes <- df_mm10Ensembl2Genesymbol[df_mm10Ensembl2Genesymbol$GeneSymbol %in% lst_genes, "EnsemblID"]
  # summarize age log2FC for all genes with median 
  return( median(df_DifAbundColonAllAges[lst_genes,"log2FoldChange"]) ) #Log2FC2to30
})

# order by log2FC
df_Gene2GO <- df_Gene2GO[order(df_Gene2GO$MedianAgeLog2FC),]
# shorten GO term labels
df_Gene2GO$DescriptionSplit <- str_simplifyDescriptions(df_Gene2GO$Description, bol_stripRomanNums = F)
df_Gene2GO$DescriptionSplit
# fix some labels that are still too long for plotting
df_Gene2GO$DescriptionSplit[3] <- "attachment of spindle to kinetochore"
df_Gene2GO$DescriptionSplit[8] <- "mitotic cell cycle spindle assembly"
df_Gene2GO$DescriptionSplit[9] <- "LDL particle remodeling"
df_Gene2GO$DescriptionSplit[14] <- "HDL particle remodeling"
df_Gene2GO$DescriptionSplit[18] <- "MHC class II antigen presentation"
df_Gene2GO$DescriptionSplit[20] <- "MHC class I antigen presentation via ER pwy."
df_Gene2GO$DescriptionSplit[21] <- "MHC class Ib antigen presentation"
df_Gene2GO$DescriptionSplit[46] <- "MHC class II antigen assembly"                             
df_Gene2GO$DescriptionSplit[47] <- "immunoglobulin prod. in immune resp."     

pdf(file = "GOplot_colonAgeAndMBassociatedGenesGO.pdf", width = 290/72, height = 540/72, pointsize = 8)
par(mar = c(4,22,1,1)+.1)
plot(x = df_Gene2GO$MedianAgeLog2FC,
     xlab = "Age log2(FC)", xlim = c(-.2,.5), xaxt = "n",
     y = 1:nrow(df_Gene2GO), ylab = "", yaxt = "n", lwd = .8,
     pch = 21, cex = scales::rescale(-log10(df_Gene2GO$pvalue), 
                                     to = c(1,3)),
     bg = brewer.pal("RdBu", n= 11)[11:1][round(scales::rescale(df_Gene2GO$MedianMBcor, from = c(-1,1), to = c(1,11)))],
     main = "Colon"
)
axis(side = 2, at = 1:nrow(df_Gene2GO), labels = Hmisc::capitalize( df_Gene2GO$DescriptionSplit ), las = 2, cex.axis = 1.1)
axis(side = 1, at = seq(from = -2,to = 5, by = 1)/10, labels = c(-2,NA,0,NA,2,NA,NA,5)/10)
abline(v = 0, col = rgb(.3,.3,.3,.7), lty = 4, lwd = 1)
dev.off()
# legend("bottomright", legend = c("pos. MB cor.", "neg. MB cor.", "size = enrich. p-val."),
#        pt.bg = c("#2166AC", "#B2182B"), pch = c(21,21,NA), bty = "n", cex = .7)

### horizontal version:
pdf(file = "GOplot_colonAgeAndMBassociatedGenesGO_horiz.pdf", width = 570/72, height = 140/72, pointsize = 6)
par(mar = c(12,10,1,1)+.1)
plot(x = 1:nrow(df_Gene2GO), 
     y = df_Gene2GO$MedianAgeLog2FC,
     ylim = c(-.2,.5), xlim = c(2,nrow(df_Gene2GO)-1),
     xlab = "", xaxt = "n", cex.lab = 1,
     ylab = "Age log2(FC)", yaxt = "n", lwd = .6,
     pch = 21, cex = scales::rescale(-log10(df_Gene2GO$pvalue), 
                                     to = c(1,3)),
     bg = brewer.pal("RdBu", n= 11)[11:1][round(scales::rescale(df_Gene2GO$MedianMBcor, from = c(-1,1), to = c(1,11)))],
     axes = F
     # main = "Colon"
)
axis(side = 1, at = 1:nrow(df_Gene2GO), labels = NA, lwd = .6, lwd.tick = .6)
text(x = 1:nrow(df_Gene2GO), y = -.32, xpd = T,
     labels = Hmisc::capitalize( df_Gene2GO$DescriptionSplit ),
     srt = 40, adj = c(1,0.5), cex = 1)
axis(side = 2, at = seq(from = -2,to = 5, by = 1)/10, labels = c(-2,NA,0,NA,2,NA,NA,5)/10,
     las = 2, cex.axis = 1, cex.lab = 1, lwd = .6, lwd.tick = .6)
abline(h = 0, col = rgb(.3,.3,.3,.7), lty = 4, lwd = .6)
dev.off()

# export GO table for supplements
export(df_Gene2GO[,c(1:6,8,9,12,13)], filename = "mtx_ColonSharedAgeMBGOPathways.csv")


###################
### Brain
## Hypergeometric test for enrichment of overlap in age-associated to MB-associated brain transcripts
phyper(# Shared genes (overlap): 260-1
  q = length(intersect(unique(df_pcorRxn2BrainRNATop$Gene),
                       na.omit(rownames(df_brainAgeDependance)[df_brainAgeDependance$padj <= 0.05]))) - 1, 
  # MB correlated genes: 888
  m = length(unique(df_pcorRxn2BrainRNATop$Gene)), 
  # total (background) - MB Correlated genes: 31406
  n = nrow(df_brainRNAnormFiltered)-length(unique(df_pcorRxn2BrainRNATop$Gene)), 
  # Age assoc. genes: 6505
  k = sum( df_brainAgeDependance$padj <= 0.05, na.rm = T ),
  lower.tail = F)
# probablity of seing an equal or bigger overlap:
# [1] 3.685676e-05
# There is a significant enrichment of MB-associated Genes in Age-associated Genes in Brain

# export results table
df_tmp <- merge(x = df_brainAgeDependance[df_brainAgeDependance$padj <= 0.05,c("log2FoldChange"),drop = F],
                y = df_pcorRxn2BrainRNATop[,c("Gene","Rxn","Spearman_rho","GeneSymbol")],
                by.x = 0, by.y = "Gene")
colnames(df_tmp) <- c("Host_EnsemblID", "Age_log2FC", "MB_ReactionID", "MB_SpearmanRho", "Host_GeneSymbol")
df_tmp <- df_tmp[,c("Host_EnsemblID", "Host_GeneSymbol", "Age_log2FC", "MB_ReactionID", "MB_SpearmanRho")]
export(df_tmp, filename = "mtx_BrainSharedAgeMBAssocGenes.csv")

#####
## age pos. log2FC (up in old) and MB Correlations Genes 
obj_Gene2GO <- enricher(gene = unique(na.omit(df_mm10Ensembl2Genesymbol[ intersect(unique(df_pcorRxn2BrainRNATop$Gene),
                                                                                   na.omit(rownames(df_brainAgeDependance)[df_brainAgeDependance$padj <= 0.05 &
                                                                                                                           df_brainAgeDependance$log2FoldChange > 0]) ),
                                                                         "GeneSymbol"])),
                        minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                        TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                        TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                        universe = unique(na.omit(df_mm10Ensembl2Genesymbol[rownames(df_brainRNAnormFiltered),"GeneSymbol"])))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GO <- obj_Gene2GO@result[order(obj_Gene2GO@result$qvalue,
                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GO@result$ID], 
                                       decreasing = c(F,F)),]
df_Gene2GOPositive <- df_Gene2GO[df_Gene2GO$p.adjust <= 0.1 &   
                                   df_Gene2GO$Count >= 3,] ## p.adj
# df_Gene2GOPositive <- df_trimGO(df_GOresultTable = df_Gene2GOPositive)

#####
## age neg. log2FC (up in old) and MB Correlations Genes 
obj_Gene2GO <- enricher(gene = unique(na.omit(df_mm10Ensembl2Genesymbol[ intersect(unique(df_pcorRxn2BrainRNATop$Gene),
                                                                                   na.omit(rownames(df_brainAgeDependance)[df_brainAgeDependance$padj <= 0.05 &
                                                                                                                           df_brainAgeDependance$log2FoldChange < 0]) ),
                                                                         "GeneSymbol"])),
                        minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                        TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                        TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                        universe = unique(na.omit(df_mm10Ensembl2Genesymbol[rownames(df_brainRNAnormFiltered),"GeneSymbol"])))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GO <- obj_Gene2GO@result[order(obj_Gene2GO@result$qvalue,
                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GO@result$ID], 
                                       decreasing = c(F,F)),]
df_Gene2GONegative <- df_Gene2GO[df_Gene2GO$p.adjust <= 0.1 &
                                   df_Gene2GO$Count >= 3,]
df_Gene2GONegative <- df_trimGO(df_GOresultTable = df_Gene2GONegative)
# combine
df_Gene2GO <- rbind(df_Gene2GOPositive, #[!df_Gene2GOPositive$softTrim,]
                    df_Gene2GONegative[!df_Gene2GONegative$softTrim,])

# add MB correlation
df_Gene2GO$MedianMBcor <- apply(X = df_Gene2GO, MARGIN = 1, FUN = function(obj_row) {
  lst_genes <-unique( unlist(strsplit(split = "/", fixed = T, x = obj_row[["geneID"]])) )
  lst_genes <- df_mm10Ensembl2Genesymbol[df_mm10Ensembl2Genesymbol$GeneSymbol %in% lst_genes, "EnsemblID"]
  # collect median MB-correlation for each gene
  df_tmp <- df_pcorRxn2BrainRNATop[df_pcorRxn2BrainRNATop$Gene %in% lst_genes,]
  df_tmp <- tapply(X = df_tmp$Spearman_rho, INDEX = df_tmp$Gene, FUN = median)
  # summarize median correlations for all genes
  return( median(df_tmp) )
})
## add Age foldChange
df_Gene2GO$MedianAgeLog2FC <- apply(X = df_Gene2GO, MARGIN = 1, FUN = function(obj_row) {
  lst_genes <-unique( unlist(strsplit(split = "/", fixed = T, x = obj_row[["geneID"]])) )
  lst_genes <- df_mm10Ensembl2Genesymbol[df_mm10Ensembl2Genesymbol$GeneSymbol %in% lst_genes, "EnsemblID"]
  # summarize age log2FC for all genes with median 
  return( median(df_DifAbundBrainAllAges[lst_genes,"log2FoldChange"]) )  #log2FoldChange
})
#
df_Gene2GO$DescriptionSplit <- str_simplifyDescriptions(df_Gene2GO$Description, bol_stripRomanNums = F)
df_Gene2GO <- df_Gene2GO[order(df_Gene2GO$MedianAgeLog2FC),]

pdf(file = "GOplot_brainAgeAndMBassociatedGenesGO.pdf", width = 290/72, height = 90/72, pointsize = 8)
par(mar = c(4,22,1,1)+.1)
plot(x = df_Gene2GO$MedianAgeLog2FC,
     xlab = "Age log2(FC)", xlim = c(-.1,.1), xaxt = "n", lwd = .8,
     y = 1:nrow(df_Gene2GO), ylab = "", yaxt = "n", ylim = c(0,nrow(df_Gene2GO)+1),
     pch = 21, cex = scales::rescale(-log10(df_Gene2GO$pvalue), 
                                     to = c(1,3)),
     bg = brewer.pal("RdBu", n= 11)[11:1][round(scales::rescale(df_Gene2GO$MedianMBcor, from = c(-1,1), to = c(1,11)))],
     main = "Brain"
)
axis(side = 2, at = 1:nrow(df_Gene2GO), labels = Hmisc::capitalize( df_Gene2GO$DescriptionSplit ), las = 2, cex.axis = 1.1)
axis(side = 1, at = c(-1, 0, 1)/10)
abline(v = 0, col = rgb(.3,.3,.3,.7), lty = 4, lwd = 1)
dev.off()
# export at 600x260 GOplot_brainAgeAndMBassociatedGenesGO

# export GO table for supplements
export(df_Gene2GO[,c(1:6,8,9,12,13)], filename = "mtx_BrainSharedAgeMBGOPathways.csv")


###################
### Liver
## Hypergeometric test for enrichment of overlap in age-associated to MB-associated liver transcripts
phyper(# Shared genes (overlap): 500-1
  q = length(intersect(unique(df_pcorRxn2LiverRNATop$Gene),
                       na.omit(rownames(df_liverAgeDependance)[df_liverAgeDependance$padj <= 0.05]))) - 1, 
  # MB correlated genes: 1265
  m = length(unique(df_pcorRxn2LiverRNATop$Gene)), 
  # total (background) - MB Correlated genes: 24796
  n = nrow(df_liverRNAnormFiltered)-length(unique(df_pcorRxn2LiverRNATop$Gene)), 
  # Age assoc. genes: 8285
  k = sum( df_liverAgeDependance$padj <= 0.05, na.rm = T ),
  lower.tail = F)
# probablity of seing an equal or bigger overlap:
# [1] 9.088653e-07
# There is a significant enrichment of MB-associated Genes in Age-associated Genes in Liver

# export results table
df_tmp <- merge(x = df_liverAgeDependance[df_liverAgeDependance$padj <= 0.05,c("log2FoldChange"),drop = F],
                y = df_pcorRxn2LiverRNATop[,c("Gene","Rxn","Spearman_rho","GeneSymbol")],
                by.x = 0, by.y = "Gene")
colnames(df_tmp) <- c("Host_EnsemblID", "Age_log2FC", "MB_ReactionID", "MB_SpearmanRho", "Host_GeneSymbol")
df_tmp <- df_tmp[,c("Host_EnsemblID", "Host_GeneSymbol", "Age_log2FC", "MB_ReactionID", "MB_SpearmanRho")]
export(df_tmp, filename = "mtx_LiverSharedAgeMBAssocGenes.csv")

#####
## age pos. log2FC (up in old) and MB Correlations Genes 
obj_Gene2GO <- enricher(gene = unique(na.omit(df_mm10Ensembl2Genesymbol[ intersect(unique(df_pcorRxn2LiverRNATop$Gene),
                                                                                   na.omit(rownames(df_liverAgeDependance)[df_liverAgeDependance$padj <= 0.05 &
                                                                                                                           df_liverAgeDependance$log2FoldChange > 0]) ),
                                                                         "GeneSymbol"])),
                        minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                        TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                        TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                        universe = unique(na.omit(df_mm10Ensembl2Genesymbol[rownames(df_liverRNAnormFiltered),"GeneSymbol"])))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GO <- obj_Gene2GO@result[order(obj_Gene2GO@result$qvalue,
                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GO@result$ID], 
                                       decreasing = c(F,F)),]
df_Gene2GOPositive <- df_Gene2GO[df_Gene2GO$p.adjust <= 0.1 &
                                   df_Gene2GO$Count >= 3,]
df_Gene2GOPositive <- df_trimGO(df_GOresultTable = df_Gene2GOPositive)

#####
## age neg. log2FC (up in old) and MB Correlations Genes 
obj_Gene2GO <- enricher(gene = unique(na.omit(df_mm10Ensembl2Genesymbol[ intersect(unique(df_pcorRxn2LiverRNATop$Gene),
                                                                                   na.omit(rownames(df_liverAgeDependance)[df_liverAgeDependance$padj <= 0.05 &
                                                                                                                           df_liverAgeDependance$log2FoldChange < 0]) ),
                                                                         "GeneSymbol"])),
                        minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH",
                        TERM2GENE = df_mouseGOannotaionBioProc[,c("GOid","GeneSymbol")],
                        TERM2NAME = df_mouseGOannotaionBioProc[,c("GOid","Description")],
                        universe = unique(na.omit(df_mm10Ensembl2Genesymbol[rownames(df_liverRNAnormFiltered),"GeneSymbol"])))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Gene2GO <- obj_Gene2GO@result[order(obj_Gene2GO@result$qvalue,
                                       table(df_mouseGOannotaionBioProc$GOid)[obj_Gene2GO@result$ID], 
                                       decreasing = c(F,F)),]
df_Gene2GONegative <- df_Gene2GO[df_Gene2GO$p.adjust <= 0.1 &
                                   df_Gene2GO$Count >= 3,]
# df_Gene2GONegative <- df_trimGO(df_GOresultTable = df_Gene2GONegative)

# combine
df_Gene2GO <- rbind(df_Gene2GOPositive[!df_Gene2GOPositive$softTrim,],
                    df_Gene2GONegative)
df_Gene2GOLiver <- df_Gene2GO
# add MB correlation
df_Gene2GO$MedianMBcor <- apply(X = df_Gene2GO, MARGIN = 1, FUN = function(obj_row) {
  lst_genes <-unique( unlist(strsplit(split = "/", fixed = T, x = obj_row[["geneID"]])) )
  lst_genes <- df_mm10Ensembl2Genesymbol[df_mm10Ensembl2Genesymbol$GeneSymbol %in% lst_genes, "EnsemblID"]
  # collect median MB-correlation for each gene
  df_tmp <- df_pcorRxn2LiverRNATop[df_pcorRxn2LiverRNATop$Gene %in% lst_genes,]
  df_tmp <- tapply(X = df_tmp$Spearman_rho, INDEX = df_tmp$Gene, FUN = median)
  # summarize median correlations for all genes
  return( median(df_tmp) )
})
# add Median Age foldChange
df_Gene2GO$MedianAgeLog2FC <- apply(X = df_Gene2GO, MARGIN = 1, FUN = function(obj_row) {
  lst_genes <-unique( unlist(strsplit(split = "/", fixed = T, x = obj_row[["geneID"]])) )
  lst_genes <- df_mm10Ensembl2Genesymbol[df_mm10Ensembl2Genesymbol$GeneSymbol %in% lst_genes, "EnsemblID"]
  # summarize age log2FC for all genes with median
  return( median(df_DifAbundLiverAllAges[lst_genes,"log2FoldChange"]) )
})
#
df_Gene2GO$DescriptionSplit <- str_simplifyDescriptions(df_Gene2GO$Description, bol_stripRomanNums = F)
df_Gene2GO <- df_Gene2GO[order(df_Gene2GO$MedianAgeLog2FC),]
#
df_Gene2GO$DescriptionSplit[5] <- "reg. of force of heart contraction"

pdf(file = "GOplot_liverAgeAndMBassociatedGenesGO.pdf", width = 290/72, height = 180/72, pointsize = 8)
par(mar = c(4,22,1,1)+.1)
plot(x = df_Gene2GO$MedianAgeLog2FC,
     xlab = "Age log2(FC)", xlim = c(-.3,1.3), lwd = .8,
     y = 1:nrow(df_Gene2GO), ylab = "", yaxt = "n", xaxt = "n",
     pch = 21, cex = scales::rescale(-log10(df_Gene2GO$pvalue), 
                                     to = c(1,3)),
     bg = brewer.pal("RdBu", n= 11)[11:1][round(scales::rescale(df_Gene2GO$MedianMBcor, from = c(-1,1), to = c(1,11)))],
     main = "Liver"
)
axis(side = 2, at = 1:nrow(df_Gene2GO), labels = Hmisc::capitalize( df_Gene2GO$DescriptionSplit ), las = 2, cex.axis = 1.1)
axis(side = 1, at = c(-2,0,4,8,12)/10)
abline(v = 0, col = rgb(.3,.3,.3,.7), lty = 4, lwd = 1)
dev.off()
# legend("topleft", legend = c("pos. MB cor.", "neg. MB cor."), 
#        col = c("#2166AC", "#B2182B"), pch = 16, bty = "n", cex = .7)
# export at 600x700 GOplot_liverAgeAndMBassociatedGenesGO

# export GO table for supplements
export(df_Gene2GO[,c(1:6,8,9,12,13)], filename = "mtx_LiverSharedAgeMBGOPathways.csv")


###################################
## MB side:
## single reactions ageing association
# df_MMRxnsByAge <- readRDS("../mouseMetagenome/df_MMRxnsByAge.rds")
rownames(df_MMRxnsByAge) <- df_MMRxnsByAge$RXN

###################
### All organs combined - any features added up from any organ (UNION)
## Hypergeometric test for enrichment of overlap in age-associated with Host-associated MB-reactions
phyper(# Shared genes (overlap): 445-1
  q = length( intersect(unique(c(df_pcorRxn2ColonRNATop$Rxn,
                                 df_pcorRxn2LiverRNATop$Rxn,
                                 df_pcorRxn2BrainRNATop$Rxn)),
                        df_MMRxnsByAge[df_MMRxnsByAge$p.adj <= 0.05,"RXN"]) ) - 1, 
  # host correlated features: 2105
  m = length( unique(c(df_pcorRxn2ColonRNATop$Rxn,
                       df_pcorRxn2LiverRNATop$Rxn,
                       df_pcorRxn2BrainRNATop$Rxn)) ),
  # total (background) - MB Correlated genes: 105 
  n = nrow(df_MMRxnsByAge)-length( unique(c(df_pcorRxn2ColonRNATop$Rxn,
                                            df_pcorRxn2LiverRNATop$Rxn,
                                            df_pcorRxn2BrainRNATop$Rxn)) ), 
  # Age assoc. genes: 466
  k = sum( df_MMRxnsByAge$p.adj <= 0.05, na.rm = T ),
  lower.tail = F)
# probablity of seing an equal or bigger overlap:
# [1] 0.5101777
# There is no significant enrichment of Host-associated microbiome functions in Age-associated microbiome functions

# export results table
df_tmp <- merge(x = df_MMRxnsByAge[df_MMRxnsByAge$p.adj <= 0.05,c("log2FC"),drop = F],
                y = rbind(df_pcorRxn2ColonRNATop,
                          df_pcorRxn2LiverRNATop,
                          df_pcorRxn2BrainRNATop)[,c("Gene","Rxn","Spearman_rho","GeneSymbol")],
                by.x = 0, by.y = "Rxn")
colnames(df_tmp) <- c("MB_ReactionID", "MB_Age_log2FC", "Host_EnsemblID", "Host_SpearmanRho", "Host_GeneSymbol")
df_tmp <- df_tmp[,c("MB_ReactionID", "MB_Age_log2FC", "Host_EnsemblID", "Host_GeneSymbol", "Host_SpearmanRho")]
export(df_tmp, filename = "mtx_MicrobiomeSharedAgeHostAssocRxns.csv")

### Ontology Enrichment for significant correlations only!
## Positive Correlations Reactions
obj_Rxn2Subsys <- enricher(gene = unique(na.omit( intersect(unique(c(df_pcorRxn2ColonRNATop$Rxn,
                                                                     df_pcorRxn2LiverRNATop$Rxn,
                                                                     df_pcorRxn2BrainRNATop$Rxn)),
                                                            df_MMRxnsByAge[df_MMRxnsByAge$p.adj <= 0.05 &
                                                                           df_MMRxnsByAge$log2FC > 0,"RXN"]) )),
                           minGSSize = 3, maxGSSize = 500, pAdjustMethod = "BH",
                           TERM2GENE = df_rxn2subsys[,c(1,2)],
                           TERM2NAME = df_rxn2subsys[,c(1,3)],
                           universe = rownames(df_colonMMrxnFiltered))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Rxn2SubsysPositive <- obj_Rxn2Subsys@result[order(obj_Rxn2Subsys@result$qvalue,
                                                     decreasing = F),]
df_Rxn2SubsysPositive <- df_Rxn2SubsysPositive[df_Rxn2SubsysPositive$p.adjust <= 0.05 &
                                                 df_Rxn2SubsysPositive$Count >= 3,]
df_checkGO(df_Rxn2SubsysPositive$ID, df_Rxn2SubsysPositive)
# drop duplicated subsystems due to identical rxns in enrichment
# [1] "Subsys:357" "Subsys:363"
df_Rxn2SubsysPositive <- df_GOresultsReference; rm(df_GOresultsReference)

#####
## Colon age neg. log2FC (up in old) and MB Correlations Genes 
obj_Rxn2Subsys <- enricher(gene = unique(na.omit( intersect(unique(c(df_pcorRxn2ColonRNATop$Rxn,
                                                                      df_pcorRxn2LiverRNATop$Rxn,
                                                                      df_pcorRxn2BrainRNATop$Rxn)),
                                                             df_MMRxnsByAge[df_MMRxnsByAge$p.adj <= 0.05 &
                                                                              df_MMRxnsByAge$log2FC < 0,"RXN"]) )),
                           minGSSize = 3, maxGSSize = 500, pAdjustMethod = "BH",
                           TERM2GENE = df_rxn2subsys[,c(1,2)],
                           TERM2NAME = df_rxn2subsys[,c(1,3)],
                           universe = rownames(df_colonMMrxnFiltered))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Rxn2SubsysNegative <- obj_Rxn2Subsys@result[order(obj_Rxn2Subsys@result$qvalue,
                                                     decreasing = F),]
df_Rxn2SubsysNegative <- df_Rxn2SubsysNegative[df_Rxn2SubsysNegative$p.adjust <= 0.05 &
                                                 df_Rxn2SubsysNegative$Count >= 3,]
df_checkGO(df_Rxn2SubsysNegative$ID, df_Rxn2SubsysNegative)
# [1] "Subsys:871" "Subsys:812"
df_Rxn2SubsysNegative <- df_GOresultsReference; rm(df_GOresultsReference)

# combine
df_Gene2GO <- rbind(df_Rxn2SubsysPositive,
                    df_Rxn2SubsysNegative)
# update descriptions
df_Gene2GO$SubsystemSplit <- str_simplifyDescriptions(df_Gene2GO$Description)
df_Gene2GO$SubsystemSplit <- Hmisc::capitalize( df_Gene2GO$SubsystemSplit )
df_Gene2GO$SubsystemSplit[4] <- "Folate synth. ..."
df_Gene2GO$SubsystemSplit[18] <- "ppGpp synth."

# merge similar terms
for(str_subsys in df_Gene2GO$SubsystemSplit) {
  if(sum(df_Gene2GO$SubsystemSplit == str_subsys) >= 2) {
    # replace genelists by shared (union) gene lists
    df_Gene2GO[df_Gene2GO$SubsystemSplit == str_subsys,"geneID"] <- paste(collapse = "/", 
                                                                          unique(unlist( strsplit(split = "/", fixed = T, 
                                                                                                  x = df_Gene2GO[df_Gene2GO$SubsystemSplit == str_subsys,"geneID"]) )))
    # update feature count to max of shared terms
    df_Gene2GO[df_Gene2GO$SubsystemSplit == str_subsys,"Count"] <- max(df_Gene2GO[df_Gene2GO$SubsystemSplit == str_subsys,"Count"])
  }
}
# order by p-val and drop duplicated ones
df_Gene2GO <- df_Gene2GO[order(df_Gene2GO$pvalue, decreasing = F),]
df_Gene2GO <- df_Gene2GO[!duplicated(df_Gene2GO$SubsystemSplit),]

# add Host correlation
df_Gene2GO$MedianMBcor <- apply(X = df_Gene2GO, MARGIN = 1, FUN = function(obj_row) {
  lst_genes <-unique( unlist(strsplit(split = "/", fixed = T, x = obj_row[["geneID"]])) )
  # lst_genes <- df_mm10Ensembl2Genesymbol[df_mm10Ensembl2Genesymbol$GeneSymbol %in% lst_genes, "EnsemblID"]
  # collect median MB-correlation for each gene
  df_tmp <- rbind(df_pcorRxn2ColonRNATop[df_pcorRxn2ColonRNATop$Rxn %in% lst_genes,],
                  df_pcorRxn2LiverRNATop[df_pcorRxn2LiverRNATop$Rxn %in% lst_genes,],
                  df_pcorRxn2BrainRNATop[df_pcorRxn2BrainRNATop$Rxn %in% lst_genes,])
  df_tmp <- tapply(X = df_tmp$Spearman_rho, INDEX = df_tmp$Rxn, FUN = median)
  # summarize median correlations for all genes
  return( median(df_tmp) )
})
# add Age foldChange
df_Gene2GO$MedianAgeLog2FC <- apply(X = df_Gene2GO, MARGIN = 1, FUN = function(obj_row) {
  lst_genes <-unique( unlist(strsplit(split = "/", fixed = T, x = obj_row[["geneID"]])) )
  # summarize age log2FC for all genes with median 
  return( median(df_MMRxnsByAge[lst_genes,"log2FC"]) )
})
# orde by age log2FC
df_Gene2GO <- df_Gene2GO[order(df_Gene2GO$MedianAgeLog2FC),]

# plot it 
pdf(file = "GOplot_AgeAndHostAssociatedRxnsSubsys.pdf", width = 290/72, height = 260/72, pointsize = 8)
par(mar = c(4,1,1,22)+.1)
plot(x = df_Gene2GO$MedianAgeLog2FC,
     xlab = "Age log2(FC)", xlim = c(-1.3,4), lwd = .8,
     y = 1:nrow(df_Gene2GO), ylab = "", yaxt = "n",
     pch = 21, cex = scales::rescale(-log10(df_Gene2GO$pvalue), 
                                     to = c(1,3)),
     bg = brewer.pal("RdBu", n= 11)[11:1][round(scales::rescale(df_Gene2GO$MedianMBcor, from = c(-1,1), to = c(1,11)))],
     main = "Microbiome"
)
axis(side = 4, at = 1:nrow(df_Gene2GO), labels = df_Gene2GO$SubsystemSplit, las = 2, cex.axis = 1.1)
abline(v = 0, col = rgb(.3,.3,.3,.7), lty = 4, lwd = 1)
dev.off()
# legend("bottomright", legend = c("pos. MB cor.", "neg. MB cor.", "size = enrich. p-val."),
#        pt.bg = c("#2166AC", "#B2182B"), pch = c(21,21,NA), bty = "n", cex = .7)

# horizontal version
pdf(file = "GOplot_AgeAndHostAssociatedRxnsSubsys_horiz.pdf", width = 320/72, height = 220/72, pointsize = 8)
par(mar = c(10,6,1,1)+.1)
plot(y = df_Gene2GO$MedianAgeLog2FC,
     ylab = "Age log2(FC)", ylim = c(-1.3,4), lwd = .6,
     x = 1:nrow(df_Gene2GO), xlab = "",
     pch = 21, cex = scales::rescale(-log10(df_Gene2GO$pvalue), 
                                     to = c(1,3)),
     bg = brewer.pal("RdBu", n= 11)[11:1][round(scales::rescale(df_Gene2GO$MedianMBcor, from = c(-1,1), to = c(1,11)))],
     axes = F,
     # main = "Microbiome"
)
axis(side = 2, at = c(-1:4), labels = c(-1,0,NA,2,NA,4), 
     lwd = .6, lwd.ticks = .6, cex.axis = 1, las = 2)
axis(side = 1, at = 1:nrow(df_Gene2GO), labels = NA, lwd.ticks = .6, lwd = .6)
text(x = 1:nrow(df_Gene2GO), y = -2, labels = df_Gene2GO$SubsystemSplit, cex = 1,
     srt = 40, adj = c(1,.5), xpd = T)
abline(h = 0, col = rgb(.3,.3,.3,.7), lty = 4, lwd = .6)
dev.off()

# export SubSys enrichment table for supplements
export(df_Gene2GO[,c(1:6,8,9,11,12)], filename = "mtx_MicrobiomeSharedAgeHostMetaCycPathways.csv")


############################

