##############################################################
##
##    Metagenome Metabolic Models FVA of ageing mice 
##
##############################################################

###
library(RColorBrewer)
library(clusterProfiler)

source("customFunctions.R")

setwd("mouseMetagenome/")


################
# get sample metadata
df_metaDataMM <- readRDS(file = "df_metaData.rds")
rownames(df_metaDataMM) <- df_metaDataMM$SampleName


##############
# load abundance data
df_rxnAbundancesMM <- readRDS(file = "./df_rxnAbundancesMM.rds")
# data is already normalized
colSums(df_rxnAbundancesMM)
df_rxnAbundancesMM <- df_rxnAbundancesMM * 1000000
# relabel columns to sample names
df_rxnAbundancesMM <- df_rxnAbundancesMM[,df_metaDataMM$SequencingID]
colnames(df_rxnAbundancesMM) <- df_metaDataMM$SampleName


#################
## Evaluate data
void_MDSplot(df_sampleMatrix = df_rxnAbundancesMM[,df_metaDataMM$SampleName], df_metaData = df_metaDataMM, bol_labels = F, 
             var_colorBy = "Age", var_shapeBy = NULL, bol_continousColor = F)


#########################################
#
df_lmMetabolicModelsRxnAbundances <- df_linearModel(df_dependantVariable = df_rxnAbundancesMM[,df_metaDataMM$SampleName],
                                                    lst_var_metaData = df_metaDataMM$Age,
                                                    int_cores = 7)
# [1] "Starting rank sum tests of 2129 observations within 52 samples! Wed May 17 16:32:16 2023"
# [1] "Finished all tests! Preparing results. Wed May 17 16:32:16 2023"
saveRDS(df_lmMetabolicModelsRxnAbundances,paste0("df_lmMetabolicModelsRxnAbundances_", Sys.Date(),".rds"))
## few diagnostic plots
# p-val Histogram
hist(df_lmMetabolicModelsRxnAbundances$p.value, breaks = 40)

# volcano plot
plot(x = df_lmMetabolicModelsRxnAbundances$est,
     y = -log10(df_lmMetabolicModelsRxnAbundances$p.value),
     pch = 4, cex = .7, #col = as.numeric(as.factor(df_tmp$ModelType)),
     xlab = "Regression Coefficient of Treatment", ylab = "-log10( p-value )",
     main = "",# xlim = c(-6,6)
); abline(h = -log10(0.05), v = c(-.1,.1), col = 2, lty = 2)

# PCoA
void_MDSplot(df_sampleMatrix = df_rxnAbundancesMM[rownames(df_lmMetabolicModelsRxnAbundances)[df_lmMetabolicModelsRxnAbundances$p.adj <= 0.05],
                                             df_metaDataMM$SampleName], 
             df_metaData = df_metaDataMM, bol_labels = F, 
             var_colorBy = "Age", var_shapeBy = NULL, bol_continousColor = F)

# Sort by p-vals
df_lmMetabolicModelsRxnAbundances <- as.data.frame(df_lmMetabolicModelsRxnAbundances[order(df_lmMetabolicModelsRxnAbundances$p.value, decreasing = F),])
df_lmMetabolicModelsRxnAbundances$RXN <- rownames(df_lmMetabolicModelsRxnAbundances)
rownames(df_lmMetabolicModelsRxnAbundances) <- NULL

# Add Fold Change of young to oldest  log2( FC ) = log2( mean(Outcome) / mean(Baseline) )
df_lmMetabolicModelsRxnAbundances$log2FC <- log2( apply(X = df_rxnAbundancesMM[df_lmMetabolicModelsRxnAbundances$RXN,
                                                                  df_metaDataMM[df_metaDataMM$Age == 30,"SampleName"]], 
                                                MARGIN = 1, FUN = mean) /
                                            apply(X = df_rxnAbundancesMM[df_lmMetabolicModelsRxnAbundances$RXN,
                                                                    df_metaDataMM[df_metaDataMM$Age == 2,"SampleName"]],
                                                  MARGIN = 1, FUN = mean) )

# volcano plot
plot(x = df_lmMetabolicModelsRxnAbundances$log2FC,
     y = -log10(df_lmMetabolicModelsRxnAbundances$p.value),
     pch = 4, cex = .7, #col = as.numeric(as.factor(df_tmp$ModelType)),
     xlab = "log2( FC ) 2 -> 30 month", ylab = "-log10( p-value )",
     main = "",
); abline(h = -log10(0.05), v = c(-.1,.1), col = 2, lty = 2)

# log2FC is negative -> Decreased from 2 month to 30 month
boxplot(as.numeric( df_rxnAbundancesMM[df_lmMetabolicModelsRxnAbundances[6,"RXN"],df_metaDataMM$SampleName] ) ~ 
          df_metaDataMM$Age); df_lmMetabolicModelsRxnAbundances[6,]
# log2FC is positive -> Increased from 2 month to 30 month
boxplot(as.numeric( df_rxnAbundancesMM[df_lmMetabolicModelsRxnAbundances[1,"RXN"],df_metaDataMM$SampleName] ) ~ 
          df_metaDataMM$Age); df_lmMetabolicModelsRxnAbundances[1,]
#####
# store for later
write.table(x = df_lmMetabolicModelsRxnAbundances[df_lmMetabolicModelsRxnAbundances$p.adj <= 0.05 &
                                                  df_lmMetabolicModelsRxnAbundances$log2FC < 0,"RXN"], 
            file = "AgeingDependentReactions_DownIn30m.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = df_lmMetabolicModelsRxnAbundances[df_lmMetabolicModelsRxnAbundances$p.adj <= 0.05 &
                                                  df_lmMetabolicModelsRxnAbundances$log2FC > 0,"RXN"], 
            file = "AgeingDependentReactions_UpIn30m.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

saveRDS(df_lmMetabolicModelsRxnAbundances, "df_MMRxnsByAge.rds")




#############################################
## Full reaction annotation for all dif. abundant reactions not only heatmap plotted ones
obj_Rxn2SubsysUp <- enricher(gene = df_lmMetabolicModelsRxnAbundances[df_lmMetabolicModelsRxnAbundances$p.adj <= 0.05 &
                                                                      df_lmMetabolicModelsRxnAbundances$est > 0, "RXN"],
                             minGSSize = 3, maxGSSize = 500, pAdjustMethod = "BH",
                             TERM2GENE = df_rxn2subsys[,c(1,2)],
                             TERM2NAME = df_rxn2subsys[,c(1,3)],
                             universe = df_lmMetabolicModelsRxnAbundances$RXN)
df_Rxn2SubsysUp <- obj_Rxn2SubsysUp@result[obj_Rxn2SubsysUp@result$p.adjust <= 0.05 &
                                           obj_Rxn2SubsysUp@result$Count >= 3,]
# # Go through Bacterial Subsystems
# lst_dropRows <- df_checkGO(lst_checkGOIDs = rownames(df_Rxn2SubsysUp),
#                            df_GOresultsReference = df_Rxn2SubsysUp )
# # drop rows
# # [1] "Subsys:172" "Subsys:355"
# df_Rxn2SubsysUp <- df_Rxn2SubsysUp[!df_Rxn2SubsysUp$ID %in% lst_dropRows,]
# df_Rxn2SubsysUp$Description <- df_GOresultsReference[df_Rxn2SubsysUp$ID,"Description"]
View(df_Rxn2SubsysUp)
#
write.table(x = df_Rxn2SubsysUp, 
            file = "RxnUpInAgeingSubsystemEnrichment.txt", 
            quote = T, sep = "\t", row.names = F, col.names = T)
saveRDS(df_Rxn2SubsysUp, "RxnUpInAgeingMiceSubsystemEnrichment.rds")

## Down regulated in ageing
obj_Rxn2SubsysDown <- enricher(gene = df_lmMetabolicModelsRxnAbundances[df_lmMetabolicModelsRxnAbundances$p.adj <= 0.05 &
                                                                        df_lmMetabolicModelsRxnAbundances$est < 0, "RXN"],
                               minGSSize = 3, maxGSSize = 500, pAdjustMethod = "BH",
                               TERM2GENE = df_rxn2subsys[,c(1,2)],
                               TERM2NAME = df_rxn2subsys[,c(1,3)],
                               universe = df_lmMetabolicModelsRxnAbundances$RXN)
df_Rxn2SubsysDown <- obj_Rxn2SubsysDown@result[obj_Rxn2SubsysDown@result$p.adjust <= 0.05 &
                                               obj_Rxn2SubsysDown@result$Count >= 3,]
# # Go through Bacterial Subsystems
# lst_dropRows <- df_checkGO(lst_checkGOIDs = rownames(df_Rxn2SubsysDown),
#                            df_GOresultsReference = df_Rxn2SubsysDown )
# # drop rows
# # [1] "Subsys:780"
# df_Rxn2SubsysDown <- df_Rxn2SubsysDown[!df_Rxn2SubsysDown$ID %in% lst_dropRows,]
# df_Rxn2SubsysDown$Description <- df_GOresultsReference[df_Rxn2SubsysDown$ID,"Description"]

View(df_Rxn2SubsysDown)
#
write.table(x = df_Rxn2SubsysDown, 
            file = "RxnDownInAgeingSubsystemEnrichment.txt", 
            quote = T, sep = "\t", row.names = F, col.names = T)
saveRDS(df_Rxn2SubsysDown, "RxnDownInAgeingMiceSubsystemEnrichment.rds")
saveRDS(df_lmMetabolicModelsRxnAbundances, "df_mouseLmMetabolicModelsRxnByAge.rds")


################ 
## Make a heatmap of mean rxn activity from each GO-group 
df_Rxn2Subsys <- rbind(df_Rxn2SubsysUp[df_Rxn2SubsysUp$p.adjust <= 0.05 &
                                       df_Rxn2SubsysUp$Count >= 3,],
                       df_Rxn2SubsysDown[df_Rxn2SubsysDown$p.adjust <= 0.05 &
                                         df_Rxn2SubsysDown$Count >= 3,]
                       )
# mtx_mydist(df_Rxn2SubsysNegative)
mtx_tmp2 <- matrix(data = NA, nrow = nrow(df_Rxn2Subsys), ncol = 5)
rownames(mtx_tmp2) <- df_Rxn2Subsys[hclust(as.dist(mtx_mydist(df_Rxn2Subsys)))$order,"Description"]
for(str_GOId in df_Rxn2Subsys$ID) {
  lst_rxns <- unlist(strsplit(x = df_Rxn2Subsys[str_GOId,"geneID"], "/", fixed = T))
  mtx_tmp <- df_rxnAbundancesMM[df_lmMetabolicModelsRxnAbundances[df_lmMetabolicModelsRxnAbundances$p.adj <= 0.05,"RXN"],
                                   df_metaDataMM$SampleName]
  mtx_tmp <- as.matrix(mtx_tmp[lst_rxns,order(df_metaDataMM[colnames(mtx_tmp),"AgeGroup"])])

  # normalise matrix so each reaction is similarily abundant (percentage of rxn-abundance)
  mtx_tmp <- mtx_tmp / rowSums(mtx_tmp) * 1000
  
  mtx_tmp2[df_Rxn2Subsys[str_GOId,"Description"],] <- c( mean( mtx_tmp[, df_metaDataMM[colnames(mtx_tmp),"AgeGroup"] == 1 ] ),
                                                         mean( mtx_tmp[, df_metaDataMM[colnames(mtx_tmp),"AgeGroup"] == 2 ] ),
                                                         mean( mtx_tmp[, df_metaDataMM[colnames(mtx_tmp),"AgeGroup"] == 3 ] ),
                                                         mean( mtx_tmp[, df_metaDataMM[colnames(mtx_tmp),"AgeGroup"] == 4 ] ),
                                                         mean( mtx_tmp[, df_metaDataMM[colnames(mtx_tmp),"AgeGroup"] == 5 ] )
  )
}

str_rowLabels <- str_simplifyDescriptions( rownames(mtx_tmp2) )
str_rowLabels <- sub(pattern = "ppGpp", replacement = "Guanosine pentaphosphate", x = str_rowLabels, fixed = T)
str_rowLabels <- Hmisc::capitalize(str_rowLabels)

rownames(mtx_tmp2) <- str_rowLabels
# drop duplicates by value
duplicated(mtx_tmp2)
mtx_tmp2 <- mtx_tmp2[rownames(mtx_tmp2) != "Cis-genanyl-CoA degr.",] # is a duplicate of "Acetate conversion to acetyl-CoA "
mtx_tmp2 <- mtx_tmp2[!duplicated(mtx_tmp2),]
# drop duplicates by name
rownames(mtx_tmp2)[duplicated(rownames(mtx_tmp2))]
mtx_tmp2 <- mtx_tmp2[!duplicated(rownames(mtx_tmp2)),]

# create a heatmap plot
pdf(file = "heatmap_MouseMeanSubSysAbundanceByAge.pdf", width = 700/72, height = 1100/72, pointsize = 14)
gplots::heatmap.2(x = mtx_tmp2,
                  dendrogram = "none",
                  Colv = NA,
                  Rowv = NA,
                  scale = "row",
                  ColSideColors = brewer.pal("Pastel1", n = 5)[1:5],
                  col = colorRampPalette(brewer.pal("Purples", n=9))(50),
                  colsep = F, rowsep = F, 
                  key = T, trace = "n", density.info = "n", 
                  cexRow = 1.6,
                  labRow = rownames(mtx_tmp2),
                  # srtCol = 1, 
                  adjCol = c(1,.5),
                  labCol = c(2,9,15,24,30),
                  margins = c(5,32)
)
dev.off()
# export at 700x1100 heatmap_MouseMeanSubSysAbundanceByAge
