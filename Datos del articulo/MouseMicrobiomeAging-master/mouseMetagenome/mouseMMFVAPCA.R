##################################################
##
##  Model differences by PCA
##
##################################################


library(vegan)
library(RColorBrewer)
library(clusterProfiler)

source("customFunctions.R")
setwd("mouseMetagenome/")


####################
## Load Models for Mice
obj_metabolicModelsMouse202305 <- readRDS("metamouse-20230510.RDS")

##Get list of all reactions contained in the models
lst_allRxns <- c()
for(obj_model in obj_metabolicModelsMouse202305){
  lst_allRxns <- union(lst_allRxns, obj_model@react_id)
}

##Build incidence matrix of reactions in each model
mtx_rxnInModels <- matrix(0,length(lst_allRxns),length(obj_metabolicModelsMouse202305))
rownames(mtx_rxnInModels) <- lst_allRxns
colnames(mtx_rxnInModels) <- names(obj_metabolicModelsMouse202305)
for(str_modelName in names(obj_metabolicModelsMouse202305)) {
  obj_model <- obj_metabolicModelsMouse202305[[str_modelName]]
  mtx_rxnInModels[obj_model@react_id,str_modelName] <- 1
}

##Normalize matrix such that the sum for each species is "1" -> species with larger genomes have 
##smaller contribution of individual reactions
colSums(mtx_rxnInModels)
mtx_rxnInModels <- prop.table(mtx_rxnInModels,2)
colSums(mtx_rxnInModels)

# Relabel
mtx_metabolicModelsMouse202305 <- mtx_rxnInModels

# # ##Or alternatively load incidence matrix of only ACTIVE reactions (more resembling transcriptomics)- predicted from FVA by Stefano
# mtx_rxnInModels <- readRDS(file = "../combineDataLayers/activeReactions_Stefano/fva_99_reactions_thr06.RDS")

# clean up
rm( mtx_rxnInModels, obj_metabolicModelsMouse202305, lst_allRxns, str_modelName, obj_model)


#####################################
## Additional metadata (like taxonomy)
df_mouseMMmetaData <- read.csv(file = "officialMAGmetaData_2023-11-02.csv", 
                               header = T, sep = "\t", stringsAsFactors = F, check.names = F, nrows = 181)
# Fix Phylum names
df_mouseMMmetaData$Phylum <- sub(pattern = "_[A-K]", replacement = "", perl = T, x = df_mouseMMmetaData$Phylum)
# Rownames
rownames(df_mouseMMmetaData) <- df_mouseMMmetaData$MAG_originalID


#####################################
## PCoA
dist_metabolicModels <- vegdist(x = t(mtx_metabolicModelsMouse202305), 
                                method = "horn", # euclidean, manhattan, horn (for abundance data)
                                na.rm = F)
obj_pca <- prcomp(x = dist_metabolicModels)
# Calculate explained variance for each principle component from standard-deviation  sdev_i²/sum(sdev²)
dbl_varianceExplained <- obj_pca$sdev^2 / sum( obj_pca$sdev^2 )
# plot(obj_pca)


################################
### associate variables to coordinates:
df_pcaMouseOnly <- obj_pca$x

colnames(df_mouseMMmetaData)
# more selective metadata
df_tmp <- df_mouseMMmetaData
df_tmp <- df_tmp[rownames(df_pcaMouseOnly),
                 c("Order", "Phylum", "16S rRNA", "# tRNAs", "Completeness [%]",
                   "Contamination [%]", "GC-content [%]", "Size [Mbp]", 
                   "# Reactions", "# Gaps")]
df_tmp$`# Reactions` <- as.numeric(sub(",","",df_tmp$`# Reactions`))
# df_tmp$`# Gaps` <- as.numeric(sub(",","",df_tmp$X..Gaps))

colnames(df_tmp) <- c("Order", "Phylum", "X16S.rRNA", "tRNAs", "Completeness", "Contamination",
                      "GC-content", "Genome Size", "Model Size", "Model Gaps")

length(unique(df_tmp$Order))
# 18

all.equal(rownames(df_tmp), rownames(df_pcaMouseOnly))
# TRUE

obj_metaFitWithOrder <- vegan::envfit(df_pcaMouseOnly, df_tmp[,c(1,3:ncol(df_tmp))])
obj_metaFitWithOrder
# weirdly enough the sum of r2 is greater 1 - 
#  which is because each varialbe is fitted seperately
obj_metaFitWithOrder$vectors$r
obj_metaFitWithOrder$factors$r

obj_metaFit <- vegan::envfit(df_pcaMouseOnly, df_tmp[,3:ncol(df_tmp)])
obj_metaFit$vectors

df_phylumCols <- readRDS("df_phylumCols.rds")

pdf(file = "pca_mouseOnlyReactionsWithMetaDataMapping.pdf", width = 300/72, height = 300/72, pointsize = 8)
par(mar = c(4,4,0,0)+.1)
# plot it
plot(df_pcaMouseOnly,
     col = df_phylumCols[df_tmp$Phylum,"Color"],
     pch = as.numeric(factor(df_tmp$Order)),
     xlab = paste0("PC1 (",round(dbl_varianceExplained[1]*100,1)," %)"),
     ylab = paste0("PC2 (",round(dbl_varianceExplained[2]*100,1)," %)"),
     xlim = c(-2.7,2.5), ylim = c(-1,1.1),
     lwd = 1.4, cex = 1.4, cex.axis = 1.4, cex.lab = 1.4)
# 16S
arrows(x0 = 0, y0 = 0,
       x1 = obj_metaFit$factors$centroids["X16S.rRNAyes",1],
       y1 = obj_metaFit$factors$centroids["X16S.rRNAyes",2],
       length = .1,
       col = 1, lwd = 1.8)
text(x = obj_metaFit$factors$centroids["X16S.rRNAyes",1]*1.8,
     y = obj_metaFit$factors$centroids["X16S.rRNAyes",2]*1.4, labels = "16S rRNA",
     cex = 1.4)
lapply(X = 1:nrow(obj_metaFit$vectors$arrows),FUN = function(int_i) {
        str_label <- rownames(obj_metaFit$vectors$arrows)[int_i]
        obj_row <- obj_metaFit$vectors$arrows[int_i,]
        
        arrows(x0 = 0, y0 = 0, x1 = obj_row[1], y1 = obj_row[2],
               length = .1, col = 1, lwd = 1.8)
        text(x = obj_row[1]*1.8, y = obj_row[2]*1.1, labels = str_label,
             cex = 1.4)
      })
# legend("bottomleft",
#        legend = unique(df_tmp$Order), title = "Microbiota Order:",
#        pch = unique(as.numeric(factor(df_tmp$Order))),
#        pt.cex = 1.2, pt.lwd = 1.5, bty = "n", cex = .9, y.intersp = .5,
#        col = df_phylumCols[df_tmp[!duplicated(df_tmp$Order),"Phylum"],"Color"],
#        title.font = 2, title.adj = c(0,0), title.cex = .95
# )
dev.off()
# export 900x800 pca_mouseOnlyReactionsWithMetaDataMapping

pdf("pcaLegend.pdf", width = 230/72, height = 230/72, pointsize = 8)
plot(1,1, pch = NA, bty = "n", yaxt = "n", xaxt = "n", ylab = "", xlab = "")
legend("center",
       legend = unique(df_tmp$Order), title = "Microbiota Order:",
       pch = unique(as.numeric(factor(df_tmp$Order))),
       pt.cex = 1.6, pt.lwd = 2, bty = "n", cex = 1.4, y.intersp = .75,
       col = df_phylumCols[df_tmp[!duplicated(df_tmp$Order),"Phylum"],"Color"],
       title.font = 2, title.adj = c(0,0), title.cex = 1.4, xpd = T
)
dev.off()

