##############################################################
##
##    Metagenomes of ageing mice baseline aging
##
##############################################################

###
library(RColorBrewer)
library(clusterProfiler)
library(coin)
library(vegan)

setwd("mouseMetagenome/")


##############
# load metadata
df_metadataMAG <- readRDS("df_metadataMAG.rds")
rownames(df_metadataMAG) <- df_metadataMAG$SampleName


##############
# load abundance data - created with Metabat jgi_summarize - normalized for Contig lengths
df_MAGabundances <- readRDS("df_MAGabundances.rds")
colSums(df_MAGabundances)
# normalize for sample depth and multiply by a million to generate TPM-style data
df_MAGabundancesTPM <- as.data.frame( t(t(df_MAGabundances)/colSums(df_MAGabundances))*1e6 )
colSums(df_MAGabundancesTPM)
df_MAGabundancesTPM <- df_MAGabundancesTPM[,df_metadataMAG$SequencingID]
colnames(df_MAGabundancesTPM) <- df_metadataMAG$SampleName


#######################
## Evaluate data
void_MDSplot(df_sampleMatrix = df_MAGabundances, df_metaData = df_metadataMAG, bol_labels = T, 
             var_colorBy = "Age", var_shapeBy = NULL, bol_continousColor = F)
# TPM
void_MDSplot(df_sampleMatrix = df_MAGabundancesTPM, df_metaData = df_metadataMAG, bol_labels = F, 
             var_colorBy = "Age", var_shapeBy = NULL, bol_continousColor = F)


############################
## Plot General MAG abundance in cohort
df_tmp <- df_MAGabundances[order(rowSums(df_MAGabundances), decreasing = T),]
lst_colors <- df_phylumCols[sub(pattern = "_[AB]", replacement = "", x = df_gtdbtkTaxonmy2022[rownames(df_tmp),"Phylum"], perl = T),"Color"]
rownames(df_tmp) <- df_updateMAGnames[rownames(df_tmp),"newID"]
par(mar = c(11,4,1,0)+.1)
barplot(rowSums(df_tmp),
        main = "MAG coverage in full cohort", ylab = "Total Read Depth", xlab = "",
        col = lst_colors, las = 2, cex.axis = .9, cex.names = .6, border = F, space = 0, horiz = F)
# boxplot(scales::rescale(rowSums(df_tmp), to = c(1,181)), horizontal = T, at = 2000, add= T )
# abline(h = quantile(rowSums(df_tmp), probs = c(.25,.50,.75)), col = 8, lty = 2)
abline(v = min(which(rowSums(df_tmp) <= quantile(rowSums(df_tmp), probs = c(.25)))), col = 8, lty = 2)
abline(v = min(which(rowSums(df_tmp) <= quantile(rowSums(df_tmp), probs = c(.5)))), col = 8, lty = 2)
# 91
abline(v = min(which(rowSums(df_tmp) <= quantile(rowSums(df_tmp), probs = c(.75)))), col = 8, lty = 2)
# export at 1800x500
## only top 91:
barplot(rowSums(df_tmp[1:91,]),
        main = "Most abundant MAG coverage in full cohort", ylab = "Total Read Depth", xlab = "",
        col = lst_colors, las = 2, cex.axis = .9, cex.names = .6, border = F, space = 0, horiz = F)

####
rm(df_tmp, lst_colors)
# # perfect correlation as it should be- so it does not matter if we look at TPM or Absolut numbers
# plot(rowSums(df_MAGabundances), rowSums(df_MAGabundancesTPM))


#########################################
#
df_lmMAGabundances <- df_linearModel(df_dependantVariable = df_MAGabundancesTPM[,df_metadataMAG$SampleName],
                                     lst_var_metaData = df_metadataMAG$Age,
                                     int_cores = 7)
# [1] "Starting rank sum tests of 181 observations within 52 samples! Tue May 16 10:22:47 2023"
# [1] "Finished all tests! Preparing results. Tue May 16 10:22:47 2023"

## few diagnostic plots
# p-val Histogram
hist(df_lmMAGabundances$p.value, breaks = 40)

# PCoA
void_MDSplot(df_sampleMatrix = df_MAGabundances[rownames(df_lmMAGabundances)[df_lmMAGabundances$p.adj <= 0.05],],
             df_metaData = df_metadataMAG, bol_labels = F, 
             var_colorBy = "Age", var_shapeBy = NULL, bol_continousColor = F)

# Sort by p-vals
df_lmMAGabundances <- as.data.frame(df_lmMAGabundances[order(df_lmMAGabundances$p.value, decreasing = F),])
df_lmMAGabundances$MAG <- rownames(df_lmMAGabundances)
rownames(df_lmMAGabundances) <- NULL

######## 
## alternative ploting method
mtx_tmp <- as.matrix(df_MAGabundancesTPM[rownames(df_MAGabundancesTPM) %in% df_lmMAGabundances[df_lmMAGabundances$p.adj <= 0.05, "MAG"],
                                      df_metadataMAG[order(df_metadataMAG$Age, decreasing = F),"SampleName"]])
mtx_tmp <- mtx_tmp[order(rowMeans(mtx_tmp[,df_metadataMAG[df_metadataMAG$Age == 30,"SampleName"]]) / 
                           rowMeans(mtx_tmp[,df_metadataMAG[df_metadataMAG$Age == 9,"SampleName"]]),decreasing = T),]
# row scaling and mean per age-group
df_tmp <- as.matrix( as.data.frame( apply(X = mtx_tmp[,df_metadataMAG$SampleName], MARGIN = 1, FUN = function(row) {
                                            return(
                                                    scale( tapply(row, df_metadataMAG$Age, mean), center = F, scale = T)
                                                  )
                                          }) ))
# update MAG names to publishable IDs
colnames(df_tmp) <- df_updateMAGnames[colnames(df_tmp),"newID"]

pdf(file = "circles_signDifMAGsMeanAbundanceByAge.pdf", width = 200/72, height = 430/72, pointsize = 8)
par(mar = c(1,14,2,1)+.1)
plot(x = rep(1:5, ncol(df_tmp)), 
     y = rep(1:ncol(df_tmp),each = 5), 
     cex = scales::rescale(as.numeric(df_tmp),
                           to = c(0.4,2.4)), 
     pch = 21, yaxt = "n", ylab = "", xaxt = "n", 
     xlim = c(0.5,5.5), xlab = "", cex.axis = 1.1,
     bg = rep( df_phylumCols[ sub(pattern = "_[A-J]", replacement = "", x = df_fullMAGmetaData[colnames(df_tmp), "Phylum"]),
                             "Color"], each = 5),
     col = 1, frame.plot = F,
     )
axis(side = 2, tick = F, 
     line = 12, hadj = 0,
     at = 1:ncol(df_tmp),
     labels = colnames(df_tmp),
     cex.axis = 1, las = 2)
# axis(side = 3, at = 1:5, tick = F, line = -2, font = 2,
#      labels = c(2,9,15,24,30), cex.axis = 1.1)
dev.off()
# export at 420x780px as circles_signDifMAGsMeanAbundanceByAge


###########################
## get alpha diversity
# rarefaction
rarecurve(t(round(df_MAGabundances*1000,0)), label = F)
min( colSums(df_MAGabundances*1000) )
# [1] 435755.1
df_MAGabundancesRare <- t( rrarefy(t(round(df_MAGabundances*1000,0)), 435000))
df_MAGabundancesRare <- df_MAGabundancesRare[order(rowSums(df_MAGabundancesRare), decreasing = T),]

dbl_alphaDiversityShannonMAGs <- vegan::diversity(x = df_MAGabundancesRare, index="shannon", MARGIN = 2)
mean(dbl_alphaDiversityShannonMAGs)
# [1] 3.730097
median(dbl_alphaDiversityShannonMAGs)
# [1] 3.774068

boxplot(dbl_alphaDiversityShannonMAGs[df_metadataMAG$SequencingID] ~ df_metadataMAG$Age,
        col = brewer.pal("Pastel1", n = 5), xlab= "Age [months]", ylab = "Shannon Index",
        main = "Alpha Diversity"
        )
DescTools::DunnTest(x = dbl_alphaDiversityShannonMAGs[df_metadataMAG$SequencingID],
                    g = df_metadataMAG$Age, method = "BH")
# mean.rank.diff   pval    
# 9-2             16.1 0.1752    
# 15-2            10.8 0.2776    
# 24-2             2.7 0.7462    
# 30-2             4.8 0.5743    
# 15-9            -5.3 0.5743    
# 24-9           -13.4 0.2401    
# 30-9           -11.3 0.2720    
# 24-15           -8.1 0.4641    
# 30-15           -6.0 0.5743    
# 30-24            2.1 0.7462 

