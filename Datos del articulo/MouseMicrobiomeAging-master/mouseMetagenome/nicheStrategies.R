##############################################################
##
##    Niche colonizing strategies by Johannes (predicted from models)
##
##############################################################

###
source("customFunctions.R")
library(RColorBrewer)
# library(clusterProfiler)
library(coin)

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


########################
# load predicted strategies
df_modelNicheStrategies <- read.csv(file = "S4.14_AdaptiveStrategiesMAGs.csv", 
                                    header = T, sep = ",", row.names = 1, stringsAsFactors = F)

# depending on how many niches a bug takes, it gets weighted down for multi-use
df_modelNicheStrategies$weight <- 1/(stringi::stri_count(str = df_modelNicheStrategies$csr, fixed = ",")+1)

# multiply niche strategies with model/MAG abundances and weight, then get ColSums  
df_NicheStrategieAbundances <- data.frame(ruderal = colSums( df_MAGabundances[rownames(df_modelNicheStrategies)[grep(pattern = "ruderal", x = df_modelNicheStrategies$csr, ignore.case = T, value = F)],]  * 
                                                              df_modelNicheStrategies[grep(pattern = "ruderal", x = df_modelNicheStrategies$csr, ignore.case = T, value = F),"weight"] ),
                                          competition = colSums( df_MAGabundances[rownames(df_modelNicheStrategies)[grep(pattern = "competition", x = df_modelNicheStrategies$csr, ignore.case = T, value = F)],] * 
                                                                  df_modelNicheStrategies[grep(pattern = "competition", x = df_modelNicheStrategies$csr, ignore.case = T, value = F),"weight"] ),
                                          stress.toleration = colSums( df_MAGabundances[rownames(df_modelNicheStrategies)[grep(pattern = "stress.toleration", x = df_modelNicheStrategies$csr, ignore.case = T, value = F)],] * 
                                                                        df_modelNicheStrategies[grep(pattern = "stress.toleration", x = df_modelNicheStrategies$csr, ignore.case = T, value = F),"weight"] )
                                          )
df_NicheStrategieAbundances <- t( df_NicheStrategieAbundances )
colSums(df_NicheStrategieAbundances)
df_NicheStrategieAbundances <- as.data.frame( t(t(df_NicheStrategieAbundances)/colSums(df_NicheStrategieAbundances))*100 )
colSums(df_NicheStrategieAbundances)
df_NicheStrategieAbundances <- df_NicheStrategieAbundances[,df_metadataMAG$SequencingID]
colnames(df_NicheStrategieAbundances) <- df_metadataMAG$SampleName


#######################
# quick box plots
pdf(file = "boxplot_nicheStrategiesByAge_gapseq.pdf")
for(int_i in 1:3) {
  print(rownames(df_NicheStrategieAbundances)[int_i])
  boxplot(as.numeric( df_NicheStrategieAbundances[int_i,df_metadataMAG$SampleName] ) ~ df_metadataMAG$Age,
          col = brewer.pal("Pastel1", n = 5),
          main = Hmisc::capitalize(rownames(df_NicheStrategieAbundances)[int_i]),
          xlab = "Age [months]", ylab = "rel. abundance [%]")
  print( 
    DescTools::DunnTest(x = as.numeric( df_NicheStrategieAbundances[int_i,df_metadataMAG$SampleName] ), 
                        g = df_metadataMAG$Age, 
                        method = "BH",
                        alternative = "two.sided")
  )
}
dev.off()
# [1] "ruderal"
# Dunn's test of multiple comparisons using rank sums : BH  
# mean.rank.diff   pval    
# 9-2        -3.100000 0.7193    
# 15-2        0.000000 1.0000    
# 24-2        9.600000 0.2611    
# 30-2       13.216667 0.1389    
# 15-9        3.100000 0.7193    
# 24-9       12.700000 0.1524    
# 30-9       16.316667 0.1192    
# 24-15       9.600000 0.2611    
# 30-15      13.216667 0.1389    
# 30-24       3.616667 0.7193    
# 
# [1] "competition"
# mean.rank.diff   pval    
# 9-2         5.200000 0.4921    
# 15-2       -1.400000 0.8363    
# 24-2      -10.300000 0.2571    
# 30-2      -19.283333 0.0148 *  
# 15-9       -6.600000 0.4127    
# 24-9      -15.500000 0.0555 .  
# 30-9      -24.483333 0.0016 ** 
# 24-15      -8.900000 0.2702    
# 30-15     -17.883333 0.0195 *  
# 30-24      -8.983333 0.2702    
# 
# [1] "stress.toleration"
# mean.rank.diff   pval    
# 9-2        -7.500000 0.4914    
# 15-2       -2.100000 0.7567    
# 24-2        5.000000 0.5119    
# 30-2       10.766667 0.2427    
# 15-9        5.400000 0.5119    
# 24-9       12.500000 0.2171    
# 30-9       18.266667 0.0488 *  
# 24-15       7.100000 0.4914    
# 30-15      12.866667 0.2171    
# 30-24       5.766667 0.5119    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#######################
# publishable figure
pdf(file = "boxplot_nicheStrategiesByAge_gapseq.pdf", width = 150/72, height = 280/72, pointsize = 10)
par(mar = c(0,4,4,0)+.1,
    mfrow = c(3,1))
int_i <- 1
boxplot(as.numeric( df_NicheStrategieAbundances[int_i,df_metadataMAG$SampleName] ) ~ df_metadataMAG$Age,
        col = brewer.pal("Pastel2", n = 5)[2],
        ylim = c(30,90), xlab = "", xaxt = "n", las = 2,
        cex.axis = 1.3, cex.lab = 1.3,
        frame = F, ylab = "Relative abundance [%]")
points(y = as.numeric( df_NicheStrategieAbundances[int_i,df_metadataMAG$SampleName] ),
       x = jitter( df_metadataMAG$AgeGroup, factor = 2 ), col = brewer.pal("Dark2", n = 5)[2], 
       pch = 4, cex = 1, lwd = 1.1)
lines(x = c(2,5), y = rep(85,2)); text(x = 3.5, y = 89, label = "**", cex = 1.5)
lines(x = c(3,5), y = rep(75,2)); text(x = 4, y = 79, label = "*", cex = 1.5)
# 30 - 9
# 30 - 15
par(mar = c(2,4,2,0)+.1)
int_i <- 2
boxplot(as.numeric( df_NicheStrategieAbundances[int_i,df_metadataMAG$SampleName] ) ~ df_metadataMAG$Age,
        col = brewer.pal("Pastel2", n = 8)[8],
        ylim = c(10,80), xlab = "", xaxt = "n", las = 2,
        cex.axis = 1.3, cex.lab = 1.3,
        frame = F, ylab = "Relative abundance [%]")
points(y = as.numeric( df_NicheStrategieAbundances[int_i,df_metadataMAG$SampleName] ),
       x = jitter( df_metadataMAG$AgeGroup, factor = 2 ), col = brewer.pal("Dark2", n = 8)[8], 
       pch = 4, cex = 1, lwd = 1.1)
lines(x = c(1,1.9), y = rep(75,2)); text(x = 1.45, y = 79, label = "*", cex = 1.5)
lines(x = c(2.1,5), y = rep(75,2)); text(x = 3.55, y = 79, label = "**", cex = 1.5)
lines(x = c(3,5), y = rep(65,2)); text(x = 4, y = 69, label = "*", cex = 1.5)
# 9 -2 *
# 30 -9 **
# 30 - 15 *
par(mar = c(4,4,0,0)+.1)
int_i <- 3
boxplot(as.numeric( df_NicheStrategieAbundances[int_i,df_metadataMAG$SampleName] ) ~ df_metadataMAG$Age,
        col = brewer.pal("Pastel2", n = 8)[5],
        ylim = c(5,45), xlab = "Age [months]", las = 1,
        cex.axis = 1.3, cex.lab = 1.3,
        frame = F, ylab = "Relative abundance [%]")
points(y = as.numeric( df_NicheStrategieAbundances[int_i,df_metadataMAG$SampleName] ),
       x = jitter( df_metadataMAG$AgeGroup, factor = 2 ), col = brewer.pal("Dark2", n = 8)[5], 
       pch = 4, cex = 1, lwd = 1.1)
lines(x = c(1,5), y = rep(42,2)); text(x = 3, y = 43.5, label = "**", cex = 1.5)
lines(x = c(1,4), y = rep(38,2)); text(x = 2.5, y = 39.5, label = "***", cex = 1.5)
lines(x = c(1,3), y = rep(34,2)); text(x = 2, y = 35.5, label = "**", cex = 1.5)
lines(x = c(1,2), y = rep(30,2)); text(x = 1.5, y = 31.5, label = "**", cex = 1.5)
# 9-2       -23.800000 0.00223 ** 
# 15-2      -22.000000 0.00389 **
# 24-2      -26.800000 0.00077 ***
# 30-2      -20.533333 0.00389 **
dev.off()

