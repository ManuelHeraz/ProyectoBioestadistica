################################################################
###
###   Grwoth Rate Prediction via Peak-to-Trough evaluation
###
################################################################

library(RColorBrewer)
library(clusterProfiler)

source("customFunctions.R")
setwd("mouseMetagenome/")


###################################
##
df_growthPredictionCOPTR <- read.table(file = "growthPredictionCOPTR_TD.csv", 
                                       header = T, sep = ",", row.names = 1, stringsAsFactors = F, check.names = F)
df_growthPredictionCOPTR <- df_growthPredictionCOPTR[,df_metadataMAG$SequencingID]
colnames(df_growthPredictionCOPTR) <- df_metadataMAG$SampleName


###################################
## Growth prediction of each MAG by Age
df_lmMAGGrowthByAge <- df_linearModel(df_dependantVariable = df_growthPredictionCOPTR[,df_metadataMAG$SampleName], 
                                      lst_var_metaData = df_metadataMAG$Age, 
                                      int_cores = 7)

par(mar = c(5,4,2,2)+.1)
hist(df_lmMAGGrowthByAge$p.value, breaks = 30)

table(df_lmMAGGrowthByAge$p.adj <= 0.05)
# FALSE  TRUE 
#  156    22 

#####
### publishable figure
## FBA community growth by age and COPTR community growth by age
pdf(file = "boxplots_communityGrowthByAge.pdf", width = 180/72, height = 280/72, pointsize = 7)
par(mar = c(1,5,3,1)+.2,
    mfrow = c(2,1))
boxplot(as.numeric(df_growthRate) ~ df_metadataFBA[colnames(df_growthRate),"Age"],
        col = brewer.pal("Pastel2",n=5)[5], xlab = "", ylab = "",
        xaxt = "n", las = 2, cex.axis = 1.3, cex.lab = 1.3, frame = F, ylim = c(0.008,0.033))
points(y = as.numeric(df_growthRate),
       x = jitter( df_metadataFBA[colnames(df_growthRate),"AgeGroup"], factor = 2 ), 
       col = brewer.pal("Dark2", n = 8)[5], 
       pch = 4, cex = 1, lwd = 1.1)
title(ylab = "Community Growth", line = 4, cex.lab = 1.3 )
lines(x = c(1,5), y = rep(32/1000,2)); text(x = 3, y = 33/1000, label = "*", cex = 1.5)
lines(x = c(1,4), y = rep(29/1000,2)); text(x = 2.5, y = 30/1000, label = "**", cex = 1.5)
lines(x = c(2,4), y = rep(26/1000,2)); text(x = 3, y = 27/1000, label = "*", cex = 1.5)
lines(x = c(3,4), y = rep(23/1000,2)); text(x = 3.5, y = 24/1000, label = "*", cex = 1.5)
#       mean.rank.diff   pval    
# 24-2           -22.7 0.0050 ** 
# 30-2           -18.3 0.0175 *  
# 24-9           -16.7 0.0260 *  
# 24-15          -18.2 0.0175 *  
par(mar = c(4,5,0,1)+.2)
###
# Median 
dbl_medianCOPTRgrowth <- apply(X = df_growthPredictionCOPTR[,df_metadataMAG$SampleName], MARGIN = 2, FUN = median, na.rm = T)
boxplot(dbl_medianCOPTRgrowth ~ df_metadataMAG$Age, col = brewer.pal("Pastel2",n=5)[3], 
        xlab = "Age [months]", ylab = "", ylim = c(.46,.82),
        xaxt = "n", las = 2, cex.axis = 1.3, cex.lab = 1.3, frame = F)
points(y = dbl_medianCOPTRgrowth,
       x = jitter( df_metadataMAG$AgeGroup, factor = 2 ), 
       col = brewer.pal("Dark2", n = 8)[3], 
       pch = 4, cex = 1, lwd = 1.1)
axis(side = 1, at = 1:5, las = 1, cex.axis = 1.3, labels = c(2,9,15,24,30))
title(ylab = "Community Growth", line = 4, cex.lab = 1.3 )
lines(x = c(1,5), y = rep(81/100,2)); text(x = 3, y = 82/100, label = "***", cex = 1.5)
lines(x = c(1,4), y = rep(78/100,2)); text(x = 2.5, y = 79/100, label = "**", cex = 1.5)
lines(x = c(1,3), y = rep(75/100,2)); text(x = 2, y = 76/100, label = "**", cex = 1.5)
lines(x = c(1,2), y = rep(72/100,2)); text(x = 1.5, y = 73/100, label = "***", cex = 1.5)
#       mean.rank.diff    pval    
# 9-2       -26.200000 0.00091 ***
# 15-2      -21.900000 0.00411 ** 
# 24-2      -20.000000 0.00792 ** 
# 30-2      -24.283333 0.00091 ***
dev.off()
#350x500


###################################
## FBA community growth for each sample (mouse) vs. COPTR derived "community growth" (median)

# build median of MAG-growth over all MAGs and correlate to FBA Community Growth rates
cor.test(x = apply(X = df_growthPredictionCOPTR, 
                   MARGIN = 2, 
                   FUN = median, na.rm = T)[colnames(df_growthRate)],
         y = as.numeric(df_growthRate[1,]), method = "spearman")
# p-value = 0.001237   rho 0.4478752 

# Let's have a small figure of it
par(mar = c(5,2,2,2)+.1)
plot(x = apply(X = df_growthPredictionCOPTR, 
               MARGIN = 2, 
               FUN = median, na.rm = T)[colnames(df_growthRate)],
     y = as.numeric(df_growthRate[1,]),
     xlab = "Community Growth Rate COPTR", ylab = "Community Growth Rate FBA", 
     main = "Median", pch = 4, col = rgb(.3,.3,.4,.8))
abline(coef = lm(as.numeric(df_growthRate[1,]) ~
                   apply(X = df_growthPredictionCOPTR, MARGIN = 2, 
                         FUN = median, na.rm = T)[colnames(df_growthRate)])$coef,
       col = 2, lty = 2)

# publishable figure:
pdf("scatter_communityGrowthMAGsToFBA.pdf", width = 170/72, height = 150/72, pointsize = 8)
par(mar = c(4,5,2,1)+.1)
plot(x = NA,
     y = NA,
     ylim = c(0.008,0.028), xlim = c(0.45,.72), yaxt = "n",
     cex.lab = 1.2, cex.axis = 1.1, pch = 21, cex = 1.4,
     xlab = "COPTR Growth", ylab = "",
     bg = brewer.pal("Pastel1", n = 5)[df_metaDataMM[colnames(df_growthRate),"AgeGroup"]],
     bty = "n", xaxt = "n"
     )
axis(side = 1, at = c(.4,.5,.6,.7,.8), cex.axis = 1.1)
axis(side = 2, at = seq(from = .5, to = 3, by = 0.5)/100, labels = c(NA,.01,NA,.02,NA,NA), las = 3, cex.axis = 1.1)
title(ylab = "FBA Growth", line = 3, cex.lab = 1.2) #font.lab = 2,
abline(coef = lm(as.numeric(df_growthRate[1,]) ~
                   apply(X = df_growthPredictionCOPTR, MARGIN = 2, 
                         FUN = median, na.rm = T)[colnames(df_growthRate)])$coef,
       col = rgb(.5,.5,.6), lty = 5, lwd = 1.6)
points(x = apply(X = df_growthPredictionCOPTR, 
                 MARGIN = 2, 
                 FUN = median, na.rm = T)[colnames(df_growthRate)],
       y = as.numeric(df_growthRate[1,]),
       pch = 21, cex = 1.4,
       bg = brewer.pal("Pastel1", n = 5)[df_metaDataMM[colnames(df_growthRate),"AgeGroup"]])

# legend("topleft", fill = brewer.pal("Pastel1", n = 5), legend = c(2,9,15,24,30),# bty = "n",
#        title = "Age [months]", cex = 1.1, title.cex = 1.2, title.font = 1, ncol = 2)
dev.off()

