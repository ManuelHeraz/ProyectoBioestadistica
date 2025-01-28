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


###############
# organize into higher level groups
df_MBpathwayGroups <- read.table(file = "../databases/MBpathwayGroups.tsv", header = T, sep = "\t", 
                                 stringsAsFactors = F, comment.char = "", as.is = T, quote = "\"")
df_MBpathwayGroups[,3] <- NULL
colnames(df_MBpathwayGroups) <- c("MB_TermID", "MB_SubSystem", "MB_Group")
str_MBGroupTerms <- sort(unique(df_MBpathwayGroups$MB_Group))

# df_MBpathwayGroups <- rbind(df_MBpathwayGroups,
#                             c("fatty acid biosynthesis type II","long chain fatty acids"))
df_hostPathwayGroups <- read.table(file = "../databases/HostPathwayGroups.tsv", header = T, sep = "\t", stringsAsFactors = F)
df_hostPathwayGroups <- df_hostPathwayGroups[,c(1,4)]
colnames(df_hostPathwayGroups) <- c("HOST_TermID","HOST_Group")
str_hostGroupTerms <- sort(unique(unlist(strsplit(split = ",", fixed = T, x = df_hostPathwayGroups$HOST_Group))))

#
df_SubSysdescriptions <- unique( df_rxn2subsys[,c("SubsysID","Subsystem")] )
rownames(df_SubSysdescriptions) <- df_SubSysdescriptions$SubsysID


############################
# summary level table for colon
df_colonAssocScoreSummary <- t(sapply(X = str_hostGroupTerms, FUN = function(str_hostGroup) {
                                    df_tmp <- df_colonAssocScoresOverview[grep(pattern = str_hostGroup, x = df_colonAssocScoresOverview$HOST_Group, fixed = T),]
                                    
                                    int_scores <- sapply(X = str_MBGroupTerms,FUN = function(str_MBgroup) {
                                                          df_tmp2 <- df_tmp[grep(pattern = str_MBgroup, x = df_tmp$MB_Group, fixed = T), c(4,7:11)]
                                                          df_tmp2 <- df_tmp2[!duplicated(df_tmp2),]
                                                          
                                                          return( sum(abs(df_tmp2$SignedLogPadj)) )
                                                        })
                                    return( int_scores )
                                  }))
write.table(x = df_colonAssocScoreSummary,
            file = "colonAssocScoreGroups.tsv",
            append = F, quote = T, sep = "\t", row.names = T, col.names = T)
df_colonAssocScoreSummary <- df_colonAssocScoreSummary[order(rowSums(df_colonAssocScoreSummary), decreasing = T),
                                                       order(colSums(df_colonAssocScoreSummary), decreasing = T)]


############################
# let's build a figure
rowSums(df_colonAssocScoreSummary) 
colSums(df_colonAssocScoreSummary)
mtx_tmp <- t( df_colonAssocScoreSummary[1:5,1:5] )
# shorten labels
rownames(mtx_tmp) <- Hmisc::capitalize( str_simplifyDescriptions(rownames(mtx_tmp)) )
rownames(mtx_tmp)[5] <- "Nucleotides"
colnames(mtx_tmp) <- Hmisc::capitalize( str_simplifyDescriptions(colnames(mtx_tmp)) )
colnames(mtx_tmp)[4] <- "Localization"

pdf(file = "contingencyPlot_colonAssocScoreOverview.pdf", width = 260/72, height = 220/72, pointsize = 7)
par(mar = c(10,14,10,10)+.1)
myDotPlot(mtx_input = mtx_tmp[,ncol(mtx_tmp):1], scaleTo = c(1,3.2), labelCex = 1, color = "#01ba3877")
# add labels
text(x = -1.5, y = ncol(mtx_tmp)+.5,
     xpd = T, cex = 1.2, font = 2, adj = c(.5,.5),
     labels = "Host")
text(x = -0.5, y = -1, xpd = T, cex = 1.2, font = 2,
     srt = 30, adj = c(1,0),
     labels = "Microbiome")
dev.off()
# 650x560 contingencyPlot_colonAssocScoreOverview


################################
# summary level table for brain
df_brainAssocScoreSummary <- t(sapply(X = str_hostGroupTerms, FUN = function(str_hostGroup) {
  df_tmp <- df_brainAssocScoresOverview[grep(pattern = str_hostGroup, x = df_brainAssocScoresOverview$HOST_Group, fixed = T),]
  
  int_scores <- sapply(X = str_MBGroupTerms,FUN = function(str_MBgroup) {
    df_tmp2 <- df_tmp[grep(pattern = str_MBgroup, x = df_tmp$MB_Group, fixed = T), c(4,7:11)]
    df_tmp2 <- df_tmp2[!duplicated(df_tmp2),]
    
    return( sum(abs(df_tmp2$SignedLogPadj)) )
  })
  return( int_scores )
}))
write.table(x = df_brainAssocScoreSummary,
            file = "brainAssocScoreGroups.tsv",
            append = F, quote = T, sep = "\t", row.names = T, col.names = T)
df_brainAssocScoreSummary <- df_brainAssocScoreSummary[order(rowSums(df_brainAssocScoreSummary), decreasing = T),
                                                       order(colSums(df_brainAssocScoreSummary), decreasing = T)]

############################
# let's build a figure
rowSums(df_brainAssocScoreSummary) 
colSums(df_brainAssocScoreSummary)
mtx_tmp <- t( df_brainAssocScoreSummary[1:5,1:5] )
# shorten labels
rownames(mtx_tmp) <- Hmisc::capitalize(str_simplifyDescriptions(rownames(mtx_tmp)))
colnames(mtx_tmp) <- Hmisc::capitalize(str_simplifyDescriptions(colnames(mtx_tmp)))
#
colnames(mtx_tmp)[1] <- "Localization"
rownames(mtx_tmp)[1] <- "Nucleotides"

pdf(file = "contingencyPlot_brainAssocScoreOverview.pdf", width = 260/72, height = 220/72, pointsize = 7)
par(mar = c(10,14,10,10)+.1)
myDotPlot(mtx_input = mtx_tmp[,ncol(mtx_tmp):1], scaleTo = c(1,3.2), labelCex = 1, color = "#f8766d77")
text(x = -1.5, y = ncol(mtx_tmp)+.5,
     xpd = T, cex = 1.2, font = 2, adj = c(.5,.5),
     labels = "Host")
text(x = -0.5, y = -1, xpd = T, cex = 1.2, font = 2,
     srt = 30, adj = c(1,0),
     labels = "Microbiome")
dev.off()
# 650x560 contingencyPlot_brainAssocScoreOverview


################################
# summary level table for liver
df_liverAssocScoreSummary <- t(sapply(X = str_hostGroupTerms, FUN = function(str_hostGroup) {
  df_tmp <- df_liverAssocScoresOverview[grep(pattern = str_hostGroup, x = df_liverAssocScoresOverview$HOST_Group, fixed = T),]
  
  int_scores <- sapply(X = str_MBGroupTerms,FUN = function(str_MBgroup) {
    df_tmp2 <- df_tmp[grep(pattern = str_MBgroup, x = df_tmp$MB_Group, fixed = T), c(4,7:11)]
    df_tmp2 <- df_tmp2[!duplicated(df_tmp2),]
    
    return( sum(abs(df_tmp2$SignedLogPadj)) )
  })
  return( int_scores )
}))
write.table(x = df_liverAssocScoreSummary,
            file = "liverAssocScoreGroups.tsv",
            append = F, quote = T, sep = "\t", row.names = T, col.names = T)
df_liverAssocScoreSummary <- df_liverAssocScoreSummary[order(rowSums(df_liverAssocScoreSummary), decreasing = T),
                                                       order(colSums(df_liverAssocScoreSummary), decreasing = T)]

############################
# let's build a figure
rowSums(df_liverAssocScoreSummary) 
colSums(df_liverAssocScoreSummary)
mtx_tmp <- t( df_liverAssocScoreSummary[1:5,1:5] )
# shorten labels
rownames(mtx_tmp) <-  Hmisc::capitalize( str_simplifyDescriptions(rownames(mtx_tmp)) )
colnames(mtx_tmp) <-  Hmisc::capitalize( str_simplifyDescriptions(colnames(mtx_tmp)) )
#
colnames(mtx_tmp)[4] <- "Multicell. proc."
rownames(mtx_tmp)[2] <- "Nucleotides"

pdf(file = "contingencyPlot_liverAssocScoreOverview.pdf", width = 260/72, height = 220/72, pointsize = 7)
par(mar = c(10,14,10,10)+.1)
myDotPlot(mtx_input = mtx_tmp[,ncol(mtx_tmp):1], scaleTo = c(1,3.2), labelCex = 1, color = "#619cff77")
text(x = -1.5, y = ncol(mtx_tmp)+.5,
     xpd = T, cex = 1.2, font = 2, adj = c(.5,.5),
     labels = "Host")
text(x = -0.5, y = -1, xpd = T, cex = 1.2, font = 2,
     srt = 30, adj = c(1,0),
     labels = "Microbiome")
dev.off()
# 650x540 contingencyPlot_liverAssocScoreOverview


##################################
## let's build a combined figure
pdf(file = "contingencyPlots.pdf", width = 800/72, height = 1200/72, pointsize = 18)
par(mfrow = c(3,1))
#####
# Colon
mtx_tmp <- df_colonAssocScoreSummary #[,1:5]
mtx_tmp <- t(mtx_tmp)
# shorten labels
rownames(mtx_tmp) <- NULL
colnames(mtx_tmp) <-  Hmisc::capitalize( str_simplifyDescriptions(colnames(mtx_tmp)) )

par(mar = c(6,14,3,3)+.1)
myDotPlot(mtx_input = mtx_tmp[,ncol(mtx_tmp):1], scaleTo = c(.8,4), labelCex = 1, color = "#b0c4de")
# host sums
barplot(height = scales::rescale(x = c(colSums(mtx_tmp)[ncol(mtx_tmp):1]), 
                                 to = c(0,2)),
        horiz = T, las = 2, axes = F, 
        # col = brewer.pal("Accent", n = ncol(mtx_tmp))[4:1], 
        col = "#b0c4de",
        add = T, offset = nrow(mtx_tmp)+.1, width = .75, space = .325, #.275
        names.arg = "")
# mb sums
barplot(height = scales::rescale(x = c(rowSums(mtx_tmp)), 
                                 to = c(0,2)),
        horiz = F, las = 2, axes = F, 
        # col = brewer.pal("Accent", n = ncol(mtx_tmp))[4:1], 
        col = "#b0c4de",
        add = T, offset = ncol(mtx_tmp)+.2, width = .75, space = .325, 
        names.arg = "")
text(x = -2.5, y = ncol(mtx_tmp)+1,
     xpd = T, cex = 1.2, font = 2, adj = c(.5,.5),
     labels = "Colon")
# 770x550

#####
mtx_tmp <- df_liverAssocScoreSummary[,colnames(df_colonAssocScoreSummary)]
mtx_tmp <- t(mtx_tmp)
# shorten labels
rownames(mtx_tmp) <- NULL
colnames(mtx_tmp) <-  Hmisc::capitalize( str_simplifyDescriptions(colnames(mtx_tmp)) )

par(mar = c(6,14,3,3)+.1)
myDotPlot(mtx_input = mtx_tmp[,ncol(mtx_tmp):1], scaleTo = c(.8,4), labelCex = 1, color = "#b0eede")
# host sums
barplot(height = scales::rescale(x = c(colSums(mtx_tmp)[ncol(mtx_tmp):1]), 
                                 to = c(0,2)),
        horiz = T, las = 2, axes = F, 
        # col = brewer.pal("Accent", n = ncol(mtx_tmp))[4:1], 
        col = "#b0eede",
        add = T, offset = nrow(mtx_tmp)+.1, width = .75, space = .325, 
        names.arg = "")
# mb sums
barplot(height = scales::rescale(x = c(rowSums(mtx_tmp)), 
                                 to = c(0,2)),
        horiz = F, las = 2, axes = F, 
        # col = brewer.pal("Accent", n = ncol(mtx_tmp))[4:1], 
        col = "#b0eede",
        add = T, offset = ncol(mtx_tmp)+.1, width = .75, space = .325, 
        names.arg = "")
text(x = -2.5, y = ncol(mtx_tmp)+1,
     xpd = T, cex = 1.2, font = 2, adj = c(.5,.5),
     labels = "Liver")
# 770x720 contingencyPlot_liverAssocScoreOverview

######
mtx_tmp <- df_brainAssocScoreSummary[,colnames(df_colonAssocScoreSummary)]
mtx_tmp <- t(mtx_tmp)
# shorten labels
rownames(mtx_tmp) <-  Hmisc::capitalize( str_simplifyDescriptions(rownames(mtx_tmp)) )
rownames(mtx_tmp)[5] <- "Nucleotides"
rownames(mtx_tmp)[9] <- "Signaling"
rownames(mtx_tmp)[18] <- "Biosynthesis"
rownames(mtx_tmp)[19] <- "Detoxification"
rownames(mtx_tmp)[22] <- "tRNA"
colnames(mtx_tmp) <-  Hmisc::capitalize( str_simplifyDescriptions(colnames(mtx_tmp)) )

par(mar = c(6,14,3,3)+.1)
myDotPlot(mtx_input = mtx_tmp[,ncol(mtx_tmp):1], scaleTo = c(.8,4), labelCex = 1, color = "#ffc4de")
# host sums
barplot(height = scales::rescale(x = c(colSums(mtx_tmp)[ncol(mtx_tmp):1]), 
                                 to = c(0,2)),
        horiz = T, las = 2, axes = F, 
        # col = brewer.pal("Accent", n = ncol(mtx_tmp))[4:1], 
        col = "#ffc4de",
        add = T, offset = nrow(mtx_tmp)+.1, width = .75, space = .325, 
        names.arg = "")
# mb sums
barplot(height = scales::rescale(x = c(rowSums(mtx_tmp)), 
                                 to = c(0,2)),
        horiz = F, las = 2, axes = F, 
        # col = brewer.pal("Accent", n = ncol(mtx_tmp))[4:1], 
        col = "#ffc4de",
        add = T, offset = ncol(mtx_tmp)+.1, width = .75, space = .325, 
        names.arg = "")
text(x = -2.5, y = ncol(mtx_tmp)+1,
     xpd = T, cex = 1.2, font = 2, adj = c(.5,.5),
     labels = "Brain")
text(x = -1, y = -1, xpd = T, cex = 1.2, font = 2,
     srt = 30, adj = c(1,0),
     labels = "Microbiome")
# 770x550
dev.off()

################################

