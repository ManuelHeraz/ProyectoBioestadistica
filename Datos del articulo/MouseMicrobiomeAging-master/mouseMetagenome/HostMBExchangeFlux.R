############################################################
###
###   Host-MB exchange fluxes from community models FBA (Samer)
###
############################################################

library(RColorBrewer)
library(gplots)
library(clusterProfiler)

source("customFunctions.R")
setwd("mouseMetagenome/")


#############
## Load matching metadata
df_metadataFBA <- readRDS("./df_metaData.rds")
rownames(df_metadataFBA) <- df_metadataFBA$SampleName


####################################
## Flux to Metabolite annotation
df_modelseedToMetabolite <- read.delim("metaboliteAnnotation.csv", 
                                       row.names=1, stringsAsFactors = F)
# keep only exchange fluxes for now
df_modelseedToMetabolite <- df_modelseedToMetabolite[order(rownames(df_modelseedToMetabolite), decreasing = T),]
df_modelseedToMetabolite <- df_modelseedToMetabolite[!duplicated(df_modelseedToMetabolite$IDs),]
# df_modelseedToMetabolite <- df_modelseedToMetabolite[grep("[e0]", rownames(df_modelseedToMetabolite), fixed = T),]
df_modelseedToMetabolite$metabolite <- sub(pattern = "-[cep]{1}0", replacement = "", x = df_modelseedToMetabolite$metabolite, perl = T)
# Add 1 missing item:
df_modelseedToMetabolite <- rbind(df_modelseedToMetabolite,
                                  c("Heptadecanoate","cpd15609"))
# relabel to match my currently used IDs
df_modelseedToMetabolite$Ex <- paste0("EX_",df_modelseedToMetabolite$IDs,"_e0_o")
rownames(df_modelseedToMetabolite) <- df_modelseedToMetabolite$Ex
# fix a few labels:
df_modelseedToMetabolite[df_modelseedToMetabolite$metabolite == "sprm","metabolite"] <- "Spermine"
df_modelseedToMetabolite[df_modelseedToMetabolite$metabolite == "ocdca","metabolite"] <- "Stearate"
df_modelseedToMetabolite[df_modelseedToMetabolite$metabolite == "ribflv","metabolite"] <- "Riboflavin" # Vitamine D
# 3240 entries


###################################
## Load Metabolite exchange fluxes within MB community
df_withinMBexchange <- read.table("mouseMBWithinCommunityExchangeFluxes_flux.csv", header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
df_withinMBexchange <- as.data.frame( t(df_withinMBexchange) )
# extract growthRate
df_growthRate <- df_withinMBexchange["growth",]
df_growthRate <- df_growthRate[,colnames(df_growthRate) %in% df_metadataFBA$SampleName]
df_withinMBexchange <- df_withinMBexchange[rownames(df_withinMBexchange) != "growth",]
# fix NAs, they should be considered as 0 according to Samer
df_withinMBexchange[is.na(df_withinMBexchange)] <- 0
# number of exchange fluxes
nrow(df_withinMBexchange)
# 232
# remove features with 0 or near 0 variance
caret::nearZeroVar(x = t(df_withinMBexchange), saveMetrics = T)
df_withinMBexchange <- df_withinMBexchange[-caret::nearZeroVar(x = t(df_withinMBexchange), saveMetrics = F),]
# 
nrow(df_withinMBexchange)
# 82

# Normalize by growth rate or not? - Yes!
df_withinMBexchange <- df_withinMBexchange / as.numeric(df_growthRate["growth",colnames(df_withinMBexchange)])
export(df_withinMBexchange)


####################################
## Load Metabolite exchange fluxes out or in of entire MB community
df_hostMBexchange <- read.table("mouseMBCommunityExchangeFluxes_fluxo.csv", header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
df_hostMBexchange <- as.data.frame( t(df_hostMBexchange) )
# fix NAs, they should be considered as 0 according to Samer
df_hostMBexchange[is.na(df_hostMBexchange)] <- 0
# number of exchange fluxes
nrow(df_hostMBexchange)
# 232
# remove features with 0 or near 0 variance
caret::nearZeroVar(x = t(df_hostMBexchange), saveMetrics = T)
df_hostMBexchange <- df_hostMBexchange[-caret::nearZeroVar(x = t(df_hostMBexchange), saveMetrics = F),]
# 
nrow(df_hostMBexchange)
# 51

# Normalize by growth rate or not? - Yes!
df_hostMBexchange <- df_hostMBexchange / as.numeric(df_growthRate["growth",colnames(df_hostMBexchange)])
export(df_hostMBexchange)


###################################
## Load Metabolite exchange fluxes within individual models averaged over all models of the community (internal reactions)
df_internalReactionFluxes <- read.table("mouseMBWithinModelFluxes_react.csv", header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
df_internalReactionFluxes <- as.data.frame( t(df_internalReactionFluxes) )
# # fix NAs, they should be considered as 0 according to Samer
# df_internalReactionFluxes[is.na(df_internalReactionFluxes)] <- 0
# number of exchange fluxes
nrow(df_internalReactionFluxes)
# 3222
# remove features with 0 or near 0 variance
caret::nearZeroVar(x = t(df_internalReactionFluxes), saveMetrics = T)
df_internalReactionFluxes <- df_internalReactionFluxes[-caret::nearZeroVar(x = t(df_internalReactionFluxes), saveMetrics = F),]
# dropped over 2000
# 
nrow(df_internalReactionFluxes)
# 1170

# Normalize by growth rate or not? - Yes!
df_internalReactionFluxes <- df_internalReactionFluxes / as.numeric(df_growthRate["growth",colnames(df_internalReactionFluxes)])
export(df_internalReactionFluxes)


#################
# GrowthRate is negatively Age correlated
cor.test(y = as.numeric(df_growthRate), x = df_metadataFBA[colnames(df_growthRate),"Age"], method = "spearman")
# p-value = 9.99e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.5223357 
par(mar = c(5,4,2,2)+.1)
boxplot(as.numeric(df_growthRate) ~ df_metadataFBA[colnames(df_growthRate),"Age"],
        col = brewer.pal("Pastel1",n=5), xlab = "Age [months]", ylab = "Community Growth Rate")
# 410x415


####################
## spearman correlation of MB-to-extern exchange fluxes to age
df_metadataFBA <- df_metadataFBA[colnames(df_hostMBexchange),]
all.equal(rownames(df_metadataFBA), colnames(df_hostMBexchange))
# TRUE

df_corHostMBExchangeToAge <- df_pcorAll2All(df_table1 = abs( df_hostMBexchange[,df_metadataFBA$SampleName] ),
                                            df_table2 = t(df_metadataFBA[,"Age",drop = F]),
                                            str_labels = c("ExFlux","Age"),
                                            var_pcorZ = NULL,
                                            int_cores = 7)
# [1] "Starting correlation of 1 Ages to 51 ExFluxs within 50 samples! Thu May 11 11:23:35 2023"
# [1] "Finished all correlations! Preparing results. Thu May 11 11:23:35 2023"

df_corHostMBExchangeToAge <- df_corHostMBExchangeToAge[order(df_corHostMBExchangeToAge$p.val, decreasing = F),]

# production and consumption
df_corHostMBExchangeToAge$MBconsumed <- rowSums( df_hostMBexchange[df_corHostMBExchangeToAge$ExFlux,] < 0, na.rm = T)
df_corHostMBExchangeToAge$MBproduced <- rowSums( df_hostMBexchange[df_corHostMBExchangeToAge$ExFlux,] > 0, na.rm = T)
df_corHostMBExchangeToAge$MBtoHost <- (df_corHostMBExchangeToAge$MBproduced - df_corHostMBExchangeToAge$MBconsumed) > 0

# export for CK
saveRDS(df_corHostMBExchangeToAge, paste0("df_corHostMBExchangeToAge_",Sys.Date(),".rds"))

df_tmp <- df_corHostMBExchangeToAge[df_corHostMBExchangeToAge$p.adj <= 0.05,]
df_tmp <- df_tmp[order(df_tmp$Spearman_rho),]

par(mar = c(10,4,4,2)+.1)
barplot(df_tmp$Spearman_rho, 
        names.arg = df_modelseedToMetabolite[df_tmp$ExFlux,"metabolite"], las = 2,
        ylab = "Spearman Rho", cex.names = .9,
        main = "Host x Microbiome exchange fluxes by Age", 
        col = brewer.pal("Set2", n = 3)[c(2,1)[df_tmp$MBtoHost+1]])
legend("top", legend = c("MB consumes", "MB produces"), 
       fill = brewer.pal("Set2", n = 3)[c(2,1)], bty = "n" )
# export at 490x410 barplot_spearmanRhoExFluxByAge.png

# publishable figure
df_tmp <- df_tmp[order(df_tmp$MBtoHost, df_tmp$Spearman_rho),]
str_names <- Hmisc::capitalize(df_modelseedToMetabolite[df_tmp$ExFlux,"metabolite"])
########
# Growth normalized:
str_names[30] <- "Hypoxanthine"
str_names[29] <- "Menaquinol 7"
str_names[17] <- "Fe3+"
str_names[3] <- "Mg2+"
# str_names[14] <- "Pyrophosphate"
# str_names[12] <- "Hexadecenoate"
########
# # Not normalized:
# str_names[21] <- "Pyridoxine"
# str_names[19] <- "Fe3+"
# str_names[13] <- "Mg2+"
# str_names[3] <- "Pyrophosphate"
# str_names[1] <- "Hexadecenoate"
# str_names[37] <- "Indole"
# str_names[38] <- "Hypoxanthine"
# str_names[39] <- "Methanethiol"
#
pdf(file = "barplot_spearmanRhoExFluxByAge.pdf", width = 390/72, height = 230/72, pointsize = 8)
par(mar = c(12,5,1,1)+.1)
a <- barplot(df_tmp$Spearman_rho, las = 2, horiz = F,
        names.arg = str_names, ylim = c(-.7, .7),
        space = .8, #width = .7,
        ylab = "Correlation to Age\n[Spearman Rho]", 
        cex.names = 1.1, xaxt = "n", cex.lab = 1.1, cex.axis = 1.1,
        # main = "Host x Microbiome exchange fluxes by Age", 
        col = brewer.pal("Paired", n = 12)[c(7,1)[df_tmp$MBtoHost+1]])
abline(h = 0, col = rgb(.5,.5,.5,.8), lty = 2)
axis(side = 1, at = a, lwd = 0, lwd.ticks = 1, labels = NA,
     col.ticks = rgb(.5,.5,.5,.8), lty = 3)
text(x = a, y = -.85, labels = str_names, xpd = T,
     srt = 60, adj = c(1,0.5), cex = 1.1)
legend("topleft", legend = c("MB consumes", "MB produces"), y.intersp = 1.15,
       fill = brewer.pal("Paired", n = 12)[c(7,1)], bty = "n", cex = 1.1 )
dev.off()
# export at 800x450 barplot_spearmanRhoExFluxByAge.png

# save as suppl. table
df_tmp$ExAnnotation <- str_names
write.table(x = df_tmp, file = "corHostMBExFluxByAge.csv", 
            quote = T, sep = ",", row.names = F, col.names = T)
rm(df_tmp, a)


############################
## Spearman correlation of MB-internal exchange fluxes to age
df_metadataFBA <- df_metadataFBA[colnames(df_withinMBexchange),]
all.equal(rownames(df_metadataFBA), colnames(df_withinMBexchange))
# TRUE
df_corWithinMBExchangeToAge <- df_pcorAll2All(df_table1 = abs( df_withinMBexchange[,df_metadataFBA$SampleName] ),
                                            df_table2 = t(df_metadataFBA[,"Age",drop = F]),
                                            str_labels = c("ExFlux","Age"),
                                            var_pcorZ = NULL,
                                            int_cores = 7)
# [1] "Starting correlation of 1 Ages to 82 ExFluxs within 50 samples! Tue May 16 09:35:59 2023"
# [1] "Finished all correlations! Preparing results. Tue May 16 09:35:59 2023"

df_corWithinMBExchangeToAge <- df_corWithinMBExchangeToAge[order(df_corWithinMBExchangeToAge$p.val, decreasing = F),]

# how many samples  
df_corWithinMBExchangeToAge$SampleCount <- rowSums( df_withinMBexchange[df_corWithinMBExchangeToAge$ExFlux,] > 0, na.rm = T)

df_tmp <- df_corWithinMBExchangeToAge[df_corWithinMBExchangeToAge$p.adj <= 0.05,]
df_tmp <- df_tmp[order(df_tmp$Spearman_rho),]

par(mar = c(4,12,4,2)+.1)
barplot(df_tmp$Spearman_rho, horiz = T,
        names.arg = df_modelseedToMetabolite[df_tmp$ExFlux,"metabolite"], las = 2,
        xlab = "Spearman Rho", cex.names = .9,
        main = "MB x Microbiome exchange fluxes by Age") # col = brewer.pal("Set2", n = 3)[c(2,1)[df_tmp$MBtoHost+1]])
# legend("bottomright", legend = c("MB consumes", "MB produces"), 
#        fill = brewer.pal("Set2", n = 3)[c(2,1)], bty = "n" )
# export at 420x760 barplot_spearmanRhoExFluxByAge.png

pdf(file = "barplot_corSummaryWithinCommunityExFluxByAge.pdf", width = 140/72, height = 150/72, pointsize = 10)
par(mar = c(5,2,2,2)+.1)
barplot(table(df_tmp$Spearman_rho > 0), horiz = F,
        names.arg = c("Reduced","Increased"), las = 1,
        xlab = "Aging correlation", cex.names = 1.1,
        cex.lab = 1.2, ylim = c(0,40), width = .5, space = 1.5,
        col = brewer.pal("Paired", n = 12)[c(1)], 
        ); abline(h = 0, col = rgb(.5,.5,.5,.8), lty = 2)
dev.off()
# 250x300


# publishable figure
str_names <- Hmisc::capitalize( df_modelseedToMetabolite[df_tmp$ExFlux,"metabolite"] )
########
# Growth normalized:
str_names[8] <- "Niacine"
str_names[3] <- "Betaine"
str_names[13] <- "Methanethiol"
str_names[16] <- "Guanosine"
str_names[23] <- "Folate"
str_names[25] <- "Deoxyguanosine"
#
pdf(file = "barplot_spearmanRhoWithinCommunityExFluxByAge.pdf", width = 800/72, height = 450/72, pointsize = 14)
par(mar = c(12,5,1,1)+.1)
a <- barplot(df_tmp$Spearman_rho, las = 2, horiz = F,
             names.arg = str_names, ylim = c(-.7, .7),
             space = .8, #width = .7,
             ylab = "Correlation to Age\n[Spearman Rho]", 
             cex.names = 1.1, xaxt = "n", cex.lab = 1.1, cex.axis = 1.1,
             # main = "Host x Microbiome exchange fluxes by Age", 
             col = brewer.pal("Paired", n = 12)[c(1)])
abline(h = 0, col = rgb(.5,.5,.5,.8), lty = 2)
axis(side = 1, at = a, lwd = 0, lwd.ticks = 1, labels = NA,
     col.ticks = rgb(.5,.5,.5,.8), lty = 3)
text(x = a, y = -.85, labels = str_names, xpd = T,
     srt = 60, adj = c(1,0.5), cex = 1.1)
dev.off()
# legend("topleft", legend = c("MB consumes", "MB produces"), y.intersp = 1.5,
#        fill = brewer.pal("Paired", n = 12)[c(7,1)], bty = "n", cex = 1.1 )
# export at 800x450 barplot_spearmanRhoWithinCommunityExFluxByAge.png

# save as suppl. table
df_tmp$ExAnnotation <- str_names
write.table(x = df_tmp, file = "corWithinCommunityExFluxByAge.csv", 
            quote = T, sep = ",", row.names = F, col.names = T)
rm(df_tmp, a)


############################
## Spearman correlation of within models reaction fluxes to age
df_tmpMetadataFBA <- df_metadataFBA[colnames(df_internalReactionFluxes),]
all.equal(rownames(df_tmpMetadataFBA), colnames(df_internalReactionFluxes))
# TRUE
df_corInternRxnFluxToAge <- df_pcorAll2All(df_table1 = abs( df_internalReactionFluxes[,df_tmpMetadataFBA$SampleName] ),
                                          df_table2 = t(df_tmpMetadataFBA[,"Age",drop = F]),
                                          str_labels = c("RxnFlux","Age"),
                                          var_pcorZ = NULL,
                                          int_cores = 7)
# [1] "Starting correlation of 1 Ages to 1170 RxnFluxs within 50 samples! Mon May 22 12:39:29 2023"
# [1] "Finished all correlations! Preparing results. Mon May 22 12:39:29 2023"

df_corInternRxnFluxToAge <- df_corInternRxnFluxToAge[order(df_corInternRxnFluxToAge$p.val, decreasing = F),]

# direction of reaction in models in % of all samples
df_corInternRxnFluxToAge$FluxDirectionPosPercent <- rowSums( df_internalReactionFluxes[df_corInternRxnFluxToAge$RxnFlux,] > 0, na.rm = T)/ncol(df_internalReactionFluxes)*100
df_corInternRxnFluxToAge$FluxDirectionNegPercent <- rowSums( df_internalReactionFluxes[df_corInternRxnFluxToAge$RxnFlux,] < 0, na.rm = T)/ncol(df_internalReactionFluxes)*100
df_corInternRxnFluxToAge$FluxDirectionPercent <- df_corInternRxnFluxToAge$FluxDirectionPos - df_corInternRxnFluxToAge$FluxDirectionNeg

df_tmp <- df_corInternRxnFluxToAge[df_corInternRxnFluxToAge$p.adj <= 0.01 &
                                   abs(df_corInternRxnFluxToAge$Spearman_rho) >= 0.5,]
df_tmp <- df_tmp[order(df_tmp$Spearman_rho),]

par(mar = c(3,20,1,1)+.1)
barplot(df_tmp$Spearman_rho, horiz = T,
        names.arg = df_rxnAnnotation[df_tmp$RxnFlux,"RxnName"], las = 2,
        xlab = "Spearman Rho", cex.names = .7,
        border = 0,
        # main = "MB x Microbiome exchange fluxes by Age"
        col = colorRampPalette(brewer.pal("RdBu", n = 11))(200)[df_tmp$FluxDirectionPercent+100] 
        )
# legend("bottomright", legend = c("MB consumes", "MB produces"), 
#        fill = brewer.pal("Set2", n = 3)[c(2,1)], bty = "n" )
# export at 650x1600 barplot_topInternalRxnFluxesByAge

# save as suppl. table
df_tmp <- df_corInternRxnFluxToAge[df_corInternRxnFluxToAge$p.adj <= 0.05,]
df_tmp$RxnAnnotation <- df_rxnAnnotation[df_tmp$RxnFlux,"RxnName"]
df_tmp <- df_tmp[,c("RxnFlux","RxnAnnotation","Spearman_rho","p.val","p.adj","FluxDirectionPosPercent","FluxDirectionNegPercent","FluxDirectionPercent")]
write.table(x = df_tmp, file = "corInternRxnFluxToAge.csv", 
            quote = T, sep = ",", row.names = F, col.names = T)
rm(df_tmp)


#####
### Subsystem enrichment
## Positive Correlations Reactions
obj_Rxn2SubsysPositive <- enricher(gene = unique(df_corInternRxnFluxToAge[df_corInternRxnFluxToAge$p.adj <= 0.05 &
                                                                          df_corInternRxnFluxToAge$Spearman_rho > 0,"RxnFlux"]),
                                   minGSSize = 3, maxGSSize = 500, pAdjustMethod = "BH",
                                   TERM2GENE = df_rxn2subsys[,c(1,2)],
                                   TERM2NAME = df_rxn2subsys[,c(1,3)],
                                   universe = rownames(df_internalReactionFluxes))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Rxn2SubsysPositive <- obj_Rxn2SubsysPositive@result[order(obj_Rxn2SubsysPositive@result$qvalue,
                                                             decreasing = F),]
# df_Rxn2SubsysPositiveColon <- df_Rxn2SubsysPositive[df_Rxn2SubsysPositive$p.adjust <= 0.05 &
#                                                df_Rxn2SubsysPositive$Count >= 3,]
df_Rxn2SubsysPositive <- df_Rxn2SubsysPositive[df_Rxn2SubsysPositive$pvalue <= 0.05 &
                                               df_Rxn2SubsysPositive$Count >= 3,]
unique( str_simplifyDescriptions( df_Rxn2SubsysPositive$Description ) )
# [1] "thymine degr."                     "uracil degr."                      "L-isoleucine synth."               "L-valine synth."                  
# [5] "purine deoxyribonucleosides degr." "adenine and adenosine salvage"     "purine ribonucleosides degr."  

## Negative Correlations Reactions
obj_Rxn2SubsysNegative <- enricher(gene = unique(df_corInternRxnFluxToAge[df_corInternRxnFluxToAge$p.adj <= 0.05 &
                                                                            df_corInternRxnFluxToAge$Spearman_rho < 0,"RxnFlux"]),
                                   minGSSize = 3, maxGSSize = 500, pAdjustMethod = "BH",
                                   TERM2GENE = df_rxn2subsys[,c(1,2)],
                                   TERM2NAME = df_rxn2subsys[,c(1,3)],
                                   universe = rownames(df_internalReactionFluxes))
# Reorder result by GO-term size from smaller (= more detailed) to larger GO-groups
df_Rxn2SubsysNegative <- obj_Rxn2SubsysNegative@result[order(obj_Rxn2SubsysNegative@result$qvalue,
                                                             decreasing = F),]
# df_Rxn2SubsysNegativeColon <- df_Rxn2SubsysNegative[df_Rxn2SubsysNegative$p.adjust <= 0.05 &
#                                                df_Rxn2SubsysNegative$Count >= 3,]
df_Rxn2SubsysNegative <- df_Rxn2SubsysNegative[df_Rxn2SubsysNegative$pvalue <= 0.05 &
                                                 df_Rxn2SubsysNegative$Count >= 3,]
unique( str_simplifyDescriptions( df_Rxn2SubsysNegative$Description ) )
# [1] "PEPTDOGLYCANSYN-PWY2"                                                                 
# [2] "CDP-diacylglycerol synth."                                                            
# [3] "peptidoglycan maturation"                                                             
# [4] "UDP-N-acetylmuramoyl-pentapeptide synth."                                             
# [5] "phosphatidylserine and phosphatidylethanolamine synth."                               
# [6] "tetrapyrrole synth."                                                                  
# [7] "UMP synth."                                                                           
# [8] "methylerythritol phosphate pwy."                                                      
# [9] "adenosylcobalamin synth. from adenosylcobinamide-GDP"                                 
# [10] "benzimidazolyl adenosylcobamide\nsynth. from adenosylcobinamide-GDP"                  
# [11] "2-methyladeninyl adenosylcobamide\nsynth. from adenosylcobinamide-GDP"                
# [12] "heme b synth."                                                                        
# [13] "L-serine synth."                                                                      
# [14] "PWY-8088"                                                                             
# [15] "5-hydroxybenzimidazolyl adenosylcobamide\nsynth. from adenosylcobinamide-GDP"         
# [16] "5-methoxy-6-methylbenzimidazolyl\nadenosylcobamide synth. from adenosylcobinamide-GDP"
# [17] "5-methoxybenzimidazolyl adenosylcobamide\nsynth. from adenosylcobinamide-GDP"         
# [18] "reductive acetyl coenzyme A pwy."                                                     
# [19] "5-methylbenzimidazolyl adenosylcobamide\nsynth. from adenosylcobinamide-GDP"          
# [20] "adeninyl adenosylcobamide synth. from adenosylcobinamide-GDP"                         
# [21] "L-lysine synth."                                                                      
# [22] "sphingosine and sphingosine-1-phosphate metab."                                       
# [23] "folate synth."                                                                        
# [24] "phytol degr."                                                                       
# [25] "anhydromuropeptides recycling"    

#####
## Make a heatmap of mean rxn activity from each GO-group 
df_Rxn2Subsys <- rbind(df_Rxn2SubsysPositive[df_Rxn2SubsysPositive$pvalue <= 0.05 &
                                             df_Rxn2SubsysPositive$Count >= 3,],
                       df_Rxn2SubsysNegative[df_Rxn2SubsysNegative$pvalue <= 0.05 &
                                             df_Rxn2SubsysNegative$Count >= 3,]
                       )
# mtx_mydist(df_Rxn2SubsysNegative)
mtx_tmp2 <- matrix(data = NA, nrow = nrow(df_Rxn2Subsys), ncol = 5)
rownames(mtx_tmp2) <- df_Rxn2Subsys[hclust(as.dist(mtx_mydist(df_Rxn2Subsys)))$order,"Description"]
for(str_GOId in df_Rxn2Subsys$ID) {
  lst_rxns <- unlist(strsplit(x = df_Rxn2Subsys[str_GOId,"geneID"], "/", fixed = T))
  # take the absolute of rxns (disregarding directions) ? 
  # mtx_tmp <- abs( df_internalReactionFluxes[df_corInternRxnFluxToAge[df_corInternRxnFluxToAge$p.adj <= 0.05,"RxnFlux"],
  #                                           df_tmpMetadataFBA$SampleName] )
  # or not absolute?
  mtx_tmp <- df_internalReactionFluxes[df_corInternRxnFluxToAge[df_corInternRxnFluxToAge$p.adj <= 0.05,"RxnFlux"],
                                       df_tmpMetadataFBA$SampleName]
  mtx_tmp <- as.matrix(mtx_tmp[lst_rxns,order(df_tmpMetadataFBA[colnames(mtx_tmp),"AgeGroup"])])
  
  # normalise matrix so each reaction is similarily abundant (percentage of rxn-abundance)
  mtx_tmp <- mtx_tmp / rowSums(mtx_tmp) * 1000
  
  mtx_tmp2[df_Rxn2Subsys[str_GOId,"Description"],] <- c( mean( mtx_tmp[, df_tmpMetadataFBA[colnames(mtx_tmp),"AgeGroup"] == 1 ] ),
                                                         mean( mtx_tmp[, df_tmpMetadataFBA[colnames(mtx_tmp),"AgeGroup"] == 2 ] ),
                                                         mean( mtx_tmp[, df_tmpMetadataFBA[colnames(mtx_tmp),"AgeGroup"] == 3 ] ),
                                                         mean( mtx_tmp[, df_tmpMetadataFBA[colnames(mtx_tmp),"AgeGroup"] == 4 ] ),
                                                         mean( mtx_tmp[, df_tmpMetadataFBA[colnames(mtx_tmp),"AgeGroup"] == 5 ] )
                                                        )
}

# export enrichment to supplements
df_tmp <- merge(x = df_Rxn2Subsys, y = mtx_tmp2,
                by.x = "Description", by.y = "row.names")
df_tmp <- df_tmp[,c(2,1,3:5,8,10:14)]
colnames(df_tmp) <- c("TermID","SubSystem","FeatureRatio","BackgroundRatio","pvalue","FeatureIDs","Mean2months","Mean9months","Mean15months","Mean24months","Mean30months")
write.table(x = df_tmp, file = "enrichmentAgeAssocInternRxnFlux.csv", 
            quote = T, sep = ",", row.names = F, col.names = T)
rm(df_tmp)


str_rowLabels <- str_simplifyDescriptions( rownames(mtx_tmp2) )
str_rowLabels <- Hmisc::capitalize(str_rowLabels)

rownames(mtx_tmp2) <- str_rowLabels
# drop duplicates by value
duplicated(mtx_tmp2)
# mtx_tmp2 <- mtx_tmp2[rownames(mtx_tmp2) != "Cis-genanyl-CoA degr.",] # is a duplicate of "Acetate conversion to acetyl-CoA "
mtx_tmp2 <- mtx_tmp2[!duplicated(mtx_tmp2),]
# drop duplicates by name
rownames(mtx_tmp2)[duplicated(rownames(mtx_tmp2))]
mtx_tmp2 <- mtx_tmp2[!duplicated(rownames(mtx_tmp2)),]

# simplify rownames a bit more:
rownames(mtx_tmp2)[1] <- "Folate synth."
rownames(mtx_tmp2)[2] <- "Sphingosine metab."
rownames(mtx_tmp2)[4] <- "Reductive acetyl-CoA pwy."
rownames(mtx_tmp2)[9] <- "Adeninyl adenosylcobamide synth."
rownames(mtx_tmp2)[13] <- "Phosphatidyl-serine & -ethanolamine synth."

# create a heatmap plot
pdf(file = "heatmap_MouseMicrobiomeModelsInternRxnFluxMeanSubSysAbundanceByAge.pdf", width = 700/72, height = 1100/72, pointsize = 14)
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

rm(df_tmpMetadataFBA, obj_Rxn2SubsysPositive, df_Rxn2SubsysPositive, obj_Rxn2SubsysNegative, df_Rxn2SubsysNegative, str_rowLabels, mtx_tmp, mtx_tmp2, df_Rxn2Subsys)


