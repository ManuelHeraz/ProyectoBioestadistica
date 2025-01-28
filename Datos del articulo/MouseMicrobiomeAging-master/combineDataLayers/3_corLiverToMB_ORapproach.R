############################################################
###
###   Combine data layers with correlations or multiOmics
###
############################################################

library(RColorBrewer)
library(DESeq2)
library(gplots)
library(clusterProfiler)
library('biomaRt')
library(parallel)

source("customFunctions.R")

setwd("combineDataLayers/")


################################
## Aging regulated host genes
# Variance stabilized RNA-Seq reads informed with Age as integer
obj_DifAbundLiverAllAges <- readRDS(file = "../mouseTranscriptome/obj_DifAbundLiverAllAges.rds")
df_liverAgeDependance <- as.data.frame( results(obj_DifAbundLiverAllAges, name = "age", alpha = 0.05) )
rm(obj_DifAbundLiverAllAges)
df_liverAgeDependance <- df_liverAgeDependance[df_liverAgeDependance$padj <= 0.05 & 
                                                 !is.na(df_liverAgeDependance$padj),]
## add GeneSymbol, req. for GO-annotation
df_liverAgeDependance$GeneSymbol <- df_mm10Ensembl2Genesymbol[rownames(df_liverAgeDependance),"GeneSymbol"]
# some are missing - help out with biomart online tools
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
df_tmp <- getBM( values = rownames(df_liverAgeDependance)[is.na(df_liverAgeDependance$GeneSymbol)],
                 filters = "ensembl_gene_id", 
                 attributes = c("ensembl_gene_id","mgi_symbol"),
                 mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl")) )
rownames(df_tmp) <- df_tmp$ensembl_gene_id
# add remaining GeneSymbols:
df_liverAgeDependance[is.na(df_liverAgeDependance$GeneSymbol),"GeneSymbol"] <- df_tmp[rownames(df_liverAgeDependance)[is.na(df_liverAgeDependance$GeneSymbol)],"mgi_symbol"]
rm(df_tmp)
# still one missing
df_liverAgeDependance["ENSMUSG00000097971","GeneSymbol"] <-  "Gm26917"
df_liverAgeDependance["ENSMUSG00000074876","GeneSymbol"] <-  "Spata5l1"
# order by p-val
df_liverAgeDependance <- df_liverAgeDependance[order(df_liverAgeDependance$pvalue, decreasing = F),]


################################
## partial Correlation results:
nrow(df_pcorRxn2LiverRNATop)
# 3425

## add GeneSymbol, req. for GO-annotation
df_pcorRxn2LiverRNATop$GeneSymbol <- df_mm10Ensembl2Genesymbol[df_pcorRxn2LiverRNATop$Gene,"GeneSymbol"]
# some are missing
unique(df_pcorRxn2LiverRNATop[is.na(df_pcorRxn2LiverRNATop$GeneSymbol),"Gene"])
# help out with biomart online tools
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
df_tmp <- getBM( values = unique(df_pcorRxn2LiverRNATop[is.na(df_pcorRxn2LiverRNATop$GeneSymbol),"Gene"]),
                 filters = "ensembl_gene_id", 
                 attributes = c("ensembl_gene_id","mgi_symbol"),
                 mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl")) )
rownames(df_tmp) <- df_tmp$ensembl_gene_id
# add remaining GeneSymbols:
df_pcorRxn2LiverRNATop[is.na(df_pcorRxn2LiverRNATop$GeneSymbol),"GeneSymbol"] <- df_tmp[df_pcorRxn2LiverRNATop[is.na(df_pcorRxn2LiverRNATop$GeneSymbol),"Gene"],"mgi_symbol"]
rm(df_tmp)
# one still missing:
df_pcorRxn2LiverRNATop[df_pcorRxn2LiverRNATop$Gene == "ENSMUSG00000115902","GeneSymbol"] <- "D730005E14Rik"


# Get GO bio processes for each gene in our correlation list
df_liverRNATopAllGOBioProcess <- df_mouseGOannotaionBioProc[df_mouseGOannotaionBioProc$GeneSymbol %in% unique(df_pcorRxn2LiverRNATop$GeneSymbol),c(1,2,6)]
df_liverRNATopAllGOBioProcess <- df_liverRNATopAllGOBioProcess[!duplicated(df_liverRNATopAllGOBioProcess),]
# concatenate genes into a single row per GO-term
df_tmp <- as.data.frame( tapply(X = df_liverRNATopAllGOBioProcess, INDEX = df_liverRNATopAllGOBioProcess$GOid, FUN = function(df_tmp) {
  return( paste0(df_tmp$GeneSymbol, collapse = ",") )
}) )
colnames(df_tmp) <- "Features"
df_liverRNATopAllGOBioProcess <- merge(x = df_tmp, y = df_liverRNATopAllGOBioProcess[,c(1,3)],
                                       by.x = "row.names", by.y = "GOid", all.x = T, all.y = F)
colnames(df_liverRNATopAllGOBioProcess) <- c("ID", "Features", "Description")
df_liverRNATopAllGOBioProcess <- df_liverRNATopAllGOBioProcess[!duplicated(df_liverRNATopAllGOBioProcess),]
rownames(df_liverRNATopAllGOBioProcess) <- df_liverRNATopAllGOBioProcess$ID
# feature count per GO-term (only features found in our correlation-data)
df_liverRNATopAllGOBioProcess$Count <- stringi::stri_count(str = df_liverRNATopAllGOBioProcess$Features, fixed = ",")+1
# keep only processes with at least 3 features:
df_liverRNATopAllGOBioProcess <- df_liverRNATopAllGOBioProcess[df_liverRNATopAllGOBioProcess$Count >= 3,]
# Drop certain terms like anything with "embryo*" in it for example
df_liverRNATopAllGOBioProcess <- df_liverRNATopAllGOBioProcess[grep(pattern = "embryo*", x = df_liverRNATopAllGOBioProcess$Description, ignore.case = T, invert = T),]
# add background size of term
df_liverRNATopAllGOBioProcess$TermSize <- tapply(X = df_mouseGOannotaionBioProc$GOid, INDEX = df_mouseGOannotaionBioProc$GOid, FUN = length)[df_liverRNATopAllGOBioProcess$ID]
# calculate overrepresentation of gene sets
int_nTotalInputGenes <- length(unique(unlist(strsplit(x = df_liverRNATopAllGOBioProcess$Features, split = ",", fixed = T))))
int_nTotalBackgroundGenes <- length(unique(df_mouseGOannotaionBioProc$GeneSymbol))
df_liverRNATopAllGOBioProcess$EnrichPval <- unlist( lapply(X = 1:nrow(df_liverRNATopAllGOBioProcess),
                                                           FUN = function(int_i) {
                                                             return( phyper(df_liverRNATopAllGOBioProcess[int_i,"Count"]-1, 
                                                                            int_nTotalInputGenes, 
                                                                            int_nTotalBackgroundGenes-int_nTotalInputGenes, 
                                                                            df_liverRNATopAllGOBioProcess[int_i,"TermSize"],
                                                                            lower.tail= FALSE) )
                                                           }))
df_liverRNATopAllGOBioProcess$EnrichPadj <- p.adjust(p = df_liverRNATopAllGOBioProcess$EnrichPval, method = "BH")
# NOTE: It makes sense to filter processes on enrichment p-val.,
#  since we don't have a defined pre-selection as we would have with metabolic modelling (MB side, or metamodel)
table(df_liverRNATopAllGOBioProcess$EnrichPadj <= 0.05)
# FALSE  TRUE 
# 618   115 
df_liverRNATopAllGOBioProcess <- df_liverRNATopAllGOBioProcess[df_liverRNATopAllGOBioProcess$EnrichPadj <= 0.05,]


## Get subsystems for each rxn in our correlation list
df_liverRNATopAllSubSys <- df_rxn2subsys[df_rxn2subsys$ReactionID %in% unique(df_pcorRxn2LiverRNATop$Rxn),c(1:3)]
df_liverRNATopAllSubSys <- df_liverRNATopAllSubSys[!duplicated(df_liverRNATopAllSubSys),]
# concatenate rxns into a single row per subsys
df_tmp <- as.data.frame( tapply(X = df_liverRNATopAllSubSys, INDEX = df_liverRNATopAllSubSys$SubsysID, FUN = function(df_tmp) {
  return( paste0(df_tmp$ReactionID, collapse = ",") )
}) )
colnames(df_tmp) <- "Features"
df_liverRNATopAllSubSys <- merge(x = df_tmp, y = df_liverRNATopAllSubSys[,c(1,3)],
                                 by.x = "row.names", by.y = "SubsysID", all.x = T, all.y = F)
colnames(df_liverRNATopAllSubSys) <- c("ID", "Features", "Description")
df_liverRNATopAllSubSys <- df_liverRNATopAllSubSys[!duplicated(df_liverRNATopAllSubSys),]
rownames(df_liverRNATopAllSubSys) <- df_liverRNATopAllSubSys$ID
# feature count per GO-term (only features found in our correlation-data)
df_liverRNATopAllSubSys$Count <- stringi::stri_count(str = df_liverRNATopAllSubSys$Features, fixed = ",")+1
# keep only processes with at least 3 features:
df_liverRNATopAllSubSys <- df_liverRNATopAllSubSys[df_liverRNATopAllSubSys$Count >= 3,]


### overrepresentation tests for each process-pair (GO and Subsys) in parallel
int_nTotalStrongCorrs <- nrow(df_pcorRxn2LiverRNATop)
obj_clusters <- makeCluster(7, type="SOCK")
clusterExport(cl = obj_clusters, varlist = c("df_liverRNATopAllGOBioProcess","df_liverRNATopAllSubSys","df_pcorRxn2LiverRNATop","int_nTotalStrongCorrs"), envir = environment())
date()
# [1] "Tue Aug  8 11:22:31 2023"
df_overRepLiverToMBprocesses <- as.data.frame(parLapply(cl = obj_clusters, 
                                                        X = 1:nrow(df_liverRNATopAllGOBioProcess), 
                                                        fun = function(int_i) {
                                                          return( lapply(X = 1:nrow(df_liverRNATopAllSubSys), FUN = function(int_j) {
                                                            str_geneList <- unlist(strsplit(x = df_liverRNATopAllGOBioProcess[int_i,"Features"], split = ",", fixed = T))
                                                            str_rxnList <- unlist(strsplit(x = df_liverRNATopAllSubSys[int_j,"Features"], split = ",", fixed = T))
                                                            
                                                            # sign. correlations overlapping in this GO processes and this SubSystem
                                                            int_countOverlapCors <- sum( (df_pcorRxn2LiverRNATop$GeneSymbol %in% str_geneList) &
                                                                                           (df_pcorRxn2LiverRNATop$Rxn %in% str_rxnList) )
                                                            # total sign. correlation found
                                                            int_countGOBackgroundCors <- sum( df_pcorRxn2LiverRNATop$GeneSymbol %in% str_geneList )
                                                            int_countSubSysBackgroundCors <- sum( df_pcorRxn2LiverRNATop$Rxn %in% str_rxnList )
                                                            
                                                            if(int_countOverlapCors > 0) {
                                                              # print(paste(int_i, int_j))
                                                              # print(c(length(dbl_rho) - 1, 
                                                              #         df_liverRNATopAllSubSys[1,"Count"],
                                                              #         int_nTotalStrongCorrs - df_liverRNATopAllSubSys[1,"Count"],
                                                              #         df_liverRNATopAllGOBioProcess[int_i,"Count"]))
                                                              
                                                              # correlation values for overlapping MB and Host features
                                                              dbl_rho <- mean( df_pcorRxn2LiverRNATop[df_pcorRxn2LiverRNATop$GeneSymbol %in% str_geneList &
                                                                                                        df_pcorRxn2LiverRNATop$Rxn %in% str_rxnList,"Spearman_rho"] )
                                                              
                                                              # overrepresentation p-value (more strong correlations in this certain subset (combination of 1 GO with 1 Subsys) than expected by chance)
                                                              dbl_pOR <-  phyper(int_countOverlapCors - 1, 
                                                                                 int_countSubSysBackgroundCors,
                                                                                 int_nTotalStrongCorrs - int_countSubSysBackgroundCors,
                                                                                 int_countGOBackgroundCors,
                                                                                 lower.tail= F,
                                                                                 log.p = T)
                                                              
                                                              # return results
                                                              return( c(df_liverRNATopAllSubSys[int_j,"ID"], # SubsysID
                                                                        df_liverRNATopAllGOBioProcess[int_i,"ID"],  # GOID
                                                                        dbl_pOR, # overrepresentation log(p-val)
                                                                        dbl_rho, # mean correlation value
                                                                        int_countOverlapCors,  # overlap size
                                                                        df_liverRNATopAllSubSys[int_j,"Count"], # group 1 background size
                                                                        df_liverRNATopAllGOBioProcess[int_i,"Count"], # group 2 background size
                                                                        int_countSubSysBackgroundCors, # total sign. cors of rxns of this Subsys
                                                                        int_countGOBackgroundCors # total sign. cors of features of this GO
                                                              ) )
                                                            }
                                                            else return(c(NA,NA,NA,NA,NA,NA,NA,NA,NA))
                                                          }) )
                                                        }) )
date()
# [1] "Tue Aug  8 11:22:42 2023"
# clean up
stopCluster(obj_clusters)
# convert to data frame and clean up
df_overRepLiverToMBprocesses <- as.data.frame(t(df_overRepLiverToMBprocesses))
df_overRepLiverToMBprocesses <- df_overRepLiverToMBprocesses[!is.na(df_overRepLiverToMBprocesses$V3),]
rownames(df_overRepLiverToMBprocesses) <- NULL
colnames(df_overRepLiverToMBprocesses) <- c("SubsysID","GOID", "logPval", "MeanRho", "Overlap", "SubSysBackground", "GOBackground","SubSysSignCors","GOSignCors")
# make numeric columns actually numeric
for(int_i in 3:9) {
  df_overRepLiverToMBprocesses[,int_i] <- as.numeric(df_overRepLiverToMBprocesses[,int_i])
}
df_overRepLiverToMBprocesses$p.adj <- p.adjust(exp(df_overRepLiverToMBprocesses$logPval), method = "BH")
df_overRepLiverToMBprocesses <- df_overRepLiverToMBprocesses[order(df_overRepLiverToMBprocesses$logPval, decreasing = F),]
# calculate association score:
df_overRepLiverToMBprocesses$Score <- df_overRepLiverToMBprocesses$Overlap / (df_overRepLiverToMBprocesses$SubSysBackground * df_overRepLiverToMBprocesses$GOBackground)
# Store results
saveRDS(df_overRepLiverToMBprocesses,"df_overRepLiverToMBprocesses.rds")

table(df_overRepLiverToMBprocesses$p.adj <= 0.01)
# FALSE  TRUE 
# 3284   428

# add GO infos
df_tmp <- merge(x = df_overRepLiverToMBprocesses, y = df_liverRNATopAllGOBioProcess[,1:3],
                by.x = "GOID", by.y = "ID", all.x = T, all.y = F)
colnames(df_tmp) <- c("GOID", "SubsysID", "logPval", "MeanRho", "Overlap", "SubSysBackground", "GOBackground", "SubSysSignCors",
                      "GOSignCors", "p.adj", "Score", "Genes", "GOTerm")
df_overRepLiverToMBprocesses <- merge(x = df_tmp, y = df_liverRNATopAllSubSys[,1:3],
                                      by.x = "SubsysID", by.y = "ID", all.x = T, all.y = F)
colnames(df_overRepLiverToMBprocesses) <- c("SubsysID", "GOID", "logPval", "MeanRho", "Overlap", "SubSysBackground", "GOBackground", "SubSysSignCors",
                                            "GOSignCors", "p.adj", "Score", "Genes", "GOTerm", "Rxns", "SubSystem")

# filter by pval
df_overRepLiverToMBprocessesTop <- df_overRepLiverToMBprocesses[df_overRepLiverToMBprocesses$p.adj <= 0.01,]

## check for duplicate entries
# simplified process names
df_overRepLiverToMBprocessesTop$GOSimple <- str_simplifyDescriptions( df_overRepLiverToMBprocessesTop$GOTerm, bol_stripRomanNums = F )
df_overRepLiverToMBprocessesTop$SubSysSimple <- str_simplifyDescriptions( df_overRepLiverToMBprocessesTop$SubSystem, bol_stripRomanNums = T )
# exact duplicates (usually due to redudandant subsystem names)
table(duplicated(df_overRepLiverToMBprocessesTop[,c(2:14,17)]))
df_overRepLiverToMBprocessesTop <- df_overRepLiverToMBprocessesTop[!duplicated(df_overRepLiverToMBprocessesTop[,c(2:14,17)]),]
# dropped 17

# same simple names, choose better logP (min)
table(duplicated(df_overRepLiverToMBprocessesTop[,c("GOSimple","SubSysSimple")]))
df_overRepLiverToMBprocessesTop <- df_overRepLiverToMBprocessesTop[order(df_overRepLiverToMBprocessesTop$logPval, decreasing = F),]
df_overRepLiverToMBprocessesTop <- df_overRepLiverToMBprocessesTop[!duplicated(df_overRepLiverToMBprocessesTop[,c("GOSimple","SubSysSimple")]),]
# dropped 62

# similar correlation results but (slightly) different subsys name 
table(duplicated(df_overRepLiverToMBprocessesTop[,c(2:14)]))
df_overRepLiverToMBprocessesTop$DupID <-  as.numeric( as.factor(unlist( apply(X = df_overRepLiverToMBprocessesTop[,c(2:14)], MARGIN = 1, function(x) {
  return( paste(collapse = "|", x) )
}))))
#  merge those into one single column
df_tmp <- as.data.frame( matrix( data = unlist( tapply(X = 1:nrow(df_overRepLiverToMBprocessesTop), INDEX = df_overRepLiverToMBprocessesTop$DupID, 
                                                       FUN = function(int_i) {
                                                         if(length(int_i) == 1) return(c(df_overRepLiverToMBprocessesTop[int_i,,drop = T],df_overRepLiverToMBprocessesTop[int_i,"SubsysID"]))
                                                         else {
                                                           # select one representative (shortest label)
                                                           df_tmp <- df_overRepLiverToMBprocessesTop[int_i,]
                                                           df_tmp <- df_tmp[order(stringi::stri_length(df_tmp$SubSysSimple), decreasing = F),]
                                                           # concatenate subsys Descriptions and Ids
                                                           df_tmp[1,"SubSysIDs"] <- paste(collapse = ",", df_tmp$SubsysID)
                                                           df_tmp[1,"SubSystem"] <- paste(collapse = ",", df_tmp$SubSystem)
                                                           
                                                           return(c(df_tmp[1,,drop=T]))
                                                         }
                                                       }, simplify = T)), ncol = 19, byrow = T) )
colnames(df_tmp) <- c(colnames(df_overRepLiverToMBprocessesTop),"SubSysIDs")
df_tmp$DupID <- NULL
df_overRepLiverToMBprocessesTop <- df_tmp; rm(df_tmp)
# make numeric columns actually numeric
for(int_i in 3:11) {
  df_overRepLiverToMBprocessesTop[,int_i] <- as.numeric(df_overRepLiverToMBprocessesTop[,int_i])
}
# merged 66 subsys

# duplicated GO terms?
table(duplicated(df_overRepLiverToMBprocessesTop[,c(1,3:14,16)]))
# -none-

# Take into account the direction of correlation
df_overRepLiverToMBprocessesTop$SignedLogP <- log10(df_overRepLiverToMBprocessesTop$p.adj) * (((df_overRepLiverToMBprocessesTop$MeanRho < 0)*2)-1)


################################
## Export to table - Liver Correlation results table S3
# View(df_overRepLiverToMBprocessesTop[,c("GOID", "SubSysIDs", "p.adj", "SignedLogP", "GOTerm", "SubSystem", "Overlap", "SubSysBackground", "GOBackground", "Genes", "Rxns" )])
# Export table in a tabular format for people to explore as excel-table
df_liverAssocScoresOverview <- df_overRepLiverToMBprocessesTop[,c("GOID", "SubSysIDs", "p.adj", "SignedLogP", "GOTerm", "SubSystem", "Overlap", "SubSysBackground", "GOBackground", "Genes", "Rxns" )]
colnames(df_liverAssocScoresOverview) <- c("HOST_TermID", "MB_TermID", "EnrichPadj", "SignedLogPadj", "HOST_BioProcess", "MB_SubSystem", 
                                           "CorrFeatures", "MB_FeatureCount", "HOST_FeatureCount", "HOST_EnrichedFeatures", "MB_EnrichedFeatureIDs")

### ADD higher level group
# host and mb higher level pathway groups
df_liverAssocScoresOverview$MB_Group <- unlist( lapply(X = df_liverAssocScoresOverview$MB_TermID, 
                                                       FUN = function(str_subsys) {
                                                         str_subsys <- unlist(strsplit(x = str_subsys, fixed = T, split = ","))
                                                         return( paste0(collapse = ",", sort( unique( df_MBpathwayGroups[df_MBpathwayGroups$MB_TermID %in% str_subsys,"MB_Group"] ))) )   
                                                       }) )
df_liverAssocScoresOverview$HOST_Group <- unlist( lapply(X = df_liverAssocScoresOverview$HOST_TermID, 
                                                         FUN = function(str_pwy) {
                                                           str_pwy <- unlist(strsplit(x = str_pwy, fixed = T, split = ","))
                                                           return( paste0(collapse = ",", sort( unique( df_hostPathwayGroups[df_hostPathwayGroups$HOST_TermID %in% str_pwy,"HOST_Group"] ))) )   
                                                         }) )
# reorder terms
df_liverAssocScoresOverview <- df_liverAssocScoresOverview[,c("HOST_TermID", "MB_TermID", "EnrichPadj", "SignedLogPadj", "HOST_BioProcess", "HOST_Group", "MB_SubSystem", "MB_Group",
                                                              "CorrFeatures", "MB_FeatureCount", "HOST_FeatureCount", "HOST_EnrichedFeatures", "MB_EnrichedFeatureIDs")]
df_liverAssocScoresOverview <- df_liverAssocScoresOverview[order(abs(df_liverAssocScoresOverview$EnrichPadj), decreasing = F),]

write.table(x = df_liverAssocScoresOverview,
            file = "enrichmentHostTranscriptsLiverToMicrobiomeReactions.tsv",
            append = F, quote = T, sep = "\t", row.names = F, col.names = T)

################

# Filter results for plotting:
table(df_overRepLiverToMBprocessesTop$p.adj <= 1e-4)
# FALSE  TRUE 
#   195    88 

## convert 3-column table to x-y-matrix
mtx_overRepLiverToMBproc <- with(df_overRepLiverToMBprocessesTop[df_overRepLiverToMBprocessesTop$p.adj <= 1e-4,],
                                 tapply(SignedLogP, list(GOID, SubsysID), sum))
dim(mtx_overRepLiverToMBproc)
# [1] 37 50
mtx_overRepLiverToMBproc[is.na(mtx_overRepLiverToMBproc)] <- 0

#### summarize processes that have the same description name:
### Row-Wise summary over GO-terms
# helper table with short Description for each GO-ID
df_tmp <- data.frame(GOID = rownames(mtx_overRepLiverToMBproc),
                     Description = str_simplifyDescriptions(df_GOdescriptions[rownames(mtx_overRepLiverToMBproc),"Description"], bol_stripRomanNums = F))
# sum up columns of GO-Ids with same description text
df_tmp2 <- matrix(data = unlist(tapply(X = df_tmp$GOID, INDEX = df_tmp$Description, 
                                       FUN = function(str_GOId) {
                                         if(length(str_GOId) == 1) return( c(str_GOId, mtx_overRepLiverToMBproc[str_GOId,]) )
                                         else return( c(str_GOId[1], colSums(mtx_overRepLiverToMBproc[str_GOId,])) )
                                       })), ncol = ncol(mtx_overRepLiverToMBproc)+1, byrow = T)
# fix row and col names and convert back to numeric
rownames( df_tmp2 ) <- df_tmp2[,1]; df_tmp2 <- df_tmp2[,-1]
class( df_tmp2 ) <- "numeric"
colnames(df_tmp2) <- colnames(mtx_overRepLiverToMBproc)
### Column-Wise summary over Subsystems:
# helper table with short Description for each SubSys:
df_tmp <- data.frame(SubsysID = colnames(mtx_overRepLiverToMBproc),
                     Description = str_simplifyDescriptions(df_SubSysdescriptions[colnames(mtx_overRepLiverToMBproc),"Subsystem"], bol_stripRomanNums = T))
# sum up columns of GO-Ids with same description text
df_tmp3 <- matrix(data = unlist(tapply(X = df_tmp$SubsysID, INDEX = df_tmp$Description, 
                                       FUN = function(str_GOId) {
                                         if(length(str_GOId) == 1) return( c(str_GOId, df_tmp2[,str_GOId]) )
                                         else return( c(str_GOId[1], rowSums(df_tmp2[,str_GOId])) )
                                       })), nrow = nrow(df_tmp2)+1, byrow = F)
# fix row and col names and convert back to numeric
colnames( df_tmp3 ) <- df_tmp3[1,]; df_tmp3 <- df_tmp3[-1,]
class( df_tmp3 ) <- "numeric"
rownames(df_tmp3) <- rownames(df_tmp2)
# save back
mtx_overRepLiverToMBproc <- df_tmp3; rm(df_tmp, df_tmp2, df_tmp3)
dim(mtx_overRepLiverToMBproc)
# [1] 37 45
######

## helper tables for home cooked clustering based on shared features (Genes, Reactions) in Process (GO, Subsystem)
# subsys geneID and ID
df_tmpSubSys <- df_overRepLiverToMBprocessesTop[,c("SubsysID", "Rxns", "SubSysBackground")]
colnames(df_tmpSubSys) <- c("ID", "geneID", "Count")
df_tmpSubSys$geneID <- gsub(pattern = ",", replacement = "/", x = df_tmpSubSys$geneID, fixed = T)
df_tmpSubSys <- df_tmpSubSys[!duplicated(df_tmpSubSys),]
rownames(df_tmpSubSys) <- df_tmpSubSys$ID
# GO geneID and ID
df_tmpGO <- df_overRepLiverToMBprocessesTop[,c("GOID", "Genes", "GOBackground")]
colnames(df_tmpGO) <- c("ID", "geneID", "Count")
df_tmpGO$geneID <- gsub(pattern = ",", replacement = "/", x = df_tmpGO$geneID, fixed = T)
df_tmpGO <- df_tmpGO[!duplicated(df_tmpGO),]
rownames(df_tmpGO) <- df_tmpGO$ID

## Select pathways to plot for heatmap
# Use all processes
mtx_tmp <- mtx_overRepLiverToMBproc
# drop empty features
mtx_tmp <- mtx_tmp[abs(rowSums(mtx_tmp)) > 0,
                   abs(colSums(mtx_tmp)) > 0, drop = F]
# # display only strong results and not single hit ones
int_clustSize = 2
while(sum(colSums(abs(mtx_tmp) > 0) < int_clustSize)+sum(rowSums(abs(mtx_tmp) > 0) < int_clustSize) > 0) {
  # print(sort(colSums(abs(mtx_tmp) > 0)))
  mtx_tmp <- mtx_tmp[,colSums(abs(mtx_tmp) > 0) >= int_clustSize]
  mtx_tmp <- mtx_tmp[rowSums(abs(mtx_tmp) > 0) >= int_clustSize,]
  # print(sort(rowSums(abs(mtx_tmp) > 0)))
}
# drop empty features
mtx_tmp <- mtx_tmp[abs(rowSums(mtx_tmp)) > 0,
                   abs(colSums(mtx_tmp)) > 0, drop = F]
## homecooked clustering
# distances based on shared features
mtx_rowDist <- mtx_mydist( df_GOresultTable = df_tmpGO[rownames(mtx_tmp),],
                           mtx_assocScores = mtx_tmp )
mtx_colDist <- mtx_mydist( df_GOresultTable = df_tmpSubSys[colnames(mtx_tmp),],
                           mtx_assocScores = mtx_tmp )
mtx_tmp <- mtx_tmp[hclust(as.dist(mtx_rowDist))$order,
                   hclust(as.dist(mtx_colDist))$order]
dim(mtx_tmp)
# [1] 17 16
max(abs(mtx_tmp))
# [1] 14.1


###############
## make a dot plot from it
# Invert row-order
mtx_tmp <- mtx_tmp[nrow(mtx_tmp):1,]
# manually abbreviate labels
str_colnames <- Hmisc::capitalize( str_simplifyDescriptions( df_SubSysdescriptions[colnames(mtx_tmp),"Subsystem"] ))
str_rownames <- Hmisc::capitalize( str_simplifyDescriptions( df_GOdescriptions[rownames(mtx_tmp),"Description"], bol_stripRomanNums = F))
str_rownames[1] <- "Mitochondrial Ca2+ transmem. transp."
str_rownames[14] <- "Neg. reg. of vascular muscle cell prolif."
str_rownames[15] <- "Detection of mechanical stimulus in pain"
#
pdf(file = "dots_liverFunctionalAssociationTerms.pdf", width = 400/72, height = 250/72, pointsize = 6, colormodel = "srgb")
par(mar = c(2,20,12,0)+.1, xpd = F, lheight=.75)
# Size of points, indicate strenght of host to microbiome associations
dbl_size <- abs(as.numeric(mtx_tmp))
dbl_size[dbl_size>0] <- scales::rescale(dbl_size[dbl_size>0], 
                                        to = c(1.3,4.0),  
                                        from = c(1e-12,50)) # max value encounter in our data is 50 for liver log10(p-val)
# color of points
lst_colors <- colorRampPalette(brewer.pal("RdBu", n = 11))(120)[120:1]
lst_colors <- lst_colors[-c(41:80)]
lst_colors <- lst_colors[ round( scales::rescale(as.numeric(mtx_tmp), to = c(1,80), from = c(-50,50)) ) ]
# empty plot
plot(y = NA,
     x = NA,
     xlim = c(1,ncol(mtx_tmp)),
     ylim = c(1,nrow(mtx_tmp)),
     cex = NA,
     axes = F, frame.plot = F,
     xlab = "", ylab = ""
)
# add grid
abline(v = 1:ncol(mtx_tmp),#-.5,
       h = 1:nrow(mtx_tmp),#-.5, 
       col = rgb(.7,.7,.7,.4), lwd = .7)
# add points
points(x = rep(1:ncol(mtx_tmp), each = nrow(mtx_tmp)),
       y = rep(1:nrow(mtx_tmp), times = ncol(mtx_tmp)),
       cex = dbl_size,
       pch = 21, lwd = .3,
       bg = c("#D6604D","#4393C3")[(as.numeric(mtx_tmp) > 0)+1]) #lst_colors) 
# MB
axis(side = 3, at = 1:ncol(mtx_tmp), line = 0.25,
     labels = rep("", ncol(mtx_tmp)),
     las = 2, cex.axis = 1.2)
text(x = 1:ncol(mtx_tmp), y = nrow(mtx_tmp)+1.8, xpd = T,
     labels = str_colnames,
     srt = 360-30, adj = c(1,0.5), cex = 1.2)
# Host
axis(side = 2,  at = 1:nrow(mtx_tmp), tick = T, line = 0, 
     labels = str_rownames, 
     las = 1, cex.axis = 1.2)
dev.off()
#
########     DONE     #########

