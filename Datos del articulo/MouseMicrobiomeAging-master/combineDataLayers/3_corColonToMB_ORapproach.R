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


##########################
# ## Generate Gene-Ontology parent process table - once
# #  requires that OverRepresentation scripts have been run for each organ
# df_hostPathwayGroups <- df_getLevel2GO( unique( unlist(strsplit(x = c(df_colonAssocScoresOverview$HOST_TermID,
#                                                                       df_brainAssocScoresOverview$HOST_TermID,
#                                                                       df_liverAssocScoresOverview$HOST_TermID), split = ",", fixed = T))) )
# df_hostPathwayGroups[df_hostPathwayGroups$QueryTerm == "",]
# # GO:0002376
# # GO:0008152
# # GO:0140972 
# df_hostPathwayGroups[df_hostPathwayGroups$QueryID == "GO:0002376",] <- c("GO:0002376","GO:0002376","immune system process","immune system process")
# df_hostPathwayGroups[df_hostPathwayGroups$QueryID == "GO:0008152",] <- c("GO:0008152","GO:0008152","metabolic process","metabolic process")
# df_hostPathwayGroups[df_hostPathwayGroups$QueryID == "GO:0140972",] <- c("GO:0140972","GO:0009987,GO:0065007","negative regulation of AIM2 inflammasome complex assembly",
#                                                                          "cellular process,biological regulation")
# write.table(x = df_hostPathwayGroups, 
#             file = "../databases/HostPathwayGroups.tsv", 
#             row.names = F, col.names = T, sep = "\t", quote = T)


##########################
## organize pathways into higher level groups (host and mb)
# Microbiome
df_MBpathwayGroups <- read.table(file = "../databases/MBpathwayGroups.tsv", header = T, sep = "\t", 
                                 stringsAsFactors = F, comment.char = "", as.is = T, quote = "\"")
df_MBpathwayGroups[,3] <- NULL
df_MBpathwayGroups[,5] <- NULL
colnames(df_MBpathwayGroups) <- c("MB_TermID", "MB_SubSystem", "MB_Group")
str_MBGroupTerms <- sort(unique(df_MBpathwayGroups$MB_Group))

# Host
df_hostPathwayGroups <- read.table(file = "../databases/HostPathwayGroups.tsv", header = T, sep = "\t", stringsAsFactors = F)
df_hostPathwayGroups <- df_hostPathwayGroups[,c(1,4)]
colnames(df_hostPathwayGroups) <- c("HOST_TermID","HOST_Group")
str_hostGroupTerms <- sort(unique(unlist(strsplit(split = ",", fixed = T, x = df_hostPathwayGroups$HOST_Group))))


##########################
## prepared subsystem annotation table
df_SubSysdescriptions <- unique( df_rxn2subsys[,c("SubsysID","Subsystem")] )
rownames(df_SubSysdescriptions) <- df_SubSysdescriptions$SubsysID


################################
## Aging regulated host genes
# Variance stabilized RNA-Seq reads informed with Age as integer
obj_DifAbundColonAllAges <- readRDS(file = "../mouseTranscriptome/obj_DifAbundColonAllAges.rds")
df_colonAgeDependance <- as.data.frame( results(obj_DifAbundColonAllAges, name = "age", alpha = 0.05) )
rm(obj_DifAbundColonAllAges)
df_colonAgeDependance <- df_colonAgeDependance[df_colonAgeDependance$padj <= 0.05 & 
                                               !is.na(df_colonAgeDependance$padj),]
## add GeneSymbol, req. for GO-annotation
df_colonAgeDependance$GeneSymbol <- df_mm10Ensembl2Genesymbol[rownames(df_colonAgeDependance),"GeneSymbol"]
# some are missing - help out with biomart online tools
library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
df_tmp <- getBM( values = rownames(df_colonAgeDependance)[is.na(df_colonAgeDependance$GeneSymbol)],
                 filters = "ensembl_gene_id", 
                 attributes = c("ensembl_gene_id","mgi_symbol"),
                 mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl")) )
rownames(df_tmp) <- df_tmp$ensembl_gene_id
# add remaining GeneSymbols:
df_colonAgeDependance[is.na(df_colonAgeDependance$GeneSymbol),"GeneSymbol"] <- df_tmp[rownames(df_colonAgeDependance)[is.na(df_colonAgeDependance$GeneSymbol)],"mgi_symbol"]
rm(df_tmp)
# order by p-val
df_colonAgeDependance <- df_colonAgeDependance[order(df_colonAgeDependance$pvalue, decreasing = F),]


################################
## partial Correlation results:
nrow(df_pcorRxn2ColonRNATop)
# 12732

## add GeneSymbol, req. for GO-annotation
df_pcorRxn2ColonRNATop$GeneSymbol <- df_mm10Ensembl2Genesymbol[df_pcorRxn2ColonRNATop$Gene,"GeneSymbol"]
# some are missing
unique(df_pcorRxn2ColonRNATop[is.na(df_pcorRxn2ColonRNATop$GeneSymbol),"Gene"])
# help out with biomart online tools
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
df_tmp <- getBM( values = unique(df_pcorRxn2ColonRNATop[is.na(df_pcorRxn2ColonRNATop$GeneSymbol),"Gene"]),
                 filters = "ensembl_gene_id", 
                 attributes = c("ensembl_gene_id","mgi_symbol"),
                 mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl")) )
rownames(df_tmp) <- df_tmp$ensembl_gene_id
# add remaining GeneSymbols:
df_pcorRxn2ColonRNATop[is.na(df_pcorRxn2ColonRNATop$GeneSymbol),"GeneSymbol"] <- df_tmp[df_pcorRxn2ColonRNATop[is.na(df_pcorRxn2ColonRNATop$GeneSymbol),"Gene"],"mgi_symbol"]
rm(df_tmp)


# Get GO bio processes for each gene in our correlation list
df_colonRNATopAllGOBioProcess <- df_mouseGOannotaionBioProc[df_mouseGOannotaionBioProc$GeneSymbol %in% unique(df_pcorRxn2ColonRNATop$GeneSymbol),c(1,2,6)]
df_colonRNATopAllGOBioProcess <- df_colonRNATopAllGOBioProcess[!duplicated(df_colonRNATopAllGOBioProcess),]
# concatenate genes into a single row per GO-term
df_tmp <- as.data.frame( tapply(X = df_colonRNATopAllGOBioProcess, INDEX = df_colonRNATopAllGOBioProcess$GOid, FUN = function(df_tmp) {
                                  return( paste0(df_tmp$GeneSymbol, collapse = ",") )
                                }) )
colnames(df_tmp) <- "Features"
df_colonRNATopAllGOBioProcess <- merge(x = df_tmp, y = df_colonRNATopAllGOBioProcess[,c(1,3)],
                                       by.x = "row.names", by.y = "GOid", all.x = T, all.y = F)
colnames(df_colonRNATopAllGOBioProcess) <- c("ID", "Features", "Description")
df_colonRNATopAllGOBioProcess <- df_colonRNATopAllGOBioProcess[!duplicated(df_colonRNATopAllGOBioProcess),]
rownames(df_colonRNATopAllGOBioProcess) <- df_colonRNATopAllGOBioProcess$ID
# feature count per GO-term (only features found in our correlation-data)
df_colonRNATopAllGOBioProcess$Count <- stringi::stri_count(str = df_colonRNATopAllGOBioProcess$Features, fixed = ",")+1
# keep only processes with at least 3 features:
df_colonRNATopAllGOBioProcess <- df_colonRNATopAllGOBioProcess[df_colonRNATopAllGOBioProcess$Count >= 3,]
# Drop certain terms like anything with "embryo*" in it for example
df_colonRNATopAllGOBioProcess <- df_colonRNATopAllGOBioProcess[grep(pattern = "embryo*", x = df_colonRNATopAllGOBioProcess$Description, ignore.case = T, invert = T),]
# add background size of term
df_colonRNATopAllGOBioProcess$TermSize <- tapply(X = df_mouseGOannotaionBioProc$GOid, INDEX = df_mouseGOannotaionBioProc$GOid, FUN = length)[df_colonRNATopAllGOBioProcess$ID]
# calculate overrepresentation of gene sets
int_nTotalInputGenes <- length(unique(unlist(strsplit(x = df_colonRNATopAllGOBioProcess$Features, split = ",", fixed = T))))
int_nTotalBackgroundGenes <- length(unique(df_mouseGOannotaionBioProc$GeneSymbol))
df_colonRNATopAllGOBioProcess$EnrichPval <- unlist( lapply(X = 1:nrow(df_colonRNATopAllGOBioProcess),
                                                   FUN = function(int_i) {
                                                      return( phyper(df_colonRNATopAllGOBioProcess[int_i,"Count"]-1, 
                                                              int_nTotalInputGenes, 
                                                              int_nTotalBackgroundGenes-int_nTotalInputGenes, 
                                                              df_colonRNATopAllGOBioProcess[int_i,"TermSize"],
                                                              lower.tail= FALSE) )
                                                    }))
df_colonRNATopAllGOBioProcess$EnrichPadj <- p.adjust(p = df_colonRNATopAllGOBioProcess$EnrichPval, method = "BH")
# NOTE: It makes sense to filter processes on enrichment p-val.,
#  since we don't have a defined pre-selection as we would have with metabolic modelling (MB side, or metamodel)
table(df_colonRNATopAllGOBioProcess$EnrichPadj <= 0.05)
# FALSE  TRUE 
#  1395   400 
df_colonRNATopAllGOBioProcess <- df_colonRNATopAllGOBioProcess[df_colonRNATopAllGOBioProcess$EnrichPadj <= 0.05,]


## Get subsystems for each rxn in our correlation list
df_colonRNATopAllSubSys <- df_rxn2subsys[df_rxn2subsys$ReactionID %in% unique(df_pcorRxn2ColonRNATop$Rxn),c(1:3)]
df_colonRNATopAllSubSys <- df_colonRNATopAllSubSys[!duplicated(df_colonRNATopAllSubSys),]
# concatenate rxns into a single row per subsys
df_tmp <- as.data.frame( tapply(X = df_colonRNATopAllSubSys, INDEX = df_colonRNATopAllSubSys$SubsysID, FUN = function(df_tmp) {
                                  return( paste0(df_tmp$ReactionID, collapse = ",") )
                                }) )
colnames(df_tmp) <- "Features"
df_colonRNATopAllSubSys <- merge(x = df_tmp, y = df_colonRNATopAllSubSys[,c(1,3)],
                                       by.x = "row.names", by.y = "SubsysID", all.x = T, all.y = F)
colnames(df_colonRNATopAllSubSys) <- c("ID", "Features", "Description")
df_colonRNATopAllSubSys <- df_colonRNATopAllSubSys[!duplicated(df_colonRNATopAllSubSys),]
rownames(df_colonRNATopAllSubSys) <- df_colonRNATopAllSubSys$ID
# feature count per GO-term (only features found in our correlation-data)
df_colonRNATopAllSubSys$Count <- stringi::stri_count(str = df_colonRNATopAllSubSys$Features, fixed = ",")+1
# keep only processes with at least 3 features:
df_colonRNATopAllSubSys <- df_colonRNATopAllSubSys[df_colonRNATopAllSubSys$Count >= 3,]
# NOTE: no enrichment filter, because subsystems are already pre-selected by the model reconstruction process


### overrepresentation tests for each process-pair (GO and Subsys) in parallel
int_nTotalStrongCorrs <- nrow(df_pcorRxn2ColonRNATop)
obj_clusters <- makeCluster(7, type="SOCK")
clusterExport(cl = obj_clusters, varlist = c("df_colonRNATopAllGOBioProcess","df_colonRNATopAllSubSys","df_pcorRxn2ColonRNATop","int_nTotalStrongCorrs"), envir = environment())
date()
# [1] "Tue Aug  8 10:34:10 2023"
df_overRepColonToMBprocesses <- as.data.frame(parLapply(cl = obj_clusters, 
                                                        X = 1:nrow(df_colonRNATopAllGOBioProcess), 
                                                        fun = function(int_i) {
  return( lapply(X = 1:nrow(df_colonRNATopAllSubSys), FUN = function(int_j) {
    str_geneList <- unlist(strsplit(x = df_colonRNATopAllGOBioProcess[int_i,"Features"], split = ",", fixed = T))
    str_rxnList <- unlist(strsplit(x = df_colonRNATopAllSubSys[int_j,"Features"], split = ",", fixed = T))
    
    # sign. correlations overlapping in this GO processes and this SubSystem
    int_countOverlapCors <- sum( (df_pcorRxn2ColonRNATop$GeneSymbol %in% str_geneList) &
                                 (df_pcorRxn2ColonRNATop$Rxn %in% str_rxnList) )
    # total sign. correlation found
    int_countGOBackgroundCors <- sum( df_pcorRxn2ColonRNATop$GeneSymbol %in% str_geneList )
    int_countSubSysBackgroundCors <- sum( df_pcorRxn2ColonRNATop$Rxn %in% str_rxnList )
    
    if(int_countOverlapCors > 0) {
      # print(paste(int_i, int_j))
      # print(c(length(dbl_rho) - 1, 
      #         df_colonRNATopAllSubSys[1,"Count"],
      #         int_nTotalStrongCorrs - df_colonRNATopAllSubSys[1,"Count"],
      #         df_colonRNATopAllGOBioProcess[int_i,"Count"]))
      
      # correlation values for overlapping MB and Host features
      dbl_rho <- mean( df_pcorRxn2ColonRNATop[df_pcorRxn2ColonRNATop$GeneSymbol %in% str_geneList &
                                              df_pcorRxn2ColonRNATop$Rxn %in% str_rxnList,"Spearman_rho"] )
      
      # overrepresentation p-value (more strong correlations in this certain subset (combination of 1 GO with 1 Subsys) than expected by chance)
      dbl_pOR <-  phyper(int_countOverlapCors - 1, 
                         int_countSubSysBackgroundCors,
                         int_nTotalStrongCorrs - int_countSubSysBackgroundCors,
                         int_countGOBackgroundCors,
                         lower.tail= F,
                         log.p = T)
      
      # return results
      return( c(df_colonRNATopAllSubSys[int_j,"ID"], # SubsysID
                df_colonRNATopAllGOBioProcess[int_i,"ID"],  # GOID
                dbl_pOR, # overrepresentation log(p-val)
                dbl_rho, # mean correlation value
                int_countOverlapCors,  # overlap size
                df_colonRNATopAllSubSys[int_j,"Count"], # group 1 background size
                df_colonRNATopAllGOBioProcess[int_i,"Count"], # group 2 background size
                int_countSubSysBackgroundCors, # total sign. cors of rxns of this Subsys
                int_countGOBackgroundCors # total sign. cors of features of this GO
                ) )
    }
    else return(c(NA,NA,NA,NA,NA,NA,NA,NA,NA))
  }) )
}) )
date()
# [1] "Tue Aug  8 10:35:37 2023" 
# clean up
stopCluster(obj_clusters)
# convert to data frame and clean up
df_overRepColonToMBprocesses <- as.data.frame(t(df_overRepColonToMBprocesses))
df_overRepColonToMBprocesses <- df_overRepColonToMBprocesses[!is.na(df_overRepColonToMBprocesses$V3),]
rownames(df_overRepColonToMBprocesses) <- NULL
colnames(df_overRepColonToMBprocesses) <- c("SubsysID","GOID", "logPval", "MeanRho", "Overlap", "SubSysBackground", "GOBackground","SubSysSignCors","GOSignCors")
# make numeric columns actually numeric
for(int_i in 3:9) {
  df_overRepColonToMBprocesses[,int_i] <- as.numeric(df_overRepColonToMBprocesses[,int_i])
}
df_overRepColonToMBprocesses$p.adj <- p.adjust(exp(df_overRepColonToMBprocesses$logPval), method = "BH")
df_overRepColonToMBprocesses <- df_overRepColonToMBprocesses[order(df_overRepColonToMBprocesses$logPval, decreasing = F),]
# calculate association score:
df_overRepColonToMBprocesses$Score <- df_overRepColonToMBprocesses$Overlap / (df_overRepColonToMBprocesses$SubSysBackground * df_overRepColonToMBprocesses$GOBackground)
# Store results
saveRDS(df_overRepColonToMBprocesses,"df_overRepColonToMBprocesses.rds")

table(df_overRepColonToMBprocesses$p.adj <= 0.01)
# FALSE  TRUE 
# 22383  2470 

# add GO infos
df_tmp <- merge(x = df_overRepColonToMBprocesses, y = df_colonRNATopAllGOBioProcess[,1:3],
                by.x = "GOID", by.y = "ID", all.x = T, all.y = F)
colnames(df_tmp) <- c("GOID", "SubsysID", "logPval", "MeanRho", "Overlap", "SubSysBackground", "GOBackground", "SubSysSignCors",
                      "GOSignCors", "p.adj", "Score", "Genes", "GOTerm")
df_overRepColonToMBprocesses <- merge(x = df_tmp, y = df_colonRNATopAllSubSys[,1:3],
                                      by.x = "SubsysID", by.y = "ID", all.x = T, all.y = F)
colnames(df_overRepColonToMBprocesses) <- c("SubsysID", "GOID", "logPval", "MeanRho", "Overlap", "SubSysBackground", "GOBackground", "SubSysSignCors",
                                            "GOSignCors", "p.adj", "Score", "Genes", "GOTerm", "Rxns", "SubSystem")

# filter by pval
df_overRepColonToMBprocessesTop <- df_overRepColonToMBprocesses[df_overRepColonToMBprocesses$p.adj <= 0.01,]

## check for duplicate entries
# simplified process names
df_overRepColonToMBprocessesTop$GOSimple <- str_simplifyDescriptions( df_overRepColonToMBprocessesTop$GOTerm, bol_stripRomanNums = F )
df_overRepColonToMBprocessesTop$SubSysSimple <- str_simplifyDescriptions( df_overRepColonToMBprocessesTop$SubSystem, bol_stripRomanNums = T )
# exact duplicates (usually due to redudandant subsystem names)
table(duplicated(df_overRepColonToMBprocessesTop[,c(2:14,17)]))
df_overRepColonToMBprocessesTop <- df_overRepColonToMBprocessesTop[!duplicated(df_overRepColonToMBprocessesTop[,c(2:14,17)]),]
# dropped 195

# same simple names, choose better logP (min)
table(duplicated(df_overRepColonToMBprocessesTop[,c("GOSimple","SubSysSimple")]))
df_overRepColonToMBprocessesTop <- df_overRepColonToMBprocessesTop[order(df_overRepColonToMBprocessesTop$logPval, decreasing = F),]
df_overRepColonToMBprocessesTop <- df_overRepColonToMBprocessesTop[!duplicated(df_overRepColonToMBprocessesTop[,c("GOSimple","SubSysSimple")]),]
# dropped 414

# similar correlation results but (slightly) different subsys name 
table(duplicated(df_overRepColonToMBprocessesTop[,c(2:14)]))
df_overRepColonToMBprocessesTop$DupID <-  as.numeric( as.factor(unlist( apply(X = df_overRepColonToMBprocessesTop[,c(2:14)], MARGIN = 1, function(x) {
                                                                                return( paste(collapse = "|", x) )
                                                                              }))))
#  merge those into one single column
df_tmp <- as.data.frame( matrix( data = unlist( tapply(X = 1:nrow(df_overRepColonToMBprocessesTop), INDEX = df_overRepColonToMBprocessesTop$DupID, 
                               FUN = function(int_i) {
                                 if(length(int_i) == 1) return(c(df_overRepColonToMBprocessesTop[int_i,,drop = T],df_overRepColonToMBprocessesTop[int_i,"SubsysID"]))
                                 else {
                                    # select one representative (shortest label)
                                    df_tmp <- df_overRepColonToMBprocessesTop[int_i,]
                                    df_tmp <- df_tmp[order(stringi::stri_length(df_tmp$SubSysSimple), decreasing = F),]
                                    # concatenate subsys Descriptions and Ids
                                    df_tmp[1,"SubSysIDs"] <- paste(collapse = ",", df_tmp$SubsysID)
                                    df_tmp[1,"SubSystem"] <- paste(collapse = ",", df_tmp$SubSystem)
                                    
                                    return(c(df_tmp[1,,drop=T]))
                                 }
                              }, simplify = T)), ncol = 19, byrow = T) )
colnames(df_tmp) <- c(colnames(df_overRepColonToMBprocessesTop),"SubSysIDs")
df_tmp$DupID <- NULL
df_overRepColonToMBprocessesTop <- df_tmp; rm(df_tmp)
# make numeric columns actually numeric
for(int_i in 3:11) {
  df_overRepColonToMBprocessesTop[,int_i] <- as.numeric(df_overRepColonToMBprocessesTop[,int_i])
}
# merged 484 subsys

# duplicated GO terms?
table(duplicated(df_overRepColonToMBprocessesTop[,c(1,3:14,16)]))
# -none-


df_overRepColonToMBprocessesTop$SignedLogP <- log10(df_overRepColonToMBprocessesTop$p.adj) * (((df_overRepColonToMBprocessesTop$MeanRho < 0)*2)-1)


################################
## Export to table - Colon Correlation results table S1
# View(df_overRepColonToMBprocessesTop[,c("GOID", "SubSysIDs", "p.adj", "SignedLogP", "GOTerm", "SubSystem", "Overlap", "SubSysBackground", "GOBackground", "Genes", "Rxns" )])
# Export table in a tabular format for people to explore as excel-table
df_colonAssocScoresOverview <- df_overRepColonToMBprocessesTop[,c("GOID", "SubSysIDs", "p.adj", "SignedLogP", "GOTerm", "SubSystem", "Overlap", "SubSysBackground", "GOBackground", "Genes", "Rxns" )]
colnames(df_colonAssocScoresOverview) <- c("HOST_TermID", "MB_TermID", "EnrichPadj", "SignedLogPadj", "HOST_BioProcess", "MB_SubSystem", 
                                           "CorrFeatures", "MB_FeatureCount", "HOST_FeatureCount", "HOST_EnrichedFeatures", "MB_EnrichedFeatureIDs")

### ADD higher level group
# host and mb higher level pathway groups
df_colonAssocScoresOverview$MB_Group <- unlist( lapply(X = df_colonAssocScoresOverview$MB_TermID, 
                                               FUN = function(str_subsys) {
                                                 str_subsys <- unlist(strsplit(x = str_subsys, fixed = T, split = ","))
                                                 return( paste0(collapse = ",", sort( unique( df_MBpathwayGroups[df_MBpathwayGroups$MB_TermID %in% str_subsys,"MB_Group"] ))) )   
                                              }) )
df_colonAssocScoresOverview$HOST_Group <- unlist( lapply(X = df_colonAssocScoresOverview$HOST_TermID, 
                                               FUN = function(str_pwy) {
                                                 str_pwy <- unlist(strsplit(x = str_pwy, fixed = T, split = ","))
                                                 return( paste0(collapse = ",", sort( unique( df_hostPathwayGroups[df_hostPathwayGroups$HOST_TermID %in% str_pwy,"HOST_Group"] ))) )   
                                               }) )
# reorder terms
df_colonAssocScoresOverview <- df_colonAssocScoresOverview[,c("HOST_TermID", "MB_TermID", "EnrichPadj", "SignedLogPadj", "HOST_BioProcess", "HOST_Group", "MB_SubSystem", "MB_Group",
                                                              "CorrFeatures", "MB_FeatureCount", "HOST_FeatureCount", "HOST_EnrichedFeatures", "MB_EnrichedFeatureIDs")]
df_colonAssocScoresOverview <- df_colonAssocScoresOverview[order(abs(df_colonAssocScoresOverview$EnrichPadj), decreasing = F),]

write.table(x = df_colonAssocScoresOverview,
            file = "enrichmentHostTranscriptsColonToMicrobiomeReactions.tsv",
            append = F, quote = T, sep = "\t", row.names = F, col.names = T)


################
# Filter results for plotting:
table(df_overRepColonToMBprocessesTop$p.adj <= 10^-10)
# FALSE  TRUE 
# 1252   125 

## convert 3-column table to x-y-matrix
mtx_overRepColonToMBproc <- with(df_overRepColonToMBprocessesTop[df_overRepColonToMBprocessesTop$p.adj <= 1e-10,],
                                 tapply(SignedLogP, list(GOID, SubsysID), sum))
dim(mtx_overRepColonToMBproc)
# [1] 36 38
mtx_overRepColonToMBproc[is.na(mtx_overRepColonToMBproc)] <- 0

#### summarize processes that have the same description name:
### Row-Wise summary over GO-terms
# helper table with short Description for each GO-ID
df_tmp <- data.frame(GOID = rownames(mtx_overRepColonToMBproc),
                     Description = str_simplifyDescriptions(df_GOdescriptions[rownames(mtx_overRepColonToMBproc),"Description"], bol_stripRomanNums = F))
# sum up columns of GO-Ids with same description text
df_tmp2 <- matrix(data = unlist(tapply(X = df_tmp$GOID, INDEX = df_tmp$Description, 
                                       FUN = function(str_GOId) {
                                         if(length(str_GOId) == 1) return( c(str_GOId, mtx_overRepColonToMBproc[str_GOId,]) )
                                         else return( c(str_GOId[1], colSums(mtx_overRepColonToMBproc[str_GOId,])) )
                                      })), ncol = ncol(mtx_overRepColonToMBproc)+1, byrow = T)
# fix row and col names and convert back to numeric
rownames( df_tmp2 ) <- df_tmp2[,1]; df_tmp2 <- df_tmp2[,-1]
class( df_tmp2 ) <- "numeric"
colnames(df_tmp2) <- colnames(mtx_overRepColonToMBproc)
### Column-Wise summary over Subsystems:
# helper table with short Description for each SubSys:
df_tmp <- data.frame(SubsysID = colnames(mtx_overRepColonToMBproc),
                     Description = str_simplifyDescriptions(df_SubSysdescriptions[colnames(mtx_overRepColonToMBproc),"Subsystem"], bol_stripRomanNums = T))
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
mtx_overRepColonToMBproc <- df_tmp3; rm(df_tmp, df_tmp2, df_tmp3)
dim(mtx_overRepColonToMBproc)
# [1] 36 37
######

## helper tables for home cooked clustering based on shared features (Genes, Reactions) in Process (GO, Subsystem)
# subsys geneID and ID
df_tmpSubSys <- df_overRepColonToMBprocessesTop[,c("SubsysID", "Rxns", "SubSysBackground")]
colnames(df_tmpSubSys) <- c("ID", "geneID", "Count")
df_tmpSubSys$geneID <- gsub(pattern = ",", replacement = "/", x = df_tmpSubSys$geneID, fixed = T)
df_tmpSubSys <- df_tmpSubSys[!duplicated(df_tmpSubSys),]
rownames(df_tmpSubSys) <- df_tmpSubSys$ID
# GO geneID and ID
df_tmpGO <- df_overRepColonToMBprocessesTop[,c("GOID", "Genes", "GOBackground")]
colnames(df_tmpGO) <- c("ID", "geneID", "Count")
df_tmpGO$geneID <- gsub(pattern = ",", replacement = "/", x = df_tmpGO$geneID, fixed = T)
df_tmpGO <- df_tmpGO[!duplicated(df_tmpGO),]
rownames(df_tmpGO) <- df_tmpGO$ID

## Select pathways to plot for heatmap
# Use all processes
mtx_tmp <- mtx_overRepColonToMBproc
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
# [1] 26 23
max(abs(mtx_tmp))
# [1] 49.48309


###############
## make a dot plot from it
# Invert row-order
mtx_tmp <- mtx_tmp[nrow(mtx_tmp):1,]
# manually abbreviate labels
str_colnames <- Hmisc::capitalize( str_simplifyDescriptions( df_SubSysdescriptions[colnames(mtx_tmp),"Subsystem"] ))
str_colnames[21] <- "Sphingosine metab."
str_rownames <- Hmisc::capitalize( str_simplifyDescriptions( df_GOdescriptions[rownames(mtx_tmp),"Description"], bol_stripRomanNums = F))
str_rownames[25] <- "mRNA proc."
str_rownames[24] <- "mRNA splicing, via spliceosome" 
str_rownames[20] <- "DNA damage resp. via p53 mediator"
str_rownames[16] <- "Chromatin remodeling"
str_rownames[15] <- "Mitochondrial complex I assembly"
str_rownames[13] <- "Mitochondrial ATP synth. by H+"
#
pdf(file = "dots_colonFunctionalAssociationTerms.pdf", width = 450/72, height = 320/72, pointsize = 6, colormodel = "srgb")
par(mar = c(2,18,12,0)+.1, xpd = F, lheight=.75)
# Size of points, indicate strenght of host to microbiome associations
dbl_size <- abs(as.numeric(mtx_tmp))
dbl_size[dbl_size>0] <- scales::rescale(dbl_size[dbl_size>0], 
                                        to = c(1.3,4.0),    
                                        from = c(1e-12,50)) # max value encounter in our data is 50 for colon log10(p-val)
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
       bg = c("#4393C3", "#D6604D")[(as.numeric(mtx_tmp) > 0)+1]) #lst_colors)
# MB
axis(side = 3, at = 1:ncol(mtx_tmp), line = 0,
     labels = rep("", ncol(mtx_tmp)),
     las = 2, cex.axis = 1.2)
text(x = 1:ncol(mtx_tmp), y = nrow(mtx_tmp)+2, xpd = T,
     labels = str_colnames,
     srt = 360-30, adj = c(1,0.5), cex = 1.2) #30
# Host
axis(side = 2,  at = 1:nrow(mtx_tmp), tick = T, line = 0, 
     labels = str_rownames, #rep("", nrow(mtx_tmp)),
     las = 1, cex.axis = 1.2)
# text(y = 1:nrow(mtx_tmp), x = -1, xpd = T,
#      labels = str_rownames,
#      srt = 360-30, adj = c(1, 0.5), cex = 1.2)
dev.off()
#
## Legend
pdf(file = "dots_legend.pdf", width = 200/72, height = 200/72, pointsize = 6, colormodel = "srgb")
par(mar = c(2,2,2,2)+.1, xpd = T, lheight=.75)
# lst_colors <- colorRampPalette(brewer.pal("RdBu", n = 11))(120)[120:1]
# lst_colors <- lst_colors[-c(51:70)]
# empty plot
plot(y = NA,
     x = NA,
     xlim = c(1,ncol(mtx_tmp)),
     ylim = c(1,nrow(mtx_tmp)),
     cex = NA,
     axes = F, frame.plot = F,
     xlab = "", ylab = ""
)
legend("center", legend = c(" -50"," -25"," -2","  50","  25","  2"), 
       horiz = F, cex = 1.2, ncol = 2,
       title = "Signed log10( FDR )", title.font = 2,
       pt.bg = rep(c("#4393C3", "#D6604D"), each = 3), #lst_colors[c(1,25,48,52,75,100)], 
       pt.cex = scales::rescale(c(50,25,2,50,25,2), 
                                to = c(1.3,4.0), 
                                from = c(1e-12,50)),
       y.intersp = 1.5, x.intersp = .8, pch = 21, pt.lwd = .3, bty = "n")
dev.off()
########     DONE     #########

