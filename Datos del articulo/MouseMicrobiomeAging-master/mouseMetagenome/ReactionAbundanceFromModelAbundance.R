##################################################
##
##  From MAG abundances to Reaction abundances 
##  
##  L. Best - 2023-01-01
##
##################################################

setwd("./mouseMetagenome/")


####################
## Load Models
# obj_metabolicModelsMouse <- readRDS("/home/lena/Dropbox/0_Work/0_Projects/18_EnrichedEnvironmentJena/metagenome/metamouse-20230510.RDS")
# 
# ##Get list of all reactions contained in the models
# lst_allRxns <- c()
# for(obj_model in obj_metabolicModelsMouse){
#   lst_allRxns <- union(lst_allRxns, obj_model@react_id)
# }
# 
# ##Build incidence matrix of reactions in each model
# mtx_rxnInModels <- matrix(0,length(lst_allRxns),length(obj_metabolicModelsMouse))
# rownames(mtx_rxnInModels) <- lst_allRxns
# colnames(mtx_rxnInModels) <- names(obj_metabolicModelsMouse)
# for(str_modelName in names(obj_metabolicModelsMouse)) {
#   obj_model <- obj_metabolicModelsMouse[[str_modelName]]
#   mtx_rxnInModels[obj_model@react_id,str_modelName] <- 1
# }
# 
# ##Normalize matrix such that the sum for each species is "1" -> species with larger genomes have 
# ##smaller contribution of individual reactions
# colSums(mtx_rxnInModels)
# mtx_rxnInModels <- prop.table(mtx_rxnInModels,2)
# colSums(mtx_rxnInModels)


####################
## Alternatively use FVA reactions
obj_metabolicModelsMouse <- readRDS("./fva_99_reactions_thr06.RDS")

##Get list of active reactions contained in the models
lst_allRxns <- c()
for(obj_model in obj_metabolicModelsMouse){
  lst_allRxns <- union(lst_allRxns, obj_model$active)
}

##Build incidence matrix of reactions in each model
mtx_rxnInModels <- matrix(0,length(lst_allRxns),length(obj_metabolicModelsMouse))
rownames(mtx_rxnInModels) <- lst_allRxns
colnames(mtx_rxnInModels) <- names(obj_metabolicModelsMouse)
for(str_modelName in names(obj_metabolicModelsMouse)) {
  obj_model <- obj_metabolicModelsMouse[[str_modelName]]
  mtx_rxnInModels[obj_model$active,str_modelName] <- 1
}

##Normalize matrix such that the sum for each species is "1" -> species with larger genomes have 
##smaller contribution of individual reactions
colSums(mtx_rxnInModels)
mtx_rxnInModels <- prop.table(mtx_rxnInModels,2)
colSums(mtx_rxnInModels)


###################
# load abundance data (jgi_summarize)
# df_hqBinAbundances <- read.csv(file = "../mouseMetagenome/hqBins/HQBinsScaffoldDepths.txt", 
#                                header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
# # drop "var" columns
# df_hqBinAbundances <- df_hqBinAbundances[,grep(".bam.var", colnames(df_hqBinAbundances), invert = T)]
# df_hqBinAbundances$contigLen <- NULL
# df_hqBinAbundances$totalAvgDepth <- NULL
# colnames(df_hqBinAbundances) <- sub(pattern = ".bam", replacement = "", fixed = T, 
#                                     x = sub(pattern = "X", replacement = "", x = colnames(df_hqBinAbundances), fixed = T))
# 
# ## Convert mean coverage depth per sample from contigs to bins (MAGs)
# df_MAGabundances <- df_hqBinAbundances
# ## load scaffolds 2 bin conversion table
# df_scaffoldsInBins <- read.csv(file = "../mouseMetagenome/hqBins/hqBinAnnotation/scaffoldsInBins", header = F, sep = "\t", stringsAsFactors = F)
# colnames(df_scaffoldsInBins) <- c("MAG","Scaffold")
# rownames(df_scaffoldsInBins) <- df_scaffoldsInBins$Scaffold
# ## Add MAG as new column
# df_MAGabundances$MAG <- df_scaffoldsInBins[rownames(df_MAGabundances),"MAG"]
# ## Sum up mean abundances over all contigs from each bin 
# # Getting confused here as to whether use sum or mean - but I think mean makes the most sense
# df_tmp <- apply(X = df_MAGabundances[,1:(ncol(df_MAGabundances)-1)], MARGIN = 2, FUN = function(col) {
#   return( tapply(X = as.numeric( col ), INDEX = df_MAGabundances$MAG, FUN = mean) )
# })
# # # validate output
# # # colSums(df_MAGabundances[df_MAGabundances$MAG == "concoct_12",1:(ncol(df_MAGabundances)-1)]) - df_tmp["concoct_12",colnames(df_MAGabundances)[1:(ncol(df_MAGabundances)-1)]]
# df_MAGabundances <- as.data.frame(df_tmp); rm(df_tmp)
# saveRDS(df_MAGabundances,"../mouseMetagenome/df_MAGabundances.rds")


#########################
## Assume that read counts are in "df_MAGabundances" (species are rows, samples in columns)
df_MAGabundances <- readRDS("./df_MAGabundances.rds")
lst_overlapMAGs <- intersect(rownames(df_MAGabundances), colnames(mtx_rxnInModels))
df_MAGabundances <- df_MAGabundances[lst_overlapMAGs,]

##Multiply reaction incidence matrix with the count data
df_rxnAbundancesMM <- t(t(df_MAGabundances[lst_overlapMAGs,]) %*% t(mtx_rxnInModels[,lst_overlapMAGs]))
# drop all zero rows
df_rxnAbundancesMM <- df_rxnAbundancesMM[rowSums(df_rxnAbundancesMM) != 0,]

##Normalize to account for differences in read depths
colSums(df_rxnAbundancesMM)
df_rxnAbundancesMM <- as.data.frame( prop.table(df_rxnAbundancesMM,2) )
colSums(df_rxnAbundancesMM)
# Reaction abundances are in df_rxnAbundancesMM
saveRDS(df_rxnAbundancesMM,"./df_rxnAbundancesMM.rds")


#########################
## Create an reaction annotation object from models
##Obtain list of reactions
rm(df_rxnAnnotation)
obj_metabolicModelsMouse <- readRDS("./metamouse-20230510.RDS")
for(obj_model in obj_metabolicModelsMouse){
  if(!exists("df_rxnAnnotation"))
    df_rxnAnnotation <- data.frame(RxnID = obj_model@react_id, RxnName = obj_model@react_name)
  else
    df_rxnAnnotation <- rbind(df_rxnAnnotation,
                              data.frame(RxnID = obj_model@react_id, RxnName = obj_model@react_name))
}
df_rxnAnnotation <- df_rxnAnnotation[!duplicated(df_rxnAnnotation$RxnID),]
rownames(df_rxnAnnotation) <- df_rxnAnnotation$RxnID
saveRDS(df_rxnAnnotation,"./df_rxnAnnotation.rds")

# clean up
rm(lst_overlapMAGs, mtx_rxnInModels, obj_metabolicModelsMouse, lst_allRxns, str_modelName, obj_model, df_scaffoldsInBins)
