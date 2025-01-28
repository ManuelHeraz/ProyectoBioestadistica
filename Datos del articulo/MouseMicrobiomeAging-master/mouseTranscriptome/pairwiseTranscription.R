###########################################################
###
###   Pairwise DESeq of RNA by Age, stratified by tissue
###
###                   L.B. 01-2022
###
###########################################################

################
##
library(DESeq2)
source("customFunctions.R")
setwd("mouseTranscriptome/")


########
# Load metadata and counn data
df_transcriptMetaData <- readRDS("df_transcriptMetaData.rds")
mtx_CountsColonAllAges <- readRDS("mtx_CountsColonAllAges.rds")
mtx_CountsLiverAllAges <- readRDS("mtx_CountsLiverAllAges.rds")
mtx_CountsBrainAllAges <- readRDS("mtx_CountsBrainAllAges.rds")


####################
## Colon:
# select sample-metadata for this organ only
df_transcriptMetaDataColon <- df_transcriptMetaData[df_transcriptMetaData$Tissue == "Colon",]
rownames(df_transcriptMetaDataColon) <- df_transcriptMetaDataColon$SampleName

# visual QC
void_MDSplot(df_sampleMatrix = mtx_CountsColonAllAges[,df_transcriptMetaDataColon$id], 
             df_metaData = df_transcriptMetaDataColon, bol_labels = T, var_colorBy = "Age", 
             var_shapeBy = NULL, bol_continousColor = F)

# DESeq-object from existing count matrix
obj_DESeqColonAllAges <- DESeqDataSetFromMatrix(countData = mtx_CountsColonAllAges[,df_transcriptMetaDataColon$SampleName], 
                                                colData = df_transcriptMetaDataColon,
                                                design = ~ Age)
obj_DifAbundColonAllAges <- DESeq(obj_DESeqColonAllAges)

## automatic testing of all age-groups vs all age-groups
for(int_i in 1:(length(sort(unique(df_transcriptMetaDataColon$Age)))-1)) {
  int_group1 <- sort(unique(df_transcriptMetaDataColon$Age))[int_i]
  for(int_j in (int_i+1):(length(sort(unique(df_transcriptMetaDataColon$Age))))) {
    int_group2 <- sort(unique(df_transcriptMetaDataColon$Age))[int_j]
    #
    print(paste(int_group1, "vs.", int_group2))

    if(!is.null(obj_DifAbundColonAllAges)) {
      df_DifAbundColonAllAges <- results(obj_DifAbundColonAllAges, independentFiltering = T, alpha = 0.05,
                                         contrast = c("Age", int_group1, int_group2))
      print( summary(df_DifAbundColonAllAges))
      df_DifAbundColonAllAges <- as.data.frame(df_DifAbundColonAllAges)
      df_DifAbundColonAllAges <- df_DifAbundColonAllAges[order(df_DifAbundColonAllAges$pvalue, decreasing = F),]
      #
      write.table(x = df_DifAbundColonAllAges,
                  file = paste0("pairwiseDESeq/DESeqColon-Age-",int_group1,"-vs-",int_group2,"-positiveFCisUpIn-",int_group1,".txt"),
                  quote = F, sep = "\t", row.names = T, col.names = T)
    }
  }
}


####################
## Liver:
# select sample-metadata for this organ only
df_transcriptMetaDataLiver <- df_transcriptMetaData[df_transcriptMetaData$Tissue == "Liver",]
rownames(df_transcriptMetaDataLiver) <- df_transcriptMetaDataLiver$SampleName

# visual QC
void_MDSplot(df_sampleMatrix = mtx_CountsLiverAllAges[,df_transcriptMetaDataLiver$id], 
             df_metaData = df_transcriptMetaDataLiver, bol_labels = F, var_colorBy = "Age", 
             var_shapeBy = "Batch", bol_continousColor = F)

# DESeq-object from existing count matrix
obj_DESeqLiverAllAges <- DESeqDataSetFromMatrix(countData = mtx_CountsLiverAllAges[,df_transcriptMetaDataLiver$SampleName], 
                                                colData = df_transcriptMetaDataLiver,
                                                design = ~ Age)
obj_DifAbundLiverAllAges <- DESeq(obj_DESeqLiverAllAges)

## automatic testing of all age-groups vs all age-groups
for(int_i in 1:(length(sort(unique(df_transcriptMetaDataLiver$Age)))-1)) {
  int_group1 <- sort(unique(df_transcriptMetaDataLiver$Age))[int_i]
  for(int_j in (int_i+1):(length(sort(unique(df_transcriptMetaDataLiver$Age))))) {
    int_group2 <- sort(unique(df_transcriptMetaDataLiver$Age))[int_j]
    #
    print(paste(int_group1, "vs.", int_group2))
    
    if(!is.null(obj_DifAbundLiverAllAges)) {
      df_DifAbundLiverAllAges <- results(obj_DifAbundLiverAllAges, independentFiltering = T, alpha = 0.05,
                                         contrast = c("Age", int_group1, int_group2))
      print( summary(df_DifAbundLiverAllAges))
      df_DifAbundLiverAllAges <- as.data.frame(df_DifAbundLiverAllAges)
      df_DifAbundLiverAllAges <- df_DifAbundLiverAllAges[order(df_DifAbundLiverAllAges$pvalue, decreasing = F),]
      #
      write.table(x = df_DifAbundLiverAllAges,
                  file = paste0("pairwiseDESeq/DESeqLiver-Age-",int_group1,"-vs-",int_group2,"-positiveFCisUpIn-",int_group1,".txt"),
                  quote = F, sep = "\t", row.names = T, col.names = T)
    }
  }
}


####################
## Brain:
# select sample-metadata for this organ only
df_transcriptMetaDataBrain <- df_transcriptMetaData[df_transcriptMetaData$Tissue == "Brain",]
rownames(df_transcriptMetaDataBrain) <- df_transcriptMetaDataBrain$SampleName

# visual QC
void_MDSplot(df_sampleMatrix = mtx_CountsBrainAllAges[,df_transcriptMetaDataBrain$id], 
             df_metaData = df_transcriptMetaDataBrain, bol_labels = T, var_colorBy = "Age", 
             var_shapeBy = "Batch", bol_continousColor = F)

# DESeq-object from existing count matrix
obj_DESeqBrainAllAges <- DESeqDataSetFromMatrix(countData = mtx_CountsBrainAllAges[,df_transcriptMetaDataBrain$SampleName], 
                                                colData = df_transcriptMetaDataBrain,
                                                design = ~ Age)
obj_DifAbundBrainAllAges <- DESeq(obj_DESeqBrainAllAges)

## automatic testing of all age-groups vs all age-groups
for(int_i in 1:(length(sort(unique(df_transcriptMetaDataBrain$Age)))-1)) {
  int_group1 <- sort(unique(df_transcriptMetaDataBrain$Age))[int_i]
  for(int_j in (int_i+1):(length(sort(unique(df_transcriptMetaDataBrain$Age))))) {
    int_group2 <- sort(unique(df_transcriptMetaDataBrain$Age))[int_j]
    #
    print(paste(int_group1, "vs.", int_group2))
    
    if(!is.null(obj_DifAbundBrainAllAges)) {
      df_DifAbundBrainAllAges <- results(obj_DifAbundBrainAllAges, independentFiltering = T, alpha = 0.05,
                                         contrast = c("Age", int_group1, int_group2))
      print( summary(df_DifAbundBrainAllAges))
      df_DifAbundBrainAllAges <- as.data.frame(df_DifAbundBrainAllAges)
      df_DifAbundBrainAllAges <- df_DifAbundBrainAllAges[order(df_DifAbundBrainAllAges$pvalue, decreasing = F),]
      #
      write.table(x = df_DifAbundBrainAllAges,
                  file = paste0("pairwiseDESeq/DESeqBrain-Age-",int_group1,"-vs-",int_group2,"-positiveFCisUpIn-",int_group1,".txt"),
                  quote = F, sep = "\t", row.names = T, col.names = T)
    }
  }
}
