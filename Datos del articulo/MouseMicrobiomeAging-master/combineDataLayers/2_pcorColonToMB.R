#############################################
###
###   Run host microbiome correlation on Cluster
###
#############################################


#################
## Load functions and resources
source("../customFunctions.R")


#############################
# Load input data
df_vstMMrxnFiltered <- readRDS("df_MMrxnFiltered.rds")
df_colonRNAvstFiltered <- readRDS("df_colonRNAnormFiltered.rds")
df_metaData <- readRDS("df_colonMetaData.rds")


#############################
## partial correlation where age is taken care of
df_pcorRxn2ColonRNA <- df_pcorAll2All(df_table1 = df_vstMMrxnFiltered[, df_metaData$SampleName],
                                      df_table2 = df_colonRNAvstFiltered[, df_metaData$SampleName], 
                                      var_pcorZ = df_metaData$Age, 
                                      str_labels = c("Rxn","Gene"), 
                                      int_cores = 32)
# [1] "Starting partial correlation of 1500 Genes to 1500 Rxns within 52 samples! Wed Jan 26 11:05:42 2022"
# 
saveRDS(df_pcorRxn2ColonRNA, "df_pcorRxn2ColonRNA.rds")


############################
## End script
quit(save = "no")

