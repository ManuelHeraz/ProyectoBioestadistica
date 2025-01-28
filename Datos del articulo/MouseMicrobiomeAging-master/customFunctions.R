############################################################
###
###   Custom functions and databases used in this project
###
############################################################

library(RColorBrewer)
library(DESeq2)
library(gplots)
library(clusterProfiler)
require(parallel)
require(car)
require(stats)
require("GO.db")

###
setwd(".")
###

#######################################
# function definitions

#####
## Partial correlation of all entries of one matrix to all entries of another matrix
## corrected for an optional variable var_pcorZ
df_pcorAll2All <- function(df_table1, df_table2, 
                           str_labels = c("Feature1", "Feature2"), 
                           var_pcorZ = NULL, int_cores = 7) {
  require(parallel)
  require(car)
  require(stats)
  # source("custom_functions.R")  <- contains pcor.test
  
  int_nFeatures1 <- nrow(df_table1)      #number of OTUs
  int_nSamples <- ncol(df_table1)  #number of Samples
  
  #make sure both tables have their sample names in the same order!
  df_table2 <- as.data.frame(df_table2)
  df_table1 <- as.data.frame(df_table1)
  # df_table2 <- df_table2[ , order(colnames(df_table2)) ]
  # df_table1 <- df_table1[ , order(colnames(df_table1)) ]
  
  #column names must be exactly matching -> break otherwise
  if( !isTRUE( all.equal( colnames(df_table2),  colnames(df_table1)) ) ) {
    print("Column names of the given data frames are not matching!")
    return(1)
  }
  
  # partial correlation?
  bol_pcor = F
  if(!is.null(var_pcorZ)) {
    bol_pcor = T
    print("Partial correlation enabled!")
  }
  
  ## retain rownames
  # add a column for the feature index (can later be used to get feature name back)
  # this column needs to be numeric as well in order to be succesfully passed to a function!
  fct_featureNames2 <- factor(x = rownames(df_table2))
  df_table2[,str_labels[2]] <- as.numeric( fct_featureNames2 )
  
  fct_featureNames1 <- factor(x = rownames(df_table1))
  df_table1[,str_labels[1]] <- as.numeric( fct_featureNames1 )
  
  # Internal correlation function
  df_one2all <- function(row_table2) {
    # extract feature id2
    int_featureID2 <- row_table2[int_nSamples + 1]
    # remove id and retain only numberic values from samples
    row_table2 <- row_table2[ 1:int_nSamples ]
    
    #
    lst_one2one <- function(row_table1) {
      # extract feature id2
      int_featureID1 <- row_table1[int_nSamples + 1]
      # remove id and retain only numberic values from samples
      row_table1 <- row_table1[ 1:int_nSamples ]
      
      # if too many NA-values (less than 3 not-NA left)
      if( sum(!is.na( row_table2 )) < 3 ||
          sum(!is.na( row_table1 )) < 3 ) {
        
        results <- results <- c(   int_featureID2,
                                   int_featureID1,
                                   NA, 
                                   NA  )  
      }
      else {
        #spearman correlation of OTU-abundance vs gene-expression
        if(bol_pcor) {
          dbl_spearmanRho <- pcor.test(x = row_table2, y = row_table1, 
                                       z = var_pcorZ,
                                       use = "mat",
                                       method = "spearman")
                                       # method = "kendall")
        }
        else {
          dbl_spearmanRho <- cor.test(x = row_table2, y = row_table1, 
                                      method = "spearman", alternative = "two.sided", exact = T, continuity = F)
        }
        results <- c( int_featureID2,
                      int_featureID1,
                      dbl_spearmanRho$estimate,
                      dbl_spearmanRho$p.value
        )
        # print(results)
      }
      
      return( results )
    }
    
    df_one2all <-  as.data.frame( t( apply(X = df_table1, MARGIN = 1, FUN = lst_one2one) ))
    colnames(df_one2all) <- c(str_labels[2],str_labels[1],"Spearman_rho","p.val") 
    
    # View(df_one2all)
    
    return( df_one2all ) 
  }
  
  # correlate each gene with each otu in parallel ...
  obj_clusters <- makeCluster(int_cores, type="SOCK")
  clusterExport(cl = obj_clusters, varlist = c("df_table1","int_nFeatures1","df_one2all","var_pcorZ","bol_pcor","pcor.test","pcor.mat"), envir = environment())
  
  #info to user
  if(bol_pcor)
    print( paste0( "Starting partial correlation of ", nrow(df_table2), " ", str_labels[2], "s to ",
                   int_nFeatures1, " ", str_labels[1], "s within ", int_nSamples, " samples! ", date() ))
  else
    print( paste0( "Starting correlation of ", nrow(df_table2), " ", str_labels[2], "s to ",
                   int_nFeatures1, " ", str_labels[1], "s within ", int_nSamples, " samples! ", date() ))
  
  # Run the actual correlations
  lst_df_all2all <- t( parApply(cl = obj_clusters, X = df_table2, MARGIN = 1, FUN = df_one2all) )
  # lst_df_genes2OTUs <- t( apply(X = df_table2, MARGIN = 1, FUN = df_oneGene2allOTUs) )
  # clean up
  stopCluster(obj_clusters)
  
  #info to user
  print( paste("Finished all correlations! Preparing results.", date()) )
  
  # melt list of data-frames into one big data frame
  df_all2all <- do.call(rbind, lst_df_all2all)
  
  # convert feature ids into strings
  df_all2all[,str_labels[2]] <-  levels(fct_featureNames2)[ df_all2all[,str_labels[2]] ]
  df_all2all[,str_labels[1]] <-  levels(fct_featureNames1)[ df_all2all[,str_labels[1]] ]
  
  # correct p-values for multiple testing
  df_all2all$p.adj <- p.adjust(p = df_all2all$p.val,method = "BH")
  
  rownames(df_all2all) <- NULL
  
  return( df_all2all )
}
## Partial correlation function - from package ppcor v1.0 - added an NA remove option
pcor.test <- function(x,y,z,use="mat",method="p",na.rm=T){
  # The partial correlation coefficient between x and y given z
  #
  # pcor.test is free and comes with ABSOLUTELY NO WARRANTY.
  #
  # x and y should be vectors
  #
  # z can be either a vector or a matrix
  #
  # use: There are two methods to calculate the partial correlation coefficient.
  #	 One is by using variance-covariance matrix ("mat") and the other is by using recursive formula ("rec").
  #	 Default is "mat".
  #
  # method: There are three ways to calculate the correlation coefficient, 
  #	    which are Pearson's ("p"), Spearman's ("s"), and Kendall's ("k") methods.
  # 	    The last two methods which are Spearman's and Kendall's coefficient are based on the non-parametric analysis.
  #	    Default is "p".
  #
  # na.rm: If na.rm is T, then all the missing samples are deleted from the whole dataset, which is (x,y,z).
  #        If not, the missing samples will be removed just when the correlation coefficient is calculated.
  #	   However, the number of samples for the p-value is the number of samples after removing 
  #	   all the missing samples from the whole dataset.
  #	   Default is "T".
  
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  
  if(use == "mat"){
    p.use <- "Var-Cov matrix"
    pcor = pcor.mat(x,y,z,method=method,na.rm=na.rm)
  }else if(use == "rec"){
    p.use <- "Recursive formula"
    pcor = pcor.rec(x,y,z,method=method,na.rm=na.rm)
  }else{
    stop("\'use\' should be either \"rec\" or \"mat\"!\n")
  }
  
  # print the method
  if(gregexpr("p",method)[[1]][1] == 1){
    p.method <- "Pearson"
  }else if(gregexpr("s",method)[[1]][1] == 1){
    p.method <- "Spearman"
  }else if(gregexpr("k",method)[[1]][1] == 1){
    p.method <- "Kendall"
  }else{
    stop("\'method\' should be \"pearson\" or \"spearman\" or \"kendall\"!\n")
  }
  
  # sample number
  n <- dim(na.omit(data.frame(x,y,z)))[1]
  
  # given variables' number
  gn <- dim(z)[2]
  
  # p-value
  if(p.method == "Kendall"){
    statistic <- pcor/sqrt(2*(2*(n-gn)+5)/(9*(n-gn)*(n-1-gn)))
    p.value <- 2*pnorm(-abs(statistic))
    
  }else{
    statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
    p.value <- 2*pnorm(-abs(statistic))
    # p.value <- 2 * pt(-abs(statistic), (n - 2 - gn))
  }
  
  data.frame(estimate=pcor,p.value=p.value,statistic=statistic,n=n,gn=gn,Method=p.method,Use=p.use)
}			
## By using var-cov matrix
pcor.mat <- function(x,y,z,method="p",na.rm=T){
  
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  
  if(dim(z)[2] == 0){
    stop("There should be given data\n")
  }
  
  data <- data.frame(x,y,z)
  
  ## Conversion of factors to integers is required for next step to work
  for(int_i in 1:ncol(data)) {
    if(is.factor(data[,int_i]))
      data[,int_i] <- as.integer(droplevels(data[,int_i]))-1
  }
  
  
  if(na.rm == T){
    data = na.omit(data)
  }
  
  xdata <- na.omit(data.frame(data[,c(1,2)]))
  Sxx <- cov(xdata,xdata,m=method)
  xzdata <- na.omit(data)
  xdata <- data.frame(xzdata[,c(1,2)])
  zdata <- data.frame(xzdata[,-c(1,2)])
  Sxz <- cov(xdata,zdata,m=method)
  zdata <- na.omit(data.frame(data[,-c(1,2)]))
  Szz <- cov(zdata,zdata,m=method)
  
  # is Szz positive definite?
  zz.ev <- eigen(Szz)$values
  if(min(zz.ev)[1]<0){
    stop("\'Szz\' is not positive definite!\n")
  }
  
  # partial correlation
  Sxx.z <- Sxx - Sxz %*% solve(Szz) %*% t(Sxz)
  
  rxx.z <- cov2cor(Sxx.z)[1,2]
  
  rxx.z
}


#####
## Distance between GO-terms based of shared genes/features
mtx_mydist <- function(df_GOresultTable, mtx_assocScores = NULL) {
  lst_str_geneIDs <- strsplit(df_GOresultTable$geneID,"/")
  
  mtx_pairwiseDist <- matrix(data = 0, nrow = length(lst_str_geneIDs), ncol = length(lst_str_geneIDs))
  colnames(mtx_pairwiseDist) <- df_GOresultTable$ID
  rownames(mtx_pairwiseDist) <- df_GOresultTable$ID
  
  for(int_i in 1:(length(lst_str_geneIDs)-1)) {
    for(int_j in (int_i+1):length(lst_str_geneIDs)) {
      mtx_pairwiseDist[int_j,int_i] <- 1 - ( sum( lst_str_geneIDs[[int_i]] %in% lst_str_geneIDs[[int_j]] ) / 
                                               min(c(length(lst_str_geneIDs[[int_i]]), length(lst_str_geneIDs[[int_j]]))) )
      mtx_pairwiseDist[int_i,int_j] <- mtx_pairwiseDist[int_j,int_i]
    }
  }
  
  # add also the euclidian distance, if values are given in mtx_assocScores
  if(!is.null(mtx_assocScores)) {
    # check if names are matching
    if( sum(!colnames(mtx_assocScores) %in% colnames(mtx_pairwiseDist)) == 0 ) {
      # if it's columns that are matching we need to get orientation right  
      mtx_assocScores <- t(mtx_assocScores)
    } else if(sum(!rownames(mtx_assocScores) %in% colnames(mtx_pairwiseDist)) != 0) {
      warning("The provided names of mtx_assocScores and df_GOresultTable are not matching!\nIgnoring euclidian distances.")
      # break here
      return( mtx_pairwiseDist )
    }
    # else
      # print(all.equal(rownames(mtx_assocScores),rownames(mtx_pairwiseDist)))
    # euclidian dist only in tenths in order to get more finegrained clustering, but still have main clustering due to related GO terms
    mtx_pairwiseDist <- scales::rescale(mtx_pairwiseDist + as.matrix(dist(mtx_assocScores))[rownames(mtx_pairwiseDist),colnames(mtx_pairwiseDist)]/10,
                                 to = c(0,1))
  }
  
  return( mtx_pairwiseDist )
}

#####
## Association scores between two tables of GO-terms, which are based on the results of partial-correlation
## Annotated features from each GO-term are checked for their correlation-results to each feature within the second tables GO-terms
mtx_getAssocScores <- function(df_GOresultTable, df_SubsysResultTable, df_correlationResultTable) {
  mtx_assocScores <- matrix(data = NA, 
                            nrow = nrow(df_SubsysResultTable), 
                            ncol = nrow(df_GOresultTable))
  rownames(mtx_assocScores) <- df_SubsysResultTable$ID
  colnames(mtx_assocScores) <- df_GOresultTable$ID
  ## Now count the features of each GO-entry which have sign. correlations with any Subsy-entry
  df_tmp <- as.data.frame(
    lapply(X = df_SubsysResultTable[rownames(mtx_assocScores),"geneID"], FUN = function(str_subsys) {
      str_featureIDs <- unlist( strsplit(str_subsys, "/") )
      
      # for each feature-list, check count of sign. correlations with host-genes
      lst_assocScores <- unlist(
        lapply(X = df_GOresultTable[colnames(mtx_assocScores),"geneID"], FUN = function(str_go) {
          # get genes from GO-term
          str_geneIDs <- unique( df_mm10Ensembl2Genesymbol[df_mm10Ensembl2Genesymbol$GeneSymbol %in% 
                                                             unlist( strsplit(str_go, "/") ),
                                                           "EnsemblID"] )
          
          # Negative Spearman Rho AND feature-IDs of this subsystem AND feature-IDs of that GO-Term
          mtx_tmp <- df_correlationResultTable[df_correlationResultTable$Spearman_rho < 0 &
                                                 df_correlationResultTable$Rxn %in% str_featureIDs &
                                                 df_correlationResultTable$Gene %in% str_geneIDs,]
          # if there are any entries left
          if(!is.null(mtx_tmp) && (dim(mtx_tmp)[1] * dim(mtx_tmp)[2]) != 0) {
            # actual sign. correlations by total number of possible sign. correlations = 
            #  percentage of sign. cor. features between this pair of GOs, value between 0 to 1
            dbl_mean <- nrow(mtx_tmp)/(length(str_featureIDs) * length(str_geneIDs))

            return(dbl_mean)
          }
          else return(0)
        })
      )
      return(lst_assocScores)
    })
  ); mtx_assocScores <- t(df_tmp); rm(df_tmp)
  rownames(mtx_assocScores) <- df_SubsysResultTable$ID
  colnames(mtx_assocScores) <- df_GOresultTable$ID
  
  # discard actually uncorrelated features
  mtx_assocScores <- mtx_assocScores[rowSums(mtx_assocScores) != 0,
                                     colSums(mtx_assocScores) != 0, drop = F]
}

#####
## Algorithm recreated according to Figure 1 in publication:
## https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-4-267
## "GO Trimming: Systematically reducing redundancy in large Gene Ontology datasets" - 2011
df_trimGO <- function(df_GOresultTable, dbl_sofTrimPercent = 1/3) {
  require("GO.db")
  
  # sort by total genes in result query - fewest to most
  df_GOresultTable <- as.data.frame( df_GOresultTable[order(df_GOresultTable$Count, decreasing = F),] )
  # make sure rownames are GOIDs:
  rownames(df_GOresultTable) <- df_GOresultTable$ID
  
  # check for unknown GO-IDs
  df_droppedGOTerms <- NULL
  if(sum(!df_GOresultTable$ID %in% keys(GOBPPARENTS)) > 0) {
    warning("Unknown GO-Term IDs detected!")
    # make a backup copy of those terms
    df_droppedGOTerms <- df_GOresultTable[!df_GOresultTable$ID %in% keys(GOBPPARENTS),,drop = F]
    df_droppedGOTerms$familyID <- NA
    df_droppedGOTerms$softTrim <- FALSE
    # and remove them from the analysis table
    df_GOresultTable <- df_GOresultTable[df_GOresultTable$ID %in% keys(GOBPPARENTS),]
  }
  
  int_queryCount <- nrow(df_GOresultTable)
  # helper variable which holds the familyIDs of GO-relationships
  lst_familyIDs <- rep(NA,int_queryCount)
  names(lst_familyIDs) <- df_GOresultTable$ID
  lst_familyIDs <- as.list( lst_familyIDs )
  # print(lst_familyIDs)
  
  if(!is.null(df_GOresultTable$familyID)) {
    # TODO load that info into a variable actually :D
    lst_familyIDs <- df_GOresultTable$familyID
    print(lst_familyIDs)
    # Step 1 was already completed on input table skipping to Step 2  
    print("Step 1 was already completed on input table skipping to Step 2")
  }
  else {
    # Phase 1: Prepare data structure
    print("Step 1: Establishing GO relationships")
    int_idCounter <- 0
    # Get through query list one item after the other
    for(int_queryItem in 1:int_queryCount) {
      # define child node
      str_childID <- df_GOresultTable[int_queryItem,"ID"]
      # obtain direcet ancestor(s
        # lst_parentIDs <- try(unlist( as.list(GOBPPARENTS[str_childID]) ), silent = F)
        # # if no parents found, skip this one
        # if(class(lst_parentIDs) == "try-error") next
          # if(str_childID == "GO:0140972") lst_parentIDs <- c("GO:0031333") #c("GO:0140970","GO:0031333","GO:0140971","GO:0141086")
          # else 
      lst_parentIDs <- unlist( as.list(GOBPPARENTS[str_childID]) )
      lst_parentIDs <- lst_parentIDs[lst_parentIDs != "all"]
      # Now loop through ancestors
      # print(str_childID)
      cat(round(int_queryItem/int_queryCount*100,0), "%  ")
      repeat {
        # use first as parent
        str_parentID <- lst_parentIDs[1]
        # is parent in query list?
        if( str_parentID %in% df_GOresultTable$ID ) {
          # child has no ID yet
          if(is.na( lst_familyIDs[str_childID] )) {
            lst_familyIDs[[str_childID]] <- paste0("ID", int_idCounter)
            int_idCounter <- int_idCounter + 1
          }
          # parent has a different ID then child?
          if( is.na(lst_familyIDs[str_parentID]) |
              0 == sum( lst_familyIDs[[str_childID]] %in% lst_familyIDs[[str_parentID]] )
          ) {
            # then add childID to parents
            lst_familyIDs[[str_parentID]] <- c(na.omit(lst_familyIDs[[str_parentID]]), lst_familyIDs[[str_childID]])
          }
          # parent has a common ID with child
          else {
            # Get next ancestor(s)
              # lst_tmpParentIDs <- try(unlist( as.list(GOBPPARENTS[str_childID]) ), silent = F)
              # # if no parents found, skip this one
              # if(class(lst_tmpParentIDs) != "try-error") lst_parentIDs <- c(lst_parentIDs, lst_tmpParentIDs)
            lst_parentIDs <- c(lst_parentIDs, unlist( as.list(GOBPPARENTS[str_parentID]) ))
            # remove already checked element
            lst_parentIDs <- lst_parentIDs[lst_parentIDs != str_parentID]
            lst_parentIDs <- lst_parentIDs[lst_parentIDs != "all"]
          }
        }
        #
        else {
          # # Get next ancestor(s)
            # lst_tmpParentIDs <- try(unlist( as.list(GOBPPARENTS[str_childID]) ), silent = F)
            # # if no parents found, skip this one
            # if(class(lst_tmpParentIDs) != "try-error") lst_parentIDs <- c(lst_parentIDs, lst_tmpParentIDs)
          lst_parentIDs <- c(lst_parentIDs, unlist( as.list(GOBPPARENTS[str_parentID]) ))
          # remove already checked element
          lst_parentIDs <- lst_parentIDs[lst_parentIDs != str_parentID]
          lst_parentIDs <- lst_parentIDs[lst_parentIDs != "all"]
        }
        if(length(lst_parentIDs) == 0){
          break
        }
      }
    }
    # clean up before reusing variables
    rm(str_childID, str_parentID, lst_parentIDs, int_idCounter, int_queryItem)
  }
  
  # Phase 2: Analyse termID-structure
  cat("\nStep 2: Selecting most relevant terms\n")  
  # helper variable to denote which term to trim
  lst_softTrimTerm <- rep(F, int_queryCount)
  names(lst_softTrimTerm) <- df_GOresultTable$ID
  # variable of discarded terms
  lst_discardedIDs <- NULL
  
  # for each item in query list but last:
  for(int_queryItem in 1:(int_queryCount-1)) {
    cat(round(int_queryItem/int_queryCount*100,0), "%  ")
    # Are there some familyIDs for this query Item?
    if( !is.na(lst_familyIDs[int_queryItem]) ) {
      # define child node
      str_childID <- names(lst_familyIDs)[int_queryItem]
      # print(str_childID)
      
      # are there items with a similar ID after current query-position
      lst_splitTrimIds <- lst_familyIDs[(int_queryItem+1):int_queryCount]
      for(int_i in 1:length(lst_splitTrimIds)) {
        if(!is.na( lst_splitTrimIds[int_i] ) & 
           1 >= sum( lst_familyIDs[[str_childID]] %in% lst_splitTrimIds[[int_i]] ) 
        ) {
          # is this term an ancestor of child Term
          lst_parentIDs <- unlist( as.list(GOBPANCESTOR[str_childID]) )
          if( !is.null(lst_parentIDs) & length(lst_parentIDs) >= 1 & 
              names(lst_familyIDs)[int_queryItem + int_i] %in% lst_parentIDs ) {
            # if ancestor use as parent
            str_parentID <- names(lst_familyIDs)[int_queryItem + int_i]
            
            # are genes of parent term exactly the same (aka included completely) in child term?
            int_sharedGenes <- sum( strsplit(df_GOresultTable[str_parentID,"geneID"],"/")[[1]] %in% 
                                      strsplit(df_GOresultTable[str_childID,"geneID"],"/")[[1]] )
            # print(int_sharedGenes)
            # print(paste(str_parentID, str_childID))
            if( 1 == (int_sharedGenes / df_GOresultTable[str_parentID,"Count"]) ) {
              # remove parent immediately
              df_GOresultTable <- df_GOresultTable[!rownames(df_GOresultTable) %in% str_parentID,]
              lst_discardedIDs <- c(lst_discardedIDs, str_parentID)
            }
            else {
              # print(paste0("parent:",df_GOresultTable[str_parentID,"Count"]-int_sharedGenes))
              # print(paste0("child:",df_GOresultTable[str_childID,"Count"]-int_sharedGenes))
              # print(paste0("shared:",int_sharedGenes))
              
              # is child NOT marked for optional trimming?
              # AND does parent not have at least 40% more genes than child?
              # AND parent unique genes to shared genes is less than 1/3
              if(!lst_softTrimTerm[str_childID] &
                 # dbl_sofTrimPercent > ( (df_GOresultTable[str_parentID,"Count"] - df_GOresultTable[str_childID,"Count"]) / df_GOresultTable[str_childID,"Count"] )
                 dbl_sofTrimPercent < ( (df_GOresultTable[str_parentID,"Count"] - int_sharedGenes - 1) / (int_sharedGenes + 1) )
              ) {
                # then mark parent for optional trimming
                lst_softTrimTerm[str_parentID] <- T
              }
            }
          }
        }
      }
    }
  }
  if(!is.null(lst_discardedIDs) & length(lst_discardedIDs) > 0) {
    cat("\nFollowing terms had to be removed: ", lst_discardedIDs, "\n")
  }
  # Add softTrim info and familyIDs as columns to result table output
  df_GOresultTable$familyID <- lst_familyIDs[df_GOresultTable$ID]
  df_GOresultTable$softTrim <- lst_softTrimTerm[df_GOresultTable$ID]
  
  # in case we've dropped unknown terms from the analysis before, now add them back
  if(!is.null(df_droppedGOTerms)) {
    df_GOresultTable <- rbind(df_GOresultTable,
                              df_droppedGOTerms)
  }
  
  cat("\nFinished with GO trimming!\n")
  return(df_GOresultTable)
}

#####
## Manually evaluate GO-terms or subsystems which are redundant
df_checkGO <- function(lst_checkGOIDs, df_GOresultsReference, bol_checkAllTerms = F) {
  # What is this function doing and how to use
  cat("Interactively currate a gene set enrichment results table:\n", 
      "(Press [Ctrl]+[D] any time to abort)\n\n")
  
  # build internal reference table
  df_GOresultsReference <- df_GOresultsReference[!duplicated(df_GOresultsReference[,c(1,2,4,8)]),]
  df_GOresultsReference <- df_GOresultsReference[df_GOresultsReference$ID %in% lst_checkGOIDs,]
  # initiate response list
  lst_dropIDs <- NULL
  
  ###
  # function that goes through duplicated GO entries
  internalCheckDuplicates <- function(lst_rowIndex) {
    if(length(lst_rowIndex) > 1) {
      bol_invalid = T
      while(bol_invalid) {
        cat("Those terms are duplicates. Select one?\n")
        # print terms with number indices
        cat( paste0(paste0("[", 2:(length(lst_rowIndex)+1), "]\t"), 
                    df_GOresultsReference[lst_rowIndex,"Description"]), 
             sep = "\n")
        # keep all option and new label option
        cat( "\n[0]\tKeep all terms unchanged\n[1]\tSelect one new label for all" )
        
        # get user selection
        chr_key <- as.integer( scan("", what = "integer", nmax = 1, nlines = 1, flush = T, quiet = T) )
        
        ## Parse user key input
        # keep all unchanged
        if(chr_key == 0) {
          cat("Unchanged\n\n")
          bol_invalid = F
          return( NULL )
        }
        # new custom label representing all terms
        else if(chr_key == 1) { 
          str_termDescription <- readline(prompt="Enter custom feature name: ")
          bol_discard <- rep(T,length(lst_rowIndex))
          bol_discard[1] <- F
          df_GOresultsReference[lst_rowIndex[1],"Description"] <- paste0(str_termDescription, " *")
          assign(x = "df_GOresultsReference", value = df_GOresultsReference, pos = parent.env(env = environment()))
          cat("Keep:\t", df_GOresultsReference[lst_rowIndex[1],"Description"],"\n\n")
          bol_invalid = F
          return( df_GOresultsReference[lst_rowIndex[ bol_discard ],"ID"] )
        }
        # select one specific term
        else if(chr_key > 1 & chr_key <= length(lst_rowIndex)+1) {
          chr_key <- chr_key-1
          cat("Keep:\t", df_GOresultsReference[lst_rowIndex[chr_key],"Description"],"\n\n")
          bol_discard <- rep(T,length(lst_rowIndex))
          bol_discard[chr_key] <- F
          bol_invalid = F
          return( df_GOresultsReference[lst_rowIndex[ bol_discard ],"ID"] )
        }
        else {
          cat("Invalid\n")
          bol_invalid = T
        }
      }
    }
    return( NULL )
  }
  # function that prompts for each term to keep or drop it
  internalPromptTerm <- function(str_intTermID) {
    bol_invalid = T
    while(bol_invalid) {
      cat("Keep this term?\n")
      # print terms with number indices
      cat( unique( df_GOresultsReference[df_GOresultsReference$ID == str_intTermID,"Description"] ), sep = "\n")
      # Get user response
      cat("[y] keep\t[x] discard\t[Ctrl] + [D]\tabort\n")
      chr_key <- stringr::str_to_lower( scan("", what = "character", nmax = 1, nlines = 1, flush = T, quiet = T) )
      
      # Parse user response
      if(chr_key == "y") {
        cat("Keep\n\n")
        return( NULL )
      }
      else if(chr_key == "x") {
        cat("Drop\n\n")
        return(str_intTermID)
      }
      else {
        cat("Invalid\n")
        bol_invalid = T
      }
    }
  }
  ###
  
  # check for duplicates first
  lst_dropIDs  <- unlist( tapply(X = 1:nrow(df_GOresultsReference),
                                 INDEX = as.factor( paste0(df_GOresultsReference$BgRatio, "_", df_GOresultsReference$geneID) ),
                                 FUN = internalCheckDuplicates) )
  names(lst_dropIDs) <- NULL
  # print(lst_dropIDs)
  
  # Drop those terms from original input list
  lst_checkGOIDs <- lst_checkGOIDs[!lst_checkGOIDs %in% lst_dropIDs]
  
  # update internal reference table
  df_GOresultsReference <- df_GOresultsReference[df_GOresultsReference$ID %in% lst_checkGOIDs,]
  df_GOresultsReference <- df_GOresultsReference[!duplicated(df_GOresultsReference[,c(1,2,4,8)]),]

  # Now go through and ask if we want to keep that term due to context
  if(bol_checkAllTerms) {
    for(str_termID in lst_checkGOIDs) {
      lst_dropIDs <- c(lst_dropIDs, internalPromptTerm(str_termID))
    }
  }
  
  # Write ref-table output to global environment
  cat("Reference table was written to \"df_GOresultsReference\"\n\n")
  df_GOresultsReference <<- df_GOresultsReference
  View(df_GOresultsReference)
  return( lst_dropIDs )
}

#####
## plotting function - dot-size represents numbers in a matrix
myDotPlot <- function(mtx_input, scaleFrom = NULL, scaleTo = c(0.7, 2.8), labelCex = .75, color = "#b0c4de77") {
  # Size of points, indicate strenght of host to microbiome associations
  if(is.null(scaleFrom)) scaleFrom <- c(min(abs(as.numeric(mtx_input))), max(abs(as.numeric(mtx_input))))
  dbl_size <- abs(as.numeric(mtx_input))
  dbl_size[dbl_size>0] <- scales::rescale(dbl_size[dbl_size>0], 
                                          to = scaleTo,
                                          from = scaleFrom)
  #
  plot(y = NA,
       x = NA,
       ylim = c(0,ncol(mtx_input)),
       xlim = c(0,nrow(mtx_input)),
       cex = NA,
       pch = 21,
       # bg = c("#D6604D","#4393C3")[(as.numeric(mtx_input) > 0)+1],
       bg = color,
       axes = F, frame.plot = F,
       xlab = "", ylab = "")
  # grid
  abline(h = 1:ncol(mtx_input)-.5,
         v = 1:nrow(mtx_input)-.5, 
         col = rgb(.5,.5,.5,.4), lwd = .7)
  # points
  points(y = rep(1:ncol(mtx_input), each = nrow(mtx_input))-.5,
         x = rep(1:nrow(mtx_input), times = ncol(mtx_input))-.5,
         cex = dbl_size,
         pch = 21,
         bg = color )
  # X-Axis
  axis(side = 1, at = 1:nrow(mtx_input)-.5, line = -.5,
       labels = rep("", nrow(mtx_input)),
       las = 2, cex.axis = labelCex)
  text(x = 1:nrow(mtx_input)-.5, y = -1, xpd = T,
       labels = rownames(mtx_input),
       srt = 30, adj = c(1,0.5), cex = labelCex)
  # Y-Axis
  axis(side = 2, at = 1:ncol(mtx_input)-.5, line = 0, 
       labels = colnames(mtx_input), 
       las = 1, cex.axis = labelCex)
}

#####
## Create MDS-Plot
## df_metaData = dataFrame with meta information, rownames reflect column names of df_sampleMatrix
## var_colorBy = integer or character to give metaData column for coloring information
## str_mainTitle = Title on the plot
## bol_labels = boolean indicating wether to print sample labels or not (default: False)
## ... = any graphics parameter
void_MDSplot <- function(df_sampleMatrix, df_metaData = NULL, bol_labels = FALSE, var_colorBy = 1, var_shapeBy = 1, bol_continousColor = T, str_mainTitle = "Multi Dimensional Scaling", str_pdfFileName = NULL, ...) {
  require("MASS")
  require("vegan")
  require("RColorBrewer")
  
  #mds uses rows as samples therefore we need to transform our data and switch rows with cols
  obj_abundanceDistances <- dist(t( df_sampleMatrix ))
  # print(obj_abundanceDistances)
  #calculate an mds plot (like pca but including 0 values) multidimensional scaling
  obj_abundanceMDS <- monoMDS(obj_abundanceDistances, k = 2)
  #replaced isoMDS with monoMDS which is a more default method
  # isoMDS(obj_abundanceDistances, k = 2)
  
  dbl_stress <- round( c(monoMDS(obj_abundanceDistances, k = 1)$stress, 
                         obj_abundanceMDS$stress,
                         monoMDS(obj_abundanceDistances, k = 3)$stress)*100, 1)
  
  # there's no such thing as explained variance in non-metric MDS plots
  # however the stress variable goes a bit in this direction
  # -> see https://stackoverflow.com/questions/49223740/cumulative-variance-explained-for-nmds-in-r
  cat( paste0( "Stress for 1-dimensional nMDS: ", dbl_stress[1], 
               "\nStress for 2 dimensions: ", dbl_stress[2],
               "\n3 dimensions: ", dbl_stress[3], "\n..." ))
  
  if(is.null(df_metaData)) {
    obj_colors <- as.factor(rep(0,nrow(obj_abundanceMDS$points)))
    obj_shapes <- factor(rep(16,nrow(obj_abundanceMDS$points)), levels = 1:16)
  }
  else {
    #select rows from metaData via columnnames of original input
    if(is.null(var_shapeBy) | isFALSE(var_shapeBy)) {
      obj_shapes <- factor(rep(16,nrow(obj_abundanceMDS$points)), levels = 1:16)
    }
    else {
      obj_shapes <- factor(df_metaData[ rownames(obj_abundanceMDS$points) , var_shapeBy], levels = 1:16)
    }
    # 
    # print(df_metaData[ rownames(obj_abundanceMDS$points) , "age"])
    # print(obj_shapes)
    # print(var_shapeBy)
    # 
    #select rows from metaData via columnnames of original input
    obj_colors <- as.factor(df_metaData[ rownames(obj_abundanceMDS$points) , var_colorBy])
    
    #in case of a continous numeric variable with more than 8 different options,
    # cast to sub-groups and stratify colors by those sub-groups
    if(length(levels(obj_colors)) > 8) {
      print("Stratifying variable has more than 8 levels and will thus be sub-grouped.")
      
      # get max limits for each sub-set
      dbl_upperLimits <- quantile(as.numeric( df_metaData[ , var_colorBy] ),
                                  probs = c(0.125,0.25,0.375,0.5,0.625,0.75,0.875,1), na.rm = T)
      
      # introduce a new metaData variable
      df_metaData[ , paste0("Subgrouped_",var_colorBy)] <- NA
      
      # in case of NA (missing values) in original metadata:
      df_metaData[is.na(df_metaData[,var_colorBy]), paste0("Subgrouped_",var_colorBy)] <- "missing"
      df_metaData[is.na(df_metaData[,var_colorBy]), var_colorBy] <- min(as.numeric( df_metaData[, var_colorBy] ), na.rm = T) - 1000
      
      # add a label for each sample based on its value
      lst_levelOrder <- c("missing")
      for(int_i in 1:8) {
        df_metaData[ is.na(df_metaData[,paste0("Subgrouped_",var_colorBy)]) &
                       df_metaData[,var_colorBy] <= dbl_upperLimits[ int_i ],
                     paste0("Subgrouped_",var_colorBy)] <- paste0( "< ", round(x = dbl_upperLimits[ int_i ], digits = 2) )
        lst_levelOrder <- c(lst_levelOrder, paste0( "< ", round(x = dbl_upperLimits[ int_i ], digits = 2) ))
      }
      
      # MetaData variable is now changed
      var_colorBy = paste0("Subgrouped_",var_colorBy)
      
      # create a new colouring variable based of the stratified sub-groups
      obj_colors <- factor(x = df_metaData[ rownames(obj_abundanceMDS$points) , var_colorBy], levels = lst_levelOrder, ordered = T)
      # discard any unused levels
      obj_colors <- droplevels(obj_colors)
    }
  }
  
  # coloring in distinct colors or one continous saturation?
  if(bol_continousColor) col_brewer <- brewer.pal(n = length(levels(obj_colors)), name = "GnBu")
  else { 
    if(length(levels(obj_colors)) > 9)
      col_brewer <- colorRampPalette(brewer.pal(n = 9, name = "Set1"))(length(levels(obj_colors)))
    else
      col_brewer <- brewer.pal(n = length(levels(obj_colors)), name = "Set1")
  }
  
  #In case a PDF-Filename is given plot to file
  if(!is.null(str_pdfFileName)) pdf(file = str_pdfFileName, width = 12, height = 9, paper='special')
  
  # print(as.numeric(obj_shapes))
  
  plot(obj_abundanceMDS$points[,c(1,2)],
       col = col_brewer[ obj_colors ],
       # pch= as.numeric( obj_shapes) +18, 
       # pch = ((as.numeric(obj_shapes)-1)*16)+1,
       pch = sort( as.numeric( unique( obj_shapes))-1 ),
       cex=.9,
       xlab = paste0("Dimension 1 (", dbl_stress[1], " %)"), 
       ylab = paste0("Dimension 2 (", dbl_stress[2], " %)"),
       main = str_mainTitle,
       ...)
  
  #print sample labels if asked for
  if(bol_labels) {
    # .. and plot it a second time with legend
    plot(obj_abundanceMDS$points[,c(1,2)],
         col = col_brewer[ obj_colors ],
         # pch= as.numeric( obj_shapes) +18, 
         # pch = ((as.numeric(obj_shapes)-1)*16)+1,
         pch = sort( as.numeric( unique( obj_shapes)) ),
         cex=.9,
         xlab = paste0("Dimension 1 (", dbl_stress[1], " %)"), 
         ylab = paste0("Dimension 2 (", dbl_stress[2], " %)"),
         main = str_mainTitle,
         ...)
    text(x = obj_abundanceMDS$points[,1], 
         y = obj_abundanceMDS$points[,2]-(max(obj_abundanceMDS$points[,2])/10), 
         labels = colnames(df_sampleMatrix),
         xpd = NA)
  }
  
  # .. and plot it a second time with legend
  plot(obj_abundanceMDS$points[,c(1,2)],
       col = col_brewer[ obj_colors ],
       # pch= as.numeric( obj_shapes) +18, 
       # pch = ((as.numeric(obj_shapes)-1)*16)+1,
       pch = sort( as.numeric( unique( obj_shapes)) ),
       cex=.9,
       xlab = paste0("Dimension 1 (", dbl_stress[1], " %)"), 
       ylab = paste0("Dimension 2 (", dbl_stress[2], " %)"),
       main = str_mainTitle,
       ...)
  #print sample labels if asked for
  if(bol_labels) text(x = obj_abundanceMDS$points[,1], 
                      y = obj_abundanceMDS$points[,2]-(max(obj_abundanceMDS$points[,2])/10), 
                      labels = colnames(df_sampleMatrix), xpd = NA)
  # text(x = obj_abundanceMDS$points[,1], 
  #      y = obj_abundanceMDS$points[,2]-(max(obj_abundanceMDS$points[,2])/10), 
  #      labels = obj_colors)
  legend("bottomright",
         legend = levels( droplevels(obj_colors) ),
         col = col_brewer[ sort(as.numeric(unique(obj_colors))) ],
         bg = "white",
         pch=19,
         title = colnames(df_metaData[var_colorBy]),
         cex = 0.8)
  
  if(!is.null(var_shapeBy) & !isFALSE(var_shapeBy)) {
    legend("topright",
           legend = levels( droplevels( obj_shapes )),
           col = 1,
           bg = "white",
           # pch = sort( unique( as.numeric( obj_shapes))),
           pch = sort( as.numeric( unique( obj_shapes)) ),
           title = colnames(df_metaData[var_shapeBy]),
           cex = 0.8)
  }
  # rm(obj_colors,obj_abundanceDistances,obj_abundanceMDS)
  if(!is.null(str_pdfFileName)) dev.off()
}

#####
## run several linear models for each row in dependentVariable dataframe and correct for multiple testing
df_linearModel <- function(df_dependantVariable, lst_var_metaData, int_cores = 7) {
  require(parallel)
  
  ## make sure dependant data is numeric
  mtx_dependantVariable <- data.matrix( frame = df_dependantVariable )
  
  ## retain rownames
  # add a column for the rownames (can later be used to get them back)
  # this column needs to be numeric as well in order to be succesfully passed to a function!
  # thus we employ factors
  int_colNum <- ncol(mtx_dependantVariable)
  str_rowNames <- rownames(df_dependantVariable)
  mtx_dependantVariable <- cbind(mtx_dependantVariable, c(1:nrow(mtx_dependantVariable)))
  
  #Internal function for calling multiple correlations in parrallel mode
  df_internalLinearModel <- function(row_dependantVariable) {
    int_rowName <- row_dependantVariable[int_colNum + 1]
    row_dependantVariable <- as.numeric( row_dependantVariable[1 : int_colNum] )
    
    # run the wilcox test
    #H0: row_dependantVariable (e.g. Expression of Gene Xyz) is independent of lst_var_metaData (e.g. Age)
    #H1: row_dependantVariable (e.g. Expression of Gene Xyz) is altered by the samples lst_var_metaData (e.g. Age)
    try( 
      obj_linearModel <- summary( lm( row_dependantVariable ~ lst_var_metaData, na.action = "na.omit") )
    )
    if( exists("obj_linearModel") & dim(obj_linearModel$coefficients)[1] >= 2 ) {
      # print(obj_linearModel)
      return( c(int_rowName, obj_linearModel$coefficients[2,]) )
    }
    else print(paste("Unable to test dependant variable row",str_rowNames[int_rowName],"!"))
    
    return( c(int_rowName,0,0,0,1) )
  }
  
  # start working in parrallel
  obj_clusters <- makeCluster(int_cores, type="SOCK")
  clusterExport(cl = obj_clusters, varlist = c("mtx_dependantVariable","int_colNum","lst_var_metaData","str_rowNames","df_internalLinearModel"), envir = environment())
  
  #info to user
  print( paste( "Starting rank sum tests of",nrow(df_dependantVariable),"observations within",int_colNum,"samples!", date() ))
  
  # run rank sum tests in parrallel
  df_resultsLinearModel <- as.data.frame( t( parApply(cl = obj_clusters, X = mtx_dependantVariable, MARGIN = 1, FUN = df_internalLinearModel) ))
  # df_resultsLinearModel <- as.data.frame( t( apply(X = mtx_dependantVariable, MARGIN = 1, FUN = df_internalLinearModel) ))
  
  # clean up
  stopCluster(obj_clusters)
  
  #info to user
  print( paste("Finished all tests! Preparing results.", date()) )
  
  colnames(df_resultsLinearModel) <- c("n","est","std.err","t.value","p.value")
  
  # correct p-values for multiple testing
  df_resultsLinearModel$p.adj <- p.adjust(p = df_resultsLinearModel$p.value,method = "BH")
  
  # get rownames
  rownames(df_resultsLinearModel) <- str_rowNames[df_resultsLinearModel$n]
  df_resultsLinearModel$n <- NULL
  
  return( df_resultsLinearModel )
}

#####
## Get level 2 Gene Ontology parent term for provided GO biological processes GO:ID
df_getLevel2GO <- function(lst_str_GOIDs) {
  # Bioconductor package GO.db has environment variables used here:
  require("GO.db")
  internalGetParentGO <- function(str_inputGO) {
    str_nextGO <- str_inputGO
    lst_GOtree <- list(L = str_nextGO)
    while(T) {
      
      lst_parentIDs <- tryCatch({
        GOBPPARENTS[str_nextGO]
      }, 
      error=function(cond){ 
        return(NA)
      })
      lst_parentIDs <- unique( unlist( as.list(lst_parentIDs) ))
      
      lst_GOtree <- c(lst_GOtree, 
                      L = list(lst_parentIDs))
      
      str_nextGO <- lst_parentIDs[lst_parentIDs != "all"]
      
      # Stop when we reach the top with "Biological Process"
      if(length(str_nextGO) == 1) {
        if(is.na(str_nextGO) | str_nextGO == "GO:0008150") {
          names(lst_GOtree) <- paste0("L", length(lst_GOtree):1)
          # # layer 2
          # print(Term( GOTERM[lst_GOtree$L2] ))
          # # layer 4
          # print(Term( GOTERM[lst_GOtree$L4] ))
          
          # remove NAs
          str_outputGO <- na.omit( lst_GOtree$L2[!lst_GOtree$L2 %in% c("all","GO:0008150")] )
          
          if(isEmpty(str_outputGO)) {
            return( c(str_inputGO,  # Input GO-ID
                      "",           # Output GO-IDs
                      "",           # Input GO-Term
                      ""            # Output GO-Terms
            )
            )
          }
          if(str_inputGO == "GO:0008150") {
            return( c(str_inputGO,  # Input GO-ID
                      "",           # Output GO-IDs
                      "",           # Input GO-Term
                      ""            # Output GO-Terms
            )
            )
          }
          else if(str_outputGO[1] == str_inputGO) {
            return( c(str_inputGO,  # Input GO-ID
                      "",           # Output GO-IDs
                      "",           # Input GO-Term
                      ""            # Output GO-Terms
            )
            )
          }
          else {
            # filter out "biological process" parent term in case it's also listed in Layer2
            return( c(str_inputGO,                                        # Input GO-ID
                      paste0(collapse = ",", str_outputGO),               # Output GO-IDs
                      Term( GOTERM[str_inputGO] ),                        # Input GO-Term
                      paste0(collapse = ",", Term( GOTERM[str_outputGO])) # Output GO-Terms
            )
            )
          }
          break();
        }
      }
    }
  }
  df_output <- t( as.data.frame( lapply(X = lst_str_GOIDs, FUN = internalGetParentGO) ))
  df_output <- as.data.frame(df_output)
  rownames(df_output) <- NULL
  colnames(df_output) <- c("QueryID", "L2ParentID", "QueryTerm", "L2ParentTerm")
  
  return( df_output )
}

#####
## Lookup function
## str_gene is a genesymbol from mouse genome
## str_organ not implemented yet
void_checkGene <- function(str_gene, str_organ = "all") {
  # sanitize gene input
  str_gene <- stringi::stri_trans_totitle( stringi::stri_trans_tolower(str_gene))
  
  str_geneID <- df_mm10Ensembl2Genesymbol[df_mm10Ensembl2Genesymbol$GeneSymbol == str_gene, "EnsemblID"]
  
  # colon
  df_tmp <- merge(x = df_pcorRxn2ColonRNATop[df_pcorRxn2ColonRNATop$Gene == str_geneID,],
                  y = df_rxnAnnotation[,c(1:2)], by.x = "Rxn", by.y = "RxnID", all.x = T, all.y = F)
  if(nrow(df_tmp) > 0) {
    df_tmp <- df_tmp[order(df_tmp$Spearman_rho, decreasing = T),]
    rownames(df_tmp) <- NULL
    print(paste0(str_gene, " in Colon:"))
    print(df_tmp[,c("RxnName","Spearman_rho","p.adj")])
    View(df_tmp[,c("RxnName","Spearman_rho","p.adj")], title = paste0("Colon_",str_gene) )
  }
  
  # liver
  df_tmp <- merge(x = df_pcorRxn2LiverRNATop[df_pcorRxn2LiverRNATop$Gene == str_geneID,],
                  y = df_rxnAnnotation[,c(1:2)], by.x = "Rxn", by.y = "RxnID", all.x = T, all.y = F)
  if(nrow(df_tmp) > 0) {
    df_tmp <- df_tmp[order(df_tmp$Spearman_rho, decreasing = T),]
    rownames(df_tmp) <- NULL
    print(paste0(str_gene, " in Liver:"))
    print(df_tmp[,c("RxnName","Spearman_rho","p.adj")])
    View(df_tmp[,c("RxnName","Spearman_rho","p.adj")], title = paste0("Liver_",str_gene) )
  }
  
  # brain
  df_tmp <- merge(x = df_pcorRxn2BrainRNATop[df_pcorRxn2BrainRNATop$Gene == str_geneID,],
                  y = df_rxnAnnotation[,c(1:2)], by.x = "Rxn", by.y = "RxnID", all.x = T, all.y = F)
  if(nrow(df_tmp) > 0) {
    df_tmp <- df_tmp[order(df_tmp$Spearman_rho, decreasing = T),]
    rownames(df_tmp) <- NULL
    print(paste0(str_gene, " in Brain:"))
    print(df_tmp[,c("RxnName","Spearman_rho","p.adj")])
    View(df_tmp[,c("RxnName","Spearman_rho","p.adj")], title = paste0("Brain_",str_gene) )
  }
}

#####
## Simplify GO-Terms or Subsys-IDs for plotting
## str_description is a list of GO-terms or subsystem-names
str_simplifyDescriptions <- function(str_description, bol_stripRomanNums = T, int_maxLength = 65, bol_cut = F) {
  str_descriptionShortened <- str_description
  
  # remove too much details in brackets from subsys names
  str_descriptionShortened<- gsub(pattern = r"{\s*\([^\)]+\)}", replacement = "", x = str_descriptionShortened, perl = T)
  str_descriptionShortened<- sub(pattern = "^[-]", replacement = "", x = str_descriptionShortened, perl = T)
  
  # will remove things like Type I, Type II ...
  if(bol_stripRomanNums) str_descriptionShortened<- gsub(pattern = "type", replacement = "", x = str_descriptionShortened, fixed = T)

  # abbreviate too long texts
  str_descriptionShortened<- gsub(pattern = "biosynthesis", replacement = "synth.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "biosynthetic", replacement = "synth.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "degradation", replacement = "degr.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "detoxification", replacement = "detox.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "metabolism", replacement = "metab.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "pathway", replacement = "pwy.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "deoxyribonucleosides", replacement = "dRNs", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "deoxyribonucleoside", replacement = "dRN", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "deoxyribonucleotides", replacement = "dNMPs", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "deoxyribonucleotide", replacement = "dNMP", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "nucleotides", replacement = "nucl.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "nucleotide", replacement = "nucl.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "metabolites", replacement = "metabs.", x = str_descriptionShortened, fixed = T)
  #
  str_descriptionShortened<- gsub(pattern = "^regulation of", replacement = "", x = str_descriptionShortened, perl = T)
  str_descriptionShortened<- gsub(pattern = "positive regulation of ", replacement = "", x = str_descriptionShortened, fixed = T)
  #
  str_descriptionShortened<- gsub(pattern = "phosphorylation", replacement = "phos.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "small subunit", replacement = "SSU", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "establishment", replacement = "estab.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "maintenance", replacement = "maint.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "localization", replacement = "loc.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "formation", replacement = "form.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "initiation", replacement = "init.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "electron", replacement = "e-", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "transfer", replacement = "transf.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "transport", replacement = "transp.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "division", replacement = "div.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "acetylation", replacement = "acet.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "termination", replacement = "term.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "transcription", replacement = "transcrpt.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "cascade", replacement = "casc.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "templated", replacement = "templ.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "dependent", replacement = "dep.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "regulation", replacement = "reg.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "negative", replacement = "neg.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "positive", replacement = "pos.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "response", replacement = "resp.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "production", replacement = "prod.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "replication", replacement = "repl.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "differentiation", replacement = "diff.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "developmental", replacement = "devel.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "development", replacement = "devel.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "synthesis", replacement = "synth.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "organization", replacement = "org.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "exchange", replacement = "ex.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "pathway", replacement = "pwy.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "system", replacement = "sys.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "processing", replacement = "proc.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "process", replacement = "proc.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "signaling", replacement = "sign.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "signal", replacement = "sign.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "toxicity", replacement = "tox.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "activity", replacement = "act.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "activation", replacement = "act.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "stimulation", replacement = "stim.", x = str_descriptionShortened, fixed = T)
  str_descriptionShortened<- gsub(pattern = "concentration", replacement = "conc.", x = str_descriptionShortened, fixed = T)
  
  # discard latin numbers
  if(bol_stripRomanNums) str_descriptionShortened<- gsub(pattern = "[IXV]", replacement = "", ignore.case = F, x = str_descriptionShortened, perl = T)
  
  # get rid of too many white space chars
  str_descriptionShortened <- trimws(str_descriptionShortened, which = "both")
  # drop double dashes or double spaces
  for(int_i in 1:10) {
    str_descriptionShortened<- gsub(pattern = "--", replacement = "-", x = str_descriptionShortened, fixed = T)
    str_descriptionShortened<- gsub(pattern = "  ", replacement = " ", x = str_descriptionShortened, fixed = T)
  }
  
  # cut lines after defined number of chars
  if(bol_cut) {
    for(int_i in 1:length(str_descriptionShortened)) {
      if(stringi::stri_length(str_descriptionShortened[int_i]) >= int_maxLength) {
        str_descriptionShortened[int_i] <- paste0( stringi::stri_sub(str = str_descriptionShortened[int_i], from = 1, to = int_maxLength), " ...")
      }
    }
  }
  # OR split too long strings over two lines
  else {
    for(int_i in 1:length(str_descriptionShortened)) {
      if(stringi::stri_length(str_descriptionShortened[int_i]) >= int_maxLength) {
        int_striMid <- floor(stringi::stri_length(str_descriptionShortened[int_i])/2)
        str_descriptionShortened[int_i] <- sub(pattern = stringi::stri_sub(str_descriptionShortened[int_i], 
                                                                           from = int_striMid-10, to = int_striMid+15), 
                                               replacement = sub(pattern = " ", replacement = "\n", 
                                                                 x = stringi::stri_sub(str_descriptionShortened[int_i],
                                                                                       from = int_striMid-10, to = int_striMid+15), fixed = T), 
                                               x = str_descriptionShortened[int_i], fixed = T)
      }
    }
  }
  
  return( str_descriptionShortened )
}


####################################
### Some annotation databases

########################
## Translation table for ENSEMBL IDs to genesymbols for mouse
df_mm10Ensembl2Genesymbol <- read.table(file = "databases/GRCm38.102.EnsemblToGeneSymbol.tsv.gz", stringsAsFactors = F)
colnames(df_mm10Ensembl2Genesymbol) <- c("EnsemblID", "GeneSymbol")
rownames(df_mm10Ensembl2Genesymbol) <- df_mm10Ensembl2Genesymbol$EnsemblID


########################
## for cluster profiler
df_rxn2subsys <- readRDS("databases/df_rxn2subsys20230510MetaMouse.rds")


#########################
# ## Mouse gene Ontology (newest version as of May 2023)
df_mouseGOannotaion <- readRDS("databases/df_mouseGOannotaion_2023-05-16.rds")
## Now let's just select GOs with meaningful interaction terms and only biological processes:
df_mouseGOannotaionBioProc <-  df_mouseGOannotaion[df_mouseGOannotaion$OntologyType == "biological process" &
                                                     rownames(df_mouseGOannotaion) %in% rownames(df_mouseGOannotaion)[grep(pattern = "negative", x = df_mouseGOannotaion$FunctionalRole, invert = T)] &
                                                     rownames(df_mouseGOannotaion) %in% rownames(df_mouseGOannotaion)[grep(pattern = "NOT", x = df_mouseGOannotaion$FunctionalRole, invert = T)],]
# manually need to add this term, as it pops up and had no description so far
df_mouseGOannotaionBioProc[df_mouseGOannotaionBioProc$GOid == "GO:0006342","Description"] <- "heterochromatin assembly"
rm(df_mouseGOannotaion)


#########################
# annotation helper table to convert GO-IDs into text
df_GOdescriptions <- na.omit(unique(df_mouseGOannotaionBioProc[,c(1,6)]))
rownames(df_GOdescriptions) <- df_GOdescriptions$GOid
####
# and for subsystems
df_SubSysdescriptions <- unique(df_rxn2subsys[,c(1,3)])
rownames(df_SubSysdescriptions) <- df_SubSysdescriptions$SubsysID
