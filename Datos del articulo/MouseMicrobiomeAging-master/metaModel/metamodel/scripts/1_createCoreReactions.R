library("data.table")
##library("vegan")
##library("ggplot2")
##library("ggrepel")
library("sva")
library("preprocessCore")
library("readxl")

## Script for creating core reactions from transcriptomic data of the individual tissues
## and metagenomic data

##Entire process

##Normalize data
##normalizeData()

##Prepare data matrices for StanDep
##processRNASeqDataStanDep()

##run StanDep script (xxx.m)

##Pull together all data and build core reaction matrix file
##combineData()


##Generate a matrix indicating the present reactions for each microbial model
getContribMatrixGapseq <- function(){
    gapseq <- readRDS("../models/metamouse-20230510.RDS")
    allReaks <- c()
    for (model in gapseq){
        allReaks <- union(allReaks,model@react_id)
    }
    contrib <- matrix(0,length(allReaks),length(gapseq))
    rownames(contrib) <- allReaks
    colnames(contrib) <- names(gapseq)
    for (modelName in names(gapseq)){
        model <- gapseq[[modelName]]
        contrib[model@react_id,modelName] <- 1
    }

    contrib
}

mapGenes <- function(data,geneMapping){
    data <- data[intersect(rownames(data),rownames(geneMapping)),]

    ##Now append missing genes
    missingGenes <- unique(setdiff(rownames(geneMapping),rownames(data)))
    missingMat <- matrix(0,length(missingGenes),ncol(data))
    rownames(missingMat) <- missingGenes
    colnames(missingMat) <- colnames(data)
    data <- rbind(data,missingMat)
    
    data
}

##Loads the individual expression data files and combines them
##Normalizes data using ComBat to remove batch effects
normalizeData <- function(){fs <- list() 
    sampleInfo <- readRDS("../rawData/df_transcriptMetaData.rds")
    
    #####################################
    ##    Process brain data (ComBat)
    #######################################
    brain <- readRDS("../rawData/mtx_brainRNAtpm.rds")
    subMeta <- sampleInfo[sampleInfo[,"location"]=="Brain",]
    rownames(subMeta) <- subMeta[,"id"]
    subMeta <- subMeta[colnames(brain),]
    
    means <- apply(brain,1,mean)
    zeros <- apply(brain,1,function(x){sum(x<=0.01)})
    brain <- brain[means>0.1 & zeros==0,]

    brainNorm <- log2(brain)
    rownames(brainNorm) <- rownames(brain)
    colnames(brainNorm) <- colnames(brain)
    brainNorm <- ComBat(brainNorm,as.character(subMeta[,"run"]),mod=model.matrix(~age,data=subMeta),par.prior=TRUE,mean.only=FALSE)  ##PERFECT
    write.csv(2^brainNorm,"../processedData/brain.csv")

    #######################################
    ##    Process liver data (ComBat)
    ######################################
    liver <- readRDS("../rawData/mtx_liverRNAtpm.rds")
    subMeta <- sampleInfo[sampleInfo[,"location"]=="Liver",]
    rownames(subMeta) <- subMeta[,"id"]
    subMeta <- subMeta[colnames(liver),]
    
    means <- apply(liver,1,mean)
    zeros <- apply(liver,1,function(x){sum(x<=0.01)})
    liver <- liver[means>0.1 & zeros==0,]

    liverNorm <- log2(liver)
    rownames(liverNorm) <- rownames(liver)
    colnames(liverNorm) <- colnames(liver)
    liverNorm <- ComBat(liverNorm,as.character(subMeta[,"run"]),mod=model.matrix(~age,data=subMeta),par.prior=TRUE,mean.only=FALSE)  ##PERFECT
    write.csv(2^liverNorm,"../processedData/liver.csv")

    #######################################
    ##    Process colon data (ComBat) (no COMBAT required since data was generated in one run)
    ######################################
    colon <- readRDS("../rawData/mtx_colonRNAtpm.rds")
    subMeta <- sampleInfo[sampleInfo[,"location"]=="Colon",]
    rownames(subMeta) <- subMeta[,"id"]
    subMeta <- subMeta[colnames(colon),]
    
    means <- apply(colon,1,mean)
    zeros <- apply(colon,1,function(x){sum(x<=0.01)})
    colon <- colon[means>0.1 & zeros==0,]

    colonNorm <- colon
    write.csv(colonNorm,"../processedData/colon.csv")
}

getMetagenomeDataGapseq <- function(metagenomeMapping="../rawData/df_MAGabundances_20210907.txt"){
    contrib <- getContribMatrixGapseq()

    meta <- read.table(metagenomeMapping,header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=1,check.names=FALSE)
    ##meta <- meta[,-1]
    ##Adjust IDs
    colnames(meta) <- gsub("_S.*","",colnames(meta))
    
    ##Normalize contribs
    contrib <- prop.table(contrib,2)
    common <- intersect(rownames(meta),colnames(contrib))
    meta <- meta[common,]
    
    ##Normalize to 1
    reaks <- t(t(meta[common,]) %*% t(contrib[,common]))
    reaks <- prop.table(reaks,2)

    reaks
}

##Appends information about net release and production of metabolites by individual organs
addExchangeData <- function(data){
    ##Load exchanges and map to Agora
    exchanges <- data.table(read_excel("../rawData/organ-exchange/exchanges.xlsx","Sheet1"))
    mapping <- data.table(read_excel("../rawData/organ-exchange/metabolite-matching.xlsx","Matching_final"))
    mapping <- mapping[Agora_name!="",]
    setkeyv(mapping,"Name")
    exchanges[,"metabolite"] <- mapping[exchanges[,metabolite],Agora_ID]
    exchanges <- exchanges[!is.na(exchanges[,metabolite])]

    ##Now add exchange information for the exchanges of each organ
    ##"1" indicates that an exchange was observed in a specific direction (uptake or secretion)
    ##and hence the reaction is included in the core list for
    ##fastcore. The opposite direction is indicated by "-1" which can not used by fastcore,
    ##but other approaches for mapping expression data such as iMAT

    ##Brain
    print("Brain exchange")
    brainUptake <- exchanges[head==-1,metabolite] ##From blood
    brainRelease <- exchanges[head==1,metabolite] ##Into blood
    for (meta in brainUptake){
        if (paste0("IEX_",meta,"(e)_Blood_Brain") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_Blood_Brain"),] <- -1
        if (paste0("IEX_",meta,"(e)_Blood_Brain_back") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_Blood_Brain_back"),] <- 1
    }
    for (meta in brainRelease){
        if (paste0("IEX_",meta,"(e)_Blood_Brain") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_Blood_Brain"),] <- 1
        if (paste0("IEX_",meta,"(e)_Blood_Brain_back") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_Blood_Brain_back"),] <- -1
    }

    ##Colon
    print("Colon exchange")
    colonUptake <- exchanges[colon==-1,metabolite] ##Via the blood
    colonRelease <- exchanges[colon==1,metabolite] ##Into DigBlood
    for (meta in colonUptake){
        if (paste0("IEX_",meta,"(e)_DigBlood_Colon") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_DigBlood_Colon"),] <- -1
        if (paste0("IEX_",meta,"(e)_Blood_Colon_back") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_Blood_Colon_back"),] <- 1
    }

    for (meta in colonRelease){
        if (paste0("IEX_",meta,"(e)_DigBlood_Colon") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_DigBlood_Colon"),] <- 1
        if (paste0("IEX_",meta,"(e)_Blood_Colon_back") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_Blood_Colon_back"),] <- -1
    }

    ##Liver
    print("Liver exchange")
    liverUptake <- exchanges[liver==-1,metabolite] ##Via DigBlood and Blood
    liverRelease <- exchanges[liver==1,metabolite] ##Into Blood
    for (meta in liverUptake){
        if (paste0("IEX_",meta,"(e)_DigBlood_Liver") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_DigBlood_Liver"),] <- -1
        if (paste0("IEX_",meta,"(e)_DigBlood_Liver_back") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_DigBlood_Liver_back"),] <- 1
        if (paste0("IEX_",meta,"(e)_Blood_Liver_back") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_Blood_Liver_back"),] <- 1
    }
    for (meta in liverRelease){
        if (paste0("IEX_",meta,"(e)_Blood_Liver") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_Blood_Liver"),] <- 1
        if (paste0("IEX_",meta,"(e)_DigBlood_Liver_back") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_DigBlood_Liver_back"),] <- -1
        if (paste0("IEX_",meta,"(e)_DigBlood_Liver_back") %in% rownames(data))
            data[paste0("IEX_",meta,"(e)_Blood_Liver_back"),] <- -1
    }

    ##Outflow through the kidneys
    print("Kidney exchange")
    kidneyUptake <- exchanges[kidney==-1,metabolite] ##Via the blood
    kidneyRelease <- exchanges[kidney==1,metabolite] ##Try to avoid excretion from blood
    for (meta in kidneyUptake){
        if (paste0("EX_",meta,"(e)_outflow") %in% rownames(data))
            data[paste0("EX_",meta,"(e)_outflow"),] <- 1
    }
    for (meta in kidneyRelease){
        if (paste0("EX_",meta,"(e)_outflow") %in% rownames(data))
            data[paste0("EX_",meta,"(e)_outflow"),] <- -1
    }
    
    data
}

##Performs some translations for the microbiome reactions
translateMicrobiomeReaks <- function(reaks){
    mapping <- fread("../mapping/SEED2VMH_translation_expanded.csv",header=FALSE)
    mapping[,V1:=gsub("EX_|\\(e\\)|\\(c\\)","",V1)]
    mapping[,V2:=gsub("EX_|\\(e\\)|\\(c\\)","",V2)]
    setkeyv(mapping,"V1")
    
    adaptedReaks <- gsub("_c0_","(c)_",reaks)
    adaptedReaks <- gsub("_e0_","(e)_",adaptedReaks)

    ##Now some adaptations for exchange compounds
    adaptedReaks <- gsub("_EX_","_IEX_",adaptedReaks)

    for (i in 1:length(adaptedReaks))
        if(grepl("_IEX_",adaptedReaks[i])){
            cpd <- gsub("G_R_IEX_|\\(.*","",adaptedReaks[i])

            if (cpd %in% mapping[,V1]){
                adaptedReaks[i] <- gsub(cpd,mapping[cpd,V2],adaptedReaks[i])
            }
        }
    adaptedReaks
}

processRNASeqDataStanDep <- function(){
    ##map sample ids
    metadata <- readRDS("../rawData/df_transcriptMetaData.rds")
    metadata[,"id_new"] <- paste(metadata[,"age"],gsub("S_","",metadata[,"id"]),sep="_")

    ##Only need data for the mice, not individual tissues
    metadata <- metadata[metadata[,"location"]=="Colon",]
    rownames(metadata) <- metadata[,"id"]
    
    ##Brain
    brain <- read.table("../processedData/brain.csv",sep=",",header=TRUE,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
    rownames(brain) <- gsub("\\.\\d+","",rownames(brain))
    colnames(brain) <- metadata[colnames(brain),"id_new"]
    
    ##Liver
    liver <- read.table("../processedData/liver.csv",sep=",",header=TRUE,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
    rownames(liver) <- gsub("\\.\\d+","",rownames(liver))
    colnames(liver) <- metadata[colnames(liver),"id_new"]
    
    ##Colon
    colon <- read.table("../processedData/colon.csv",sep=",",header=TRUE,row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
    rownames(colon) <- gsub("\\.\\d+","",rownames(colon))
    colnames(colon) <- metadata[colnames(colon),"id_new"]
    
    geneMapping <- read.table("../mapping/mappingGenes.csv",header=TRUE,row.names=3,sep=",",stringsAsFactors=FALSE)

    ##Determine samples for which all data is available
    commonSamples <- intersect(intersect(colnames(colon),colnames(liver)),colnames(brain))
    
    colon <- mapGenes(colon[,commonSamples],geneMapping)
    liver <- mapGenes(liver[,commonSamples],geneMapping)
    brain <- mapGenes(brain[,commonSamples],geneMapping)

    ##Adjust gene names to the model
    rownames(colon) <- geneMapping[rownames(colon),1]
    rownames(liver) <- geneMapping[rownames(liver),1]
    rownames(brain) <- geneMapping[rownames(brain),1]

    genes <- intersect(intersect(rownames(colon),rownames(liver)),rownames(brain))

    colon <- colon[,commonSamples]
    brain <- brain[,commonSamples]
    liver <- liver[,commonSamples]
    
    colnames(colon) <- paste0("colon_",colnames(colon))
    colnames(brain) <- paste0("brain_",colnames(brain))
    colnames(liver) <- paste0("liver_",colnames(liver))
    
    fullData <- cbind(colon[genes,],liver[genes,],brain[genes,])
    fullData[fullData==0] <- NA
   
    fullData <- fullData[apply(fullData,1,function(x){sum(is.na(x))})<ncol(fullData),]
    fullData[is.na(fullData)] <- 0
    
    write.table(fullData,"../processedData/standep/expression-host.tsv",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    write.table(gsub("_\\d+$","",colnames(fullData)),"../processedData/standep/conditions-host.tsv",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    write.table(colnames(fullData),"../processedData/standep/sampleNames-host.tsv",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    write.table(rownames(fullData),"../processedData/standep/genes-host.tsv",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    
    ###################################
    ##   Process metagenome data
    ###################################
    meta <- getMetagenomeDataGapseq()
    rownames(meta) <- paste0("G_R_",rownames(meta),"_Microbiome")
    
    ##Divide by 100 to get into normal range of TPMs/FPKMs
    meta <- meta[,commonSamples]
    tmp <- normalize.quantiles(as.matrix(meta))##/100
    rownames(tmp) <- rownames(meta)
    colnames(tmp) <- colnames(meta)
    meta <- tmp

    ##Remove biomass exchange (from expression data)
    meta <- meta[!grepl("bio1|cpd11416",rownames(meta)),]
    
    adaptedReaks <- gsub("_c0_","(c)_",rownames(meta))
    adaptedReaks <- gsub("_e0_","(e)_",adaptedReaks)
    adaptedReaks <- gsub("_EX_","_IEX_",adaptedReaks)

    ##Write data
    write.table(meta,"../processedData/standep/expression-microbiome.tsv",row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    write.table(paste0("microbiome_",gsub("_\\d+$","",colnames(meta))),
                "../processedData/standep/conditions-microbiome.tsv",
                row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    write.table(translateMicrobiomeReaks(rownames(meta)),
                "../processedData/standep/genes-microbiome.tsv",
                row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    write.table(paste0("microbiome_",colnames(meta)),
                "../processedData/standep/samples-microbiome.tsv",
                row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
}

addMicrobiomeData <- function(data){
    microData <- fread("../processedData/standep/coreRxns-microbiome.csv")
    sampleNames <- fread("../processedData/standep/samples-microbiome.tsv",header=FALSE)[,V1]
    colnames(microData) <- gsub("microbiome_","",sampleNames)
    reaks <- fread("../processedData/standep/reactions-metamodel.csv")[,Var1]
    microData[,reaks:=reaks]

    for (sample in colnames(data)){
        reaks <- microData[microData[[sample]]==1,reaks]
        data[reaks,sample] <- 1
    }
    data
}

combineData <- function(){
    ##Derive matrix with core reactions for each metamodel
    coreRxns <- fread("../processedData/standep/coreRxns-host.csv")
    sampleNames <- fread("../processedData/standep/sampleNames-host.tsv",header=FALSE)[,V1]
    colnames(coreRxns) <- sampleNames
    coreRxns[,reaks:=fread("../processedData/standep/reactions-host.csv")]
    samples <- unique(gsub("colon_|brain_|liver_","",sampleNames))
    
    ##Convert into matrix for metamodel
    metamodelReaks <- fread("../processedData/standep/reactions-metamodel.csv")[,Var1]
    coreReakMat <- matrix(0,length(metamodelReaks),length(samples))
    colnames(coreReakMat) <- samples
    rownames(coreReakMat) <- metamodelReaks
    
    ##Collect data from all samples
    for (sample in samples){
        ## Process core reactions
        ##Colon
        colon <- coreRxns[coreRxns[[paste0("colon_",sample)]]==1,reaks]
        colon <- paste0(colon,"_Colon")
        
        ##Transport reactions are named differently
        rem <- setdiff(colon,metamodelReaks)
        colon <- setdiff(colon,rem)
        set1 <- gsub("_Colon","_Blood_Colon",rem)
        set2 <- gsub("_Colon","_DigBlood_Colon",rem)
        colon <- union(union(colon,set1),set2)
        coreReakMat[colon,sample] <- 1

        ##Liver
        liver <- coreRxns[coreRxns[[paste0("liver_",sample)]]==1,reaks]
        liver <- paste0(liver,"_Liver")

        ##Transport reactions are named differently
        rem <- setdiff(liver,metamodelReaks)
        liver <- setdiff(liver,rem)
        set1 <- gsub("_Liver","_Blood_Liver",rem)
        set2 <- gsub("_Liver","_DigBlood_Liver",rem)
        liver <- union(union(liver,set1),set2)
        coreReakMat[liver,sample] <- 1

        ##Brain
        brain <- coreRxns[coreRxns[[paste0("brain_",sample)]]==1,reaks]
        brain <- paste0(brain,"_Brain")

        ##Transport reactions are named differently
        rem <- setdiff(brain,metamodelReaks)
        brain <- setdiff(brain,rem)
        set1 <- gsub("_Brain","_Blood_Brain",rem)
        brain <- union(brain,set1)
        coreReakMat[brain,sample] <- 1

        ## Process non-core reactions
        ##Colon
        colon <- coreRxns[coreRxns[[paste0("colon_",sample)]]==-1,reaks]
        colon <- paste0(colon,"_Colon")
        
        ##Transport reactions are named differently
        rem <- setdiff(colon,metamodelReaks)
        colon <- setdiff(colon,rem)
        set1 <- gsub("_Colon","_Blood_Colon",rem)
        set2 <- gsub("_Colon","_DigBlood_Colon",rem)
        colon <- union(union(colon,set1),set2)
        coreReakMat[colon,sample] <- -1

        ##Liver
        liver <- coreRxns[coreRxns[[paste0("liver_",sample)]]==-1,reaks]
        liver <- paste0(liver,"_Liver")

        ##Transport reactions are named differently
        rem <- setdiff(liver,metamodelReaks)
        liver <- setdiff(liver,rem)
        set1 <- gsub("_Liver","_Blood_Liver",rem)
        set2 <- gsub("_Liver","_DigBlood_Liver",rem)
        liver <- union(union(liver,set1),set2)
        coreReakMat[liver,sample] <- -1

        ##Brain
        brain <- coreRxns[coreRxns[[paste0("brain_",sample)]]==-1,reaks]
        brain <- paste0(brain,"_Brain")

        ##Transport reactions are named differently
        rem <- setdiff(brain,metamodelReaks)
        brain <- setdiff(brain,rem)
        set1 <- gsub("_Brain","_Blood_Brain",rem)
        brain <- union(brain,set1)
        coreReakMat[brain,sample] <- -1
    }

    ##Add exchange reactions from metabolomic data
    coreReakMat <- addExchangeData(coreReakMat)
    coreReakMat <- data.frame(addMicrobiomeData(coreReakMat))
    colnames(coreReakMat) <- gsub("^X","",colnames(coreReakMat))

    ##Properly add reaction names
    coreReakMat$reaction_id <- rownames(coreReakMat)

    coreReakMat <- coreReakMat[,c("reaction_id",samples)]
    write.table(coreReakMat,"../processedData/coreReakMat.tsv",row.names=FALSE,sep="\t",quote=FALSE)
}
