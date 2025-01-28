library("data.table")
##Functions to summarize fastcore and post-processing results

#############################
##    Helper functions
#############################

##Load Recon 2.2
model <- readRDS("../models/recon22.RDS")

##Simplify reaction names from the metamodel into the original names
cleanReakNamesVec <- function(reaks){
    reaks <- gsub("_Lumen|_Colon|_Brain|_Liver|_DigBlood|_Blood|_back|_Microbiome|\\(c\\)|\\(e\\)","",reaks)
    reaks <- gsub("^IEX_","EX_",reaks)
    reaks
}

##Translates a list of genes to the corresponding reactions in Recon 2.2
genesToReactions <- function(genes){
    reaks <- c()
    if (sum(model@allGenes %in% genes)>=1){
        reaks <- model@react_id[apply(model@rxnGeneMat[,model@allGenes %in% genes,drop=FALSE],1,sum)>0]
        reaks <- unique(reaks)
    }
    reaks
}

##Load metabolic dependencies
##The threshold indicates the factor by which flux is decreased in germfree conditions
loadDeps <- function(threshold=10){
    depList <- list()
    for (f in list.files(path="../results/mouseModels/",".fva-mic_corrected.csv",full.name=TRUE)){
        tab <- read.table(f,header=TRUE,stringsAsFactors=FALSE,sep=",")
        name <- gsub(".fva-mic_corrected.csv|../results/mouseModels/","",f)
        dep <- tab[(-tab[,"minFluxFull"]>1e-10 | tab[,"maxFluxFull"]>1e-10) &
                   (tab[,"maxFluxGF"]-tab[,"minFluxGF"])*threshold<=tab[,"maxFluxFull"]-tab[,"minFluxFull"],"V1"]
        brain <- dep[grepl("_Brain",dep)]
        colon <- dep[grepl("_Colon",dep)]
        liver <- dep[grepl("_Liver",dep)]
        maxLength <- max(length(brain),max(length(colon),length(liver)))
        brain <- c(brain,rep("",maxLength-length(brain)))
        colon <- c(colon,rep("",maxLength-length(colon)))
        liver <- c(liver,rep("",maxLength-length(liver)))
        
        resTab <- data.frame(depBrainMicro=brain,depLiverMicro=liver,depColonMicro=colon)

        depList[[name]] <- resTab
    }

    depList
}

##Summarizes dependencies of a tissue
evalDeps <- function(deps,model){
    summaryTmp <- list()
    for (i in 1:ncol(deps[[1]]))
        summaryTmp[[colnames(deps[[1]])[i]]] <- c()

    for (dep in deps)
        for (depName in colnames(dep)){
            depList <- dep[,depName]
            depList <- depList[depList!=""]
            summaryTmp[[depName]] <- c(summaryTmp[[depName]],as.character(depList))
        }

    summary <- list()
    for (depSum in names(summaryTmp))
        summary[[depSum]] <- sort(table(summaryTmp[[depSum]]))

    summary
}

##############################
##      Analysis code
##############################

##Determine metabolic exchanges between the microbiota and individual host tissues
##cutoff indicates by which factor a flux through a reaction must be greater
##during presence of the microbiome compared to microbiome absence to be
##considered microbiome-dependent
summarizeInteractions <- function(cutoff=10){
    ##Go through all tissues
    tissues <- c("Brain","Colon","Liver")
    toHost <- list()
    toMicro <- list()
    for (tissue in tissues){
        ##Save metabolites going to host and microbiota
        toHostTmp <- c()
        toMicroTmp <- c()

        ##Load individual files
        fList <- list.files(path="../results/mouseModels",pattern=".fva-mic_corrected.csv",full.name=TRUE)
        for (f in fList){
            data <- fread(f)

            ##Filter reactions for tissue and exchanges
            data <- data[grepl(tissue,V1) & grepl("^IEX_",V1)]

            ##Determine direction of exchange
            data[,toTissue:=TRUE]
            data[!grepl("_back",V1),toTissue:=FALSE]
            
            ##Determine metabolite that is being exchanged
            data[,metabolite:=gsub("EX_|\\(e\\)|\\(c\\)","",cleanReakNamesVec(V1))]

            ##If there is forward and backward direction, only
            ##consider the one with higher flux (will happen further down,
            ##for now just sorting)
            data <- data[order(maxFluxFull,decreasing=TRUE)]
            
            ##Separate considerations for different directions

            ###################
            ##  To Host
            ###################
            ##Filter appropriate reaction directions
            sub <- data[grepl("_Lumen_Colon|_DigBlood_Liver|_Blood_Brain",V1)]
            sub <- sub[!duplicated(metabolite)]

            ##Get microbiome-dependent exchanges
            sub <- sub[maxFluxFull/maxFluxGF>=cutoff]
            toHostTmp <- c(toHostTmp,sub[toTissue==TRUE,metabolite])

            ###################
            ##  To Microbiome
            ###################
            ##Filter appropriate reaction directions
            sub <- data[grepl("_Lumen_Colon|_Blood_Liver|_Blood_Brain",V1)]
            sub <- sub[!duplicated(metabolite)]

            ##Get microbiome-dependent exchanges
            sub <- sub[maxFluxFull/maxFluxGF>=cutoff]
            toMicroTmp <- c(toMicroTmp,sub[toTissue==FALSE,metabolite])
        }

        toHost[[tissue]] <- sort(table(toHostTmp))
        toMicro[[tissue]] <- sort(table(toMicroTmp))
    }

    ##Now summarize everything into a matrix
    ##Collect all metabolites
    metas <- c()
    for (tissue in tissues)
        metas <- c(metas,names(toHost[[tissue]]),names(toMicro[[tissue]]))
    metas <- unique(metas)
    inter <- data.frame(ColonToMicro=rep(0,length(metas)),
                        MicroToColon=rep(0,length(metas)),
                        LiverToMicro=rep(0,length(metas)),
                        MicroToLiver=rep(0,length(metas)),
                        BrainToMicro=rep(0,length(metas)),
                        MicroToBrain=rep(0,length(metas)),
                        row.names=metas
                        )
    for (tissue in tissues){
        inter[names(toMicro[[tissue]]),paste0(tissue,"ToMicro")] <- toMicro[[tissue]]
        inter[names(toHost[[tissue]]),paste0("MicroTo",tissue)] <- toHost[[tissue]]
    }

    ##Some sorting
    inter <- inter[order(apply(inter,1,mean),decreasing=TRUE),]
    
    ##Save
    write.csv(inter,paste0("../results/mouseModelsSummary/inter_",cutoff,".csv"))

    inter
}

##Determine host-dependent microbiome reactions
getMicrobiomeDepsReactions <- function(cutoff=10){
    ##Get list of all reactions
    fluxes <- read.table("../results/mouseModelsSummary/fluxesSummary.csv",header=TRUE,row.names=1,sep="\t")
    allReaks <- rownames(fluxes)
    allReaks <- gsub("_LPAREN_","(",allReaks)
    allReaks <- gsub("_RPAREN_",")",allReaks)
    
    fList <- list.files(path="../results/mouseModels",pattern=".fva-mic_corrected.csv",full.name=TRUE)
    fluxes <- matrix(0,length(allReaks),length(fList))
    rownames(fluxes) <- allReaks
    colnames(fluxes) <- fList

    for (f in fList){
        ##Load microbiome data
        data <- fread(gsub("mic","host",f))
        vec <- (data[,maxFluxFull]-data[,minFluxFull])/(data[,maxFluxHF]-data[,minFluxHF])
        ##Correct numeric inaccuracies
        vec[data[,maxFluxFull]-data[,minFluxFull]<=1e-6] <- 0
        fluxes[data[,V1],f] <- vec
    }
    fluxes[fluxes==0] <- NA
    fluxes[is.infinite(fluxes)] <- NA

    fluxes <- fluxes[grepl("Microbiome",rownames(fluxes)),]
    nas <- apply(fluxes,1,function(x){sum(is.na(x))})

    ##Remove fluxes with too many NAs
    fluxes <- fluxes[nas<=ncol(fluxes)*0.5,]

    ##Determine fluxes below the cutoff
    fluxes[fluxes<cutoff] <- -1
    fluxes[fluxes>cutoff&fluxes>=0] <- 1
    fluxes[fluxes==-1] <- 0

    rownames(fluxes) <- gsub("_Microbiome","",rownames(fluxes))
    rownames(fluxes) <- gsub("^IEX_","EX_",rownames(fluxes))

    res <- sort(apply(fluxes,1,sum)/apply(fluxes,1,function(x){sum(!is.na(x))}),decreasing=TRUE)
    reakNames <- gsub("_back","",names(res))
    res <- res[!duplicated(reakNames)]
    names(res) <- reakNames[!duplicated(reakNames)]
    
    write.csv(res,paste0("../results/mouseModelsSummary/microDep_reactions_",cutoff,".csv"))
}

##determines the frequency at which each reaction is dependent on the microbiome
##cutoff: foldchange in flux upon removal of the microbiota for a reaction
##to be considered microbiome-dependent
##minTissueCount: minimal number of mice in which the reaction should be present in the
##respective tissue
getHostDepsReactions <- function(cutoff=10,minTissueCount=40){
    ##Count how often each reaction occurs in each tissue
    deps <- evalDeps(loadDeps(0.9))
    base <- matrix(0,length(model@react_id),3)
    rownames(base) <- model@react_id
    colnames(base) <- c("Brain","Colon","Liver")
    for (dep in names(deps)){
        depName <- gsub("dep|Micro","",dep)
        depList <- deps[[dep]]
        reakList <- names(depList)
        cleanedReaks <- gsub("_Lumen|_Colon|_Brain|_Liver|_DigBlood|_Blood|_back","",reakList)
        cleanedReaks <- gsub("^IEX_","EX_",cleanedReaks)
        common <- cleanedReaks %in% rownames(base)
        base[cleanedReaks[common],depName] <- depList[reakList[common]]
    }

    ##Count reaction dependencies
    deps <- evalDeps(loadDeps(cutoff))
    summary <- matrix(0,length(model@react_id),3)
    rownames(summary) <- model@react_id
    colnames(summary) <- c("Brain","Colon","Liver")
    for (dep in names(deps)){
        depName <- gsub("dep|Micro","",dep)
        depList <- deps[[dep]]
        reakList <- names(depList)
        cleanedReaks <- gsub("_Lumen|_Colon|_Brain|_Liver|_DigBlood|_Blood","",reakList)
        cleanedReaks <- gsub("^IEX_","EX_",cleanedReaks)
        common <- cleanedReaks %in% rownames(summary)
        ##common <- intersect(cleanedReaks,rownames(summary))
        summary[cleanedReaks[common],depName] <- depList[reakList[common]]
    }
    summary <- summary/base
    summary <- summary
    summary[base<minTissueCount] <- NA
    
    write.csv(summary,paste0("../results/mouseModelsSummary/hostDep_reactions_",cutoff,".csv"))
}

##determines the frequency at which each host metabolite is dependent on the microbiome
##cutoff: fold-change in maximal production of a metabolite upon removal of the microbiota
##to be considered microbiome-dependent
##minTissueCount: minimal number of mice in which the metabolite should be present in the
##respective tissue
getHostDepsMetabolites <- function(cutoff=10,tissueCutoff=40){
    compounds <- c()
    gfCompounds <- c()

    ##Get microbiome-dependent compounds for each dataset
    for (f in list.files(path="../results/mouseModels/",pattern="\\.prod\\.",full.names=TRUE)){
        data <- fread(f)
        data <- data[abs(wt)>1e-6,]
        data[,fac:=abs(wt)/abs(gf)]  ##Need to take absolutes due to numerical inaccuracies

        compounds <- c(compounds,data[,V1])
        gfCompounds <- c(gfCompounds,data[fac>=cutoff,V1])
    }
    
    ##Minimal number of required occurrences of each compound
    compounds <- table(compounds)
    allowedCompounds <- names(compounds)[compounds>=tissueCutoff]
    compounds <- compounds[compounds>=tissueCutoff]
    
    gfCompounds <- gfCompounds[gfCompounds %in% allowedCompounds]
    gfCompounds <- table(gfCompounds)

    res <- rep(0,length(compounds))
    names(res) <- names(compounds)
    res[names(gfCompounds)] <- gfCompounds/compounds[names(gfCompounds)]

    ##Remove microbiome compounds since they are not informative
    res <- res[!grepl("Microbiome_",names(res))]

    ##Summarize per tissue
    resTable <- data.table(metabolite=names(res),dep=res)
    resTable[,tissue:=gsub("_M.*","",metabolite)]
    resTable[,tissue:=gsub("Brain_|Colon_|Liver_","",tissue)]
    resTable[,metabolite:=gsub("Liver_|Brain_|Colon_|Lumen_|DigBlood_|Blood_|M_","",metabolite)]
    resTable <- resTable[order(dep,decreasing=TRUE)]
    resTable <- resTable[!duplicated(paste0(metabolite,tissue))]
    resTable <- dcast(resTable,metabolite ~ tissue,value.var="dep")
    write.csv(resTable,paste0("../results/mouseModelsSummary/hostDep_metabolites_",cutoff,".csv"))
}

##Checks for a tissue for which indicator reactions the sampled elementary flux modes are enriched for
##aging-regulated genes.
##pCutoff: Cut-off p-value for a gene to be assumed to be aging-regulated
##interCutoff: Frequency at which a reaction needs to occur in elementary flux modes of an indicator reaction
##to be assumed to be part of the pathway of that reaction
##dir: which direction of aging-regulation to consider
enrichAgingPathwaysTissues <- function(tissue="Colon",pCutoff=0.05,interCutoff=0.2,dir="down"){
    ##Load interaction matrix
    inter <- readRDS(paste0("../results/interaction/host_microbiome_interactions_",tolower(tissue),".RDS"))

    ##For mapping to human genes
    mapping <- read.table("../mapping/mappingGenes.csv",header=TRUE,sep=",",row.names=3,stringsAsFactors=FALSE)
    
    ##Generate translation table for reactions in the metamodel since they can be from several tissues
    reakMap <- data.table(index=1:(nrow(inter)+ncol(inter)),reaks=c(rownames(inter),colnames(inter)))
    reakMap[,organ:="Colon"]
    reakMap[grepl("Brain",reaks),organ:="Brain"]
    reakMap[grepl("Liver",reaks),organ:="Liver"]
    reakMap[grepl("Microbiome",reaks),organ:="Microbiome"]
    reakMap[,redName:=cleanReakNamesVec(reaks)]
    
    ##Load aging-regulated genes and translate them to reactions
    agingReaks <- c()

    tis <- tissue
    degs <- fread(paste0("../processedData/agingRegulation/DESeq",tis,"AllAges.txt"))

    genesUp <- degs[padj<=pCutoff & log2FoldChange>0,V1]
    genesDown <- degs[padj<=pCutoff & log2FoldChange<0,V1]

    reaksUp <- genesToReactions(na.omit(mapping[genesUp,"HGNC.ID"]))
    reaksDown <- genesToReactions(na.omit(mapping[genesDown,"HGNC.ID"]))

    common <- intersect(reaksUp,reaksDown)
    reaksUp <- setdiff(reaksUp,common)
    reaksDown <- setdiff(reaksDown,common)

    reaks <- reaksUp
    if (dir == "down")
        reaks <- reaksDown
    
    ##Get corresponding reactions from the metabolic model
    sel <- reakMap[,redName] %in% reaks
    agingReaks <- c(agingReaks,reakMap[reakMap[,organ]==tis & sel,reaks])

    ##Add microbiome reactions (only for counting, not for enrichment)
    ##Now load microbiome reactions (internal fluxes)
    micAging <- data.table(readRDS("../processedData/agingRegulation/df_corInternRxnFluxToAge.RDS"))
    micAging[,RxnFlux:=gsub("_c0|_e0","",RxnFlux)]
    reaks <- micAging[p.adj<=0.05,RxnFlux]
    sel <- reakMap[,redName] %in% reaks
    agingReaks <- c(agingReaks,reakMap[reakMap[,organ]=="Microbiome" & sel,reaks])
    microAgingReaks <- agingReaks[grepl("_Microbiome",agingReaks)]

    ##Microbiome exchange reactions    
    micAgingExHost <- data.table(readRDS("../processedData/agingRegulation/df_corHostMBExchangeToAge.RDS"))
    upEx <- micAgingExHost[Spearman_rho>0 & p.adj<=pCutoff,ExFlux]
    downEx <- micAgingExHost[Spearman_rho<0 & p.adj<=pCutoff,ExFlux]

    micAgingExInternal <- fread("../processedData/agingRegulation/corWithinCommunityExFluxByAge.csv")
    upEx <- c(upEx,micAgingExInternal[Spearman_rho>0 & p.adj<=pCutoff,ExFlux])
    downEx <- c(downEx,micAgingExInternal[Spearman_rho<0 & p.adj<=pCutoff,ExFlux])

    ##Translate to AGORA
    ##Translate to Recon2.2 IDs as far as possible
    mappingEx <- fread("../mapping/SEED2VMH_translation_expanded.csv",header=FALSE)
    mappingEx[,V1:=gsub("\\(e\\)|\\(c\\)|EX_","",V1)]
    mappingEx[,V2:=gsub("\\(e\\)|\\(c\\)","",V2)]

    ##Some manual adjustment to the mapping
    mappingEx <- rbind(mappingEx,data.table(V1="cpd23430",V2="EX_isochol"))
    mappingEx <- mappingEx[!duplicated(V1) & !duplicated(V2)]
    setkeyv(mappingEx,"V1")
    upEx <- gsub("EX_","",mappingEx[gsub("EX_|_e0|_c0","",upEx),V2])
    downEx <- gsub("EX_","",mappingEx[gsub("EX_|_e0|_c0","",downEx),V2])
    
    ##Get set of reactions occuring in at least five mice (for ground-truth)
    fluxes <- read.table("../results/mouseModelsSummary/individual_mousemodels.csv",header=TRUE,row.names=1,sep="\t")
    allReaks <- rownames(fluxes)[apply(fluxes,1,sum)>=5]
    allReaks <- gsub("_LPAREN_","(",allReaks)
    allReaks <- gsub("_RPAREN_",")",allReaks)
    allReaks <- allReaks[!grepl("_back",allReaks)]

    agingReaks <- intersect(agingReaks,allReaks)
    
    ##Now go through the individually sampled reactions
    inter[inter>=interCutoff] <- 1
    inter[inter<interCutoff] <- 0

    overlaps <- data.table(reaks=colnames(inter),
                           hostReakNum=0,
                           microReakNum=0,
                           brainReakNum=0,
                           colonReakNum=0,
                           liverReakNum=0,
                           overlap=0,
                           microFrac=0,
                           p.value=1)
    setkeyv(overlaps,"reaks")
    for (i in 1:ncol(inter)){
        if (i %% 100==0) print(i)
        interReaks <- rownames(inter)[inter[,i]==1]
        interReaks <- intersect(interReaks,allReaks)
        overlaps[colnames(inter)[i],hostReakNum:=sum(!grepl("Microbiome",interReaks))]
        overlaps[colnames(inter)[i],microReakNum:=sum(grepl("Microbiome",interReaks))]
        overlaps[colnames(inter)[i],brainReakNum:=sum(grepl("Brain",interReaks))]
        overlaps[colnames(inter)[i],colonReakNum:=sum(grepl("Colon",interReaks))]
        overlaps[colnames(inter)[i],liverReakNum:=sum(grepl("Liver",interReaks))]
        overlaps[colnames(inter)[i],overlap:=length(intersect(agingReaks,interReaks))/hostReakNum]
        overlaps[colnames(inter)[i],microFrac:=length(intersect(microAgingReaks,interReaks))/microReakNum]

        agingReaksTmp <- agingReaks[!grepl("Microbiome",agingReaks)]
        
        mat <- matrix(c(length(intersect(interReaks,agingReaksTmp)),
                        length(interReaks)-length(intersect(interReaks,agingReaksTmp)),
                        length(intersect(agingReaksTmp,allReaks)),
                        length(allReaks)-length(intersect(agingReaksTmp,allReaks))
                        ),2,2)
        overlaps[colnames(inter)[i],p.value:=fisher.test(t(mat))$p.value]
    }

    ##Filter significantly enriched pathways
    overlaps[,p.adj:=p.adjust(p.value,method="fdr")]
    overlaps <- overlaps[p.adj<=0.05]
    overlaps <- overlaps[order(p.adj)]
    write.csv(overlaps,paste0("../results/agingPathwayEnrichment/agingEnrichment_",tissue,"_",dir,".csv"))
    overlaps
}

##Identifies target reactions whose elementary flux modes are enriched for aging-regulated reactions
enrichAgingPathways <- function(){
    ##Enrich aging pathways for all tissues
    enrichAgingPathwaysTissues(tissue="Colon",pCutoff=0.05,interCutoff=0.2,dir="down")
    enrichAgingPathwaysTissues(tissue="Colon",pCutoff=0.05,interCutoff=0.2,dir="up")

    enrichAgingPathwaysTissues(tissue="Liver",pCutoff=0.05,interCutoff=0.2,dir="down")
    enrichAgingPathwaysTissues(tissue="Liver",pCutoff=0.05,interCutoff=0.2,dir="up")

    enrichAgingPathwaysTissues(tissue="Brain",pCutoff=0.05,interCutoff=0.2,dir="down")
    enrichAgingPathwaysTissues(tissue="Brain",pCutoff=0.05,interCutoff=0.2,dir="up")

    ##Filter for shared reactions between up and down for each tissue
    for (tissue in c("Brain","Colon","Liver")){
        up <- fread(paste0("../results/agingPathwayEnrichment/agingEnrichment_",tissue,"_up.csv"))
        down <- fread(paste0("../results/agingPathwayEnrichment/agingEnrichment_",tissue,"_down.csv"))

        common <- intersect(up[,reaks],down[,reaks])
        up <- up[!(reaks %in% common)]
        down <- down[!(reaks %in% common)]

        write.csv(up,paste0("../results/agingPathwayEnrichment/agingEnrichment_",tissue,"_up.csv"),row.names=FALSE)
        write.csv(down,paste0("../results/agingPathwayEnrichment/agingEnrichment_",tissue,"_down.csv"),row.names=FALSE)
    }

}
