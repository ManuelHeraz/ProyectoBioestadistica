library("sybil")
library("cplexAPI")
library("sybilSBML")
library("data.table")
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
sybil::SYBIL_SETTINGS("SOLVER_CTRL_PARM",data.frame("CPX_PARAM_THREADS"=1L,"CPXPARAM_TimeLimit"=100000,"CPX_PARAM_LPMETHOD"=1L)); ok <- 1

##Functions for post-processing and further analysis of the generated metamodels for each mouse
##Order of execution:
##
## ## Correct FVA results for exchange reactions
## correctExchanges()
## combineFVAResults()
##
## ## Test production capacity of metabolites
## testMetaboliteProductionHost()
## testMetaboliteProductionMicro()

##In the basic FVA run after fastcore exchange reactions are not considered separately even though they have
##been split into irreversible forward and backward direction. Here FVA for exchange reactions is run again
##for each direction of a exchange reaction while blocking the other direction
correctExchanges <- function(){
    ##Load list of active reactions for each metamodel
    fluxes <- fread("../results/mouseModelsSummary/fluxesSummary.csv")

    ##Test the capacity of exchange reactions but always block the reverse direction

    ##Determine active reactions
    for (modInd in 2:ncol(fluxes))
        if(!file.exists(paste0("../results/mouseModels/",colnames(fluxes)[modInd],".exchange.csv"))){
            print(modInd)
            model <- readRDS("../models/metamodel_20230510.RDS")
            
            mod <- colnames(fluxes)[modInd]
            reaks <- fluxes[fluxes[[modInd]]!=0,Reak]
            reaks <- gsub("_LPAREN_","(",reaks)
            reaks <- gsub("_RPAREN_",")",reaks)
            
            model@lowbnd[!(model@react_id %in% reaks)] <- 0
            model@uppbnd[!(model@react_id %in% reaks)] <- 0
            
            ##Germfree variant
            modelGF <- model
            modelGF@lowbnd[grepl("_Microbiome",model@react_id)] <- 0
            modelGF@uppbnd[grepl("_Microbiome",model@react_id)] <- 0
            
            ##Hostfree variant
            modelHF <- model
            modelHF@lowbnd[grepl("_Colon|_Liver|_Brain",model@react_id)] <- 0
            modelHF@uppbnd[grepl("_Colon|_Liver|_Brain",model@react_id)] <- 0
            
            ##Determine associated reactions
            ##exReaks <- model@react_id[grepl("^IEX_",model@react_id) & grepl("Microbiome",model@react_id)]
            exReaks <- model@react_id[grepl("^IEX_",model@react_id)]
            exReaks <- intersect(reaks,exReaks)
            
            res <- matrix(0,length(exReaks),2) ##No minimal flux since exchanges reactions are irreversible
            write.csv(res,paste0("../results/mouseModels/",mod,".exchange.csv"))
            
            rownames(res) <- exReaks
            colnames(res) <- c("maxWT","maxGF_HF")
            
            for (i in 1:length(exReaks)){
                reak <- exReaks[i]
                testModel <- model
                
                ##Block reverse directions
                reverse <- ""
                if (grepl("_back",reak)){
                    reverse <- gsub("_back","",reak)
                }else{
                    reverse <- paste0(reak,"_back")
                }
                                          
                testModel@uppbnd[reverse==testModel@react_id] <- 0
                
                testModel@obj_coef[] <- 0
                testModel@obj_coef[reak==model@react_id] <- 1
                
                opt <- optimizeProb(testModel)
                if (opt@lp_stat==1)
                    res[i,1] <- opt@lp_obj
                
                ##print(paste(meta,paste(res[meta,],collapse=" ")))
                
                ##Now without microbiome or host
                woModel <- modelGF
                if (grepl("_Microbiome",reak))
                    woModel <- modelHF
                
                woModel@uppbnd[reverse==woModel@react_id] <- 0
                
                woModel@obj_coef[] <- 0
                woModel@obj_coef[reak==woModel@react_id] <- 1
                                          
                opt <- optimizeProb(woModel)
                if (opt@lp_stat==1)
                    res[i,"maxGF_HF"] <- opt@lp_obj
                                          
                ##Write intermediate results
                if (i%%10==0)
                    write.csv(res,paste0("../results/mouseModels/",mod,".exchange.csv"))
            }
            
            ##Write full file
            write.csv(res,paste0("../results/mouseModels/",mod,".exchange.csv"))
        }
}

##Create corrected FVA results
combineFVAResults <- function(){
    files <- list.files(path="../results/mouseModels",pattern=".exchange.csv",full.names=TRUE)
    for (f in files){
        fvaHostFile <- gsub(".exchange.",".fva-mic.",f)
        fvaMicFile <- gsub(".exchange.",".fva-host.",f)

        ##Load corrected exchange data
        exchanges <- fread(f)
        setkeyv(exchanges,"V1")

        ##Correct host side (lower bound is zero -> need only to correct maximum)
        fvaHost <- fread(fvaHostFile)
        fvaHost[,V1:=gsub("_LPAREN_","(",V1)]
        fvaHost[,V1:=gsub("_RPAREN_",")",V1)]
        setkeyv(fvaHost,"V1")
        fvaHost[exchanges[,V1],"maxFluxFull"] <- exchanges[,maxWT]
        fvaHost[exchanges[,V1],"maxFluxGF"] <- exchanges[,maxGF_HF]
        write.table(fvaHost,gsub(".fva-mic.",".fva-mic_corrected.",fvaHostFile),
                    sep=",",quote=FALSE,row.names=FALSE)
        
        ##Correct microbiome side
        fvaMic <- fread(fvaMicFile)
        fvaMic[,V1:=gsub("_LPAREN_","(",V1)]
        fvaMic[,V1:=gsub("_RPAREN_",")",V1)]
        setkeyv(fvaMic,"V1")
        fvaMic[exchanges[,V1],"maxFluxFull"] <- exchanges[,maxWT]
        fvaMic[exchanges[,V1],"maxFluxHF"] <- exchanges[,maxGF_HF]

        fvaMic[V1 %in% exchanges[,V1],maxFluxFull:=exchanges[V1,maxWT]]
        fvaMic[V1 %in% exchanges[,V1],maxFluxHF:=exchanges[V1,maxGF_HF]]
        write.table(fvaMic,gsub(".fva-host.",".fva-host_corrected.",fvaMicFile),
                    sep=",",quote=FALSE,row.names=FALSE)
    }
}

##Determines production capacity for each host metabolite with and without the microbiota
testMetaboliteProductionHost <- function(){
    ##Load list of active reactions
    fluxes <- fread("../results/mouseModelsSummary/fluxesSummary.csv")

    ##Determine active reactions
    for (modInd in 2:ncol(fluxes))
        if(!file.exists(paste0("../results/mouseModels/",colnames(fluxes)[modInd],".prod.csv"))){
            print(modInd)
            model <- readRDS("../models/metamodel_20230510.RDS")
            
            mod <- colnames(fluxes)[modInd]
            reaks <- fluxes[fluxes[[modInd]]!=0,Reak]
            reaks <- gsub("_LPAREN_","(",reaks)
            reaks <- gsub("_RPAREN_",")",reaks)
            
            model@lowbnd[!(model@react_id %in% reaks)] <- 0
            model@uppbnd[!(model@react_id %in% reaks)] <- 0
            
            ##Determine associated metabolites
            metas <- model@met_id[apply(abs(model@S[,model@react_id %in% reaks]),1,max)>0]
            
            micro <- which(grepl("_Microbiome",model@react_id))
            
            res <- matrix(0,length(metas),2)
            
            rownames(res) <- metas
            colnames(res) <- c("wt","gf")

            ##Create a temporary file
            write.csv(res,paste0("../results/mouseModels/",mod,".prod.csv"))
            
            for (i in 1:length(metas)){##Go through all metabolites
                meta <- metas[i]

                ##Add a reaction consuming that metabolite
                testModel <- addReact(model,met=meta,Scoef=c(-1),id="testReact")
                testModel@obj_coef[] <- 0
                testModel@obj_coef[testModel@react_num] <- 1

                ##Maximize outflow of metabolite = maximal production capacity
                opt <- optimizeProb(testModel)
                if (opt@lp_stat==1)
                    res[i,1] <- opt@lp_obj

                ##Maximal production without the microbiome
                testModel@lowbnd[micro] <- 0
                testModel@uppbnd[micro] <- 0
                opt <- optimizeProb(testModel)
                if (opt@lp_stat==1)
                    res[i,2] <- opt@lp_obj
                
                ##Write intermediate results
                if (i%%1000==0)
                    write.csv(res,paste0("../results/mouseModels/",mod,".prod.csv"))
            }
            
            ##Write full file
            write.csv(res,paste0("../results/mouseModels/",mod,".prod.csv"))
        }
}

##Determines production capacity for each microbial metabolite with and without the host
testMetaboliteProductionMicro <- function(){
    ##Load list of active reactions
    fluxes <- fread("../results/mouseModelsSummary/fluxesSummary.csv")

    ##Determine active reactions
    for (modInd in 2:ncol(fluxes))
        if(!file.exists(paste0("../results/mouseModels/",colnames(fluxes)[modInd],".mic-prod.csv"))){
            print(modInd)
            model <- readRDS("../models/metamodel_20230510.RDS")
            
            mod <- colnames(fluxes)[modInd]
            reaks <- fluxes[fluxes[[modInd]]!=0,Reak]
            reaks <- gsub("_LPAREN_","(",reaks)
            reaks <- gsub("_RPAREN_",")",reaks)

            model@lowbnd[!(model@react_id %in% reaks)] <- 0
            model@uppbnd[!(model@react_id %in% reaks)] <- 0

            ##Determine associated metabolites
            metas <- model@met_id[apply(abs(model@S[,model@react_id %in% reaks]),1,max)>0]

            ##Only microbial metabolites are of interest
            metas <- metas[grepl("Microbiome",metas)]
            
            ##metas <- c("Liver_M_nadh[c]","Liver_M_adp[c]","Liver_M_amp[c]")
            host <- which(grepl("_Liver|_Colon|_Brain",model@react_id))

            res <- matrix(0,length(metas),2)
            write.csv(res,paste0("../results/mouseModels/",mod,".mic-prod.csv"))
            
            rownames(res) <- metas
            colnames(res) <- c("wt","gf")
            for (i in 1:length(metas)){##Go through all microbial metabolites
                meta <- metas[i]

                ##Add outflow reaction for metabolite
                testModel <- addReact(model,met=meta,Scoef=c(-1),id="testReact")
                testModel@obj_coef[] <- 0
                testModel@obj_coef[testModel@react_num] <- 1

                ##Optimize outflow
                opt <- optimizeProb(testModel)
                if (opt@lp_stat==1)
                    res[i,1] <- opt@lp_obj

                ##Now without the microbiome
                testModel@lowbnd[host] <- 0
                testModel@uppbnd[host] <- 0
                opt <- optimizeProb(testModel)
                if (opt@lp_stat==1)
                    res[i,2] <- opt@lp_obj

                ##Write intermediate results
                if (i%%10==0)
                    write.csv(res,paste0("../results/mouseModels/",mod,".mic-prod.csv"))
            }

            ##Write full file
            write.csv(res,paste0("../results/mouseModels/",mod,".mic-prod.csv"))
        }
}
