###############################################
###
###     Rxn-ID enrichment to subsystems
###
###############################################

source("customFunctions.R")
setwd("databases/")


################
# function definition

## subsystem annotations from gapseq models by CK Feb2022
# Load mapping from "pathway IDs" to "clear names"
df_ModelSubsysNameMapping <- read.table(file = "./gapseqPathwayIDsToNames.csv",
                                        header = T, stringsAsFactors = F, sep = "\t", quote = "", fill = TRUE)
void_extractSubsys <- function(lst_obj_models, str_fileName = "subsysMap"){
  mapping <- df_ModelSubsysNameMapping
  colnames(mapping)[1:2] <- c("ID","Name")
  rownames(mapping) <- mapping[,"ID"]
  # setkeyv(mapping,"ID")
  
  ##Obtain list of reactions
  allReaks <- c()
  for (model in lst_obj_models){
    allReaks <- union(allReaks,model@react_id)
  }
  
  ##Create lists with subsystems
  subsysMap <- list()
  for (reak in allReaks)
    subsysMap[[reak]] <- list()
  
  count <- 1
  for (model in lst_obj_models){
    print(count)
    subsys <- apply(model@subSys,1,function(x){colnames(model@subSys)[x]})
    for (i in 1:length(subsys)){
      subsysMap[[model@react_id[i]]] <- union(subsysMap[[model@react_id[i]]],gsub("\\|","",subsys[[i]]))
    }
    count <- count + 1
  }
  
  ##Only retain pathways with minimal number of reactions
  reakCount <- table(unlist(subsysMap))
  keepPathways <- names(reakCount)[reakCount>2]
  
  finalMap <- list()
  for (reak in names(subsysMap)){ 
    subsysList <- intersect(subsysMap[[reak]],keepPathways)
    subsysList <- c(subsysList[!(subsysList %in% mapping[,"ID"])],
                    na.omit(mapping[mapping[,"ID"] %in% subsysList,"Name"]))
    finalMap[[reak]] <- subsysList
  }
  
  saveRDS(finalMap, paste0(str_fileName, ".rds"))
}


################
## Create subsystem annotation by extracting reactions and corresponding subsys from models
obj_metabolicModelsMouse202305 <- readRDS("../databases/obj_metamouse-2023-05-10.rds")
void_extractSubsys(lst_obj_models = obj_metabolicModelsMouse202305, str_fileName = "rxns2subsysMetaMouse20230510")
# 181 models


###############
## Prepare subsystem annotation for cluster profiler
lst_rxn2subsys <- readRDS("rxns2subsysMetaMouse20230510.rds")
# convert list to data frame
df_rxn2subsys <- data.frame(ReactionID = rep(names(lst_rxn2subsys), times = unlist(lapply(X = lst_rxn2subsys, FUN = length))),
                            Subsystem = unlist(lst_rxn2subsys), row.names = NULL)
# Revert HTML-special chars
df_rxn2subsys$Subsystem <- textutils::HTMLdecode(df_rxn2subsys$Subsystem)
# remove html tags
df_rxn2subsys$Subsystem <- gsub("<.*?>", "", df_rxn2subsys$Subsystem)
# discard non-bacteria stuff
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "engineered", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "yeast", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "fungi", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "animal", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "plant", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "eukaryot", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "mammal", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "invertebrate", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "Plasmodium", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "archae", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "Wollastonia", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "Spartina", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "Firefly", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
# drop plant specific pathways
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "photo", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "humulone", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "tomatine", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "linamarin", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "coumarin", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "lotaustralin", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "linustatin", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "neolinustatin", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "rosmarinic", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "suberin", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "sporopollenin", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "cutin", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "estolide", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "genanyl", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "chitin", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "wood", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]
df_rxn2subsys <- df_rxn2subsys[unique(c(grep(pattern = "geranial", x = df_rxn2subsys$Subsystem, invert = T, ignore.case = T),
                                        grep(pattern = "bacteria", x = df_rxn2subsys$Subsystem))),]

# Add subsystem ID
df_rxn2subsys$SubsysID <- paste0("Subsys:",as.numeric(as.factor(df_rxn2subsys$Subsystem)))
# relabel certain subsystem names
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY66-430","Subsystem"] <- "fatty acid biosynthesis type II"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8148","Subsystem"] <- "NADP+ biosynthesis"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8130","Subsystem"] <- "5'-deoxyadenosine degradation I"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8142","Subsystem"] <- "L-ascorbate biosynthesis VI"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-5966-1","Subsystem"] <- "fatty acid biosynthesis initiation I"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8175","Subsystem"] <- "iso-branched-chain fatty acid biosynthesis" # even
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8174","Subsystem"] <- "iso-branched-chain fatty acid biosynthesis" # odd
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8173","Subsystem"] <- "anteiso-branched-chain fatty acid biosynthesis"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY66-429","Subsystem"] <- "fatty acid biosynthesis initiation"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY66-385","Subsystem"] <- "dTMP de novo biosynthesis"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY1A0-6120","Subsystem"] <- "streptorubin B biosynthesis"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY1ZNC-1","Subsystem"] <- "assimilatory sulfate reduction IV"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8088","Subsystem"] <- "dipicolinate biosynthesis"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8073","Subsystem"] <- "lipid IV biosynthesis"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8017-2","Subsystem"] <- "L-tryptophan degradation"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8099","Subsystem"] <- "tetrahydropteridine recycling"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8131","Subsystem"] <- "5'-deoxyadenosine degradation II"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8178","Subsystem"] <- "pentose phosphate pathway (non-oxidative branch) II"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8183","Subsystem"] <- "L-valine degradation III"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8184","Subsystem"] <- "L-isoleucine degradation III"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8185","Subsystem"] <- "L-leucine degradation V"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8188","Subsystem"] <- "L-alanine degradation VI"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8200","Subsystem"] <- "backdoor pathway of androgen biosynthesis"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY-8112","Subsystem"] <- "coenzyme F420 biosynthesis"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY1YI0-3","Subsystem"] <- "succinate to cytochrome c oxidase via plastocyanin"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY1YI0-6","Subsystem"] <- "succinate to cytochrome c oxidase"
df_rxn2subsys[df_rxn2subsys$Subsystem == "PWY3O-1743","Subsystem"] <- "D-mannose degradation II"
df_rxn2subsys[df_rxn2subsys$Subsystem == "ARGDEGRAD-PWY2","Subsystem"] <- "arginine degradation"

# df_rxn2subsys[df_rxn2subsys$Subsystem == "EC-STARCH-N27","Subsystem"]
# df_rxn2subsys[df_rxn2subsys$Subsystem == "PEPTDOGLYCANSYN-PWY2","Subsystem"] <- "Peptidoglycan biosynthesis II"
# df_rxn2subsys[df_rxn2subsys$Subsystem == "TRPSYN-PWY2","Subsystem"] <- "L-tryptophan biosynthesis II"
# df_rxn2subsys[df_rxn2subsys$Subsystem == "FOLATE-SYN-GS1","Subsystem"] <- "folate biosynthesis"

## Create a short version of subsystem labels
df_rxn2subsys$SubsystemSplit <- str_simplifyDescriptions(df_rxn2subsys$Subsystem)
# oder by SubsysID, and reorder columns
df_rxn2subsys <- df_rxn2subsys[order(df_rxn2subsys$SubsysID),c("SubsysID", "ReactionID", "Subsystem", "SubsystemSplit")]
# save
saveRDS(df_rxn2subsys,"./df_rxn2subsys20230510MetaMouse.rds")
rm(lst_rxn2subsys, df_ModelSubsysNameMapping, obj_metabolicModelsMouse202305)


