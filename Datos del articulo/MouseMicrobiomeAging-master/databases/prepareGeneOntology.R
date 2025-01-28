###################################################################
###
###   Prepare a GeneOntology database for use in clusterProfiler
###
###################################################################


setwd("databases/")


###############
## Mouse gene Ontology (newest version as of May 2023)
# http://purl.obolibrary.org/obo/go/go-basic.obo  <- descriptions of GOids in a weird format
# http://current.geneontology.org/annotations/mgi.gaf.gz  <- mouse specific translation table of gene to GOid
## Load mouse specific gene to GO table
df_mouseGOannotaion <- read.csv(file = "~/Downloads/mgi.gaf.gz",
                                header = F, sep = "\t", comment.char = "!", strip.white = T,
                                blank.lines.skip = T, stringsAsFactors = F)
# drop non usefull columns
df_mouseGOannotaion <- df_mouseGOannotaion[,c(3:5,9:10)]
colnames(df_mouseGOannotaion) <- c("GeneSymbol", "FunctionalRole", "GOid", "OntologyType", "GeneName")
# convert to factors, where it's meaningfull
df_mouseGOannotaion$FunctionalRole <- as.factor(df_mouseGOannotaion$FunctionalRole)
df_mouseGOannotaion$OntologyType <- factor(x = df_mouseGOannotaion$OntologyType,
                                           levels = c("P", "F", "C"),
                                           labels = c("biological process", "molecular function", "cellular component"))
## load table for annotating GO-ids
# Note this file was created with bash-script: "databases/getGOannotation.sh"
#  In case of interesting terms with missing descriptions, consult original obo file or online GO database
df_GOannotation <- read.csv(file = "GOBioProcessAnnotationSimple.csv.gz",
                            header = F, sep = "\t", strip.white = T, stringsAsFactors = F)
colnames(df_GOannotation) <- c("GOid", "Description", "OntologyType")
# convert that one to a factor as well
df_GOannotation$OntologyType <- factor(x = trimws(sub(pattern = "_", replacement = " ", fixed = T, x = df_GOannotation$OntologyType)),
                                       levels = c("biological process", "molecular function", "cellular component"))
# merge the two tables, because why would you keep that info seperated in the first place?
df_tmp <- merge(x = df_mouseGOannotaion, y = df_GOannotation[,1:2], by = "GOid", all.x = T)
df_mouseGOannotaion <- df_tmp; rm(df_tmp, df_GOannotation)
# since we are not using a bunch of columns, we have duplicate entries - get rid of those
df_mouseGOannotaion <- df_mouseGOannotaion[!duplicated(df_mouseGOannotaion),]
# and keep for later, thanks!
saveRDS(df_mouseGOannotaion, paste0("df_mouseGOannotaion_", Sys.Date(),".rds"))

