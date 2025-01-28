###############################################
###
###     Use this script to explore single genes 
###   for their correlation to microbiome functions
###
###############################################

#################
## Load functions and resources
source("../customFunctions.R")


#################
## Load strong host-microbiome correlation results (spearman-rho >= 0.55 & FDR <= 0.1)
df_pcorRxn2ColonRNATop <- readRDS("df_pcorRxn2ColonRNATop.rds")
df_pcorRxn2LiverRNATop <- readRDS("df_pcorRxn2LiverRNATop.rds")
df_pcorRxn2BrainRNATop <- readRDS("df_pcorRxn2BrainRNATop.rds")


#######################
#### Checking out Succinate related genes
lst_succinateRelatedGenesFromMM
# [1] "Slc13a3" "Sirt1"   "Sirt2"   "Sirt3"   "Sirt4"   "Sirt5"   "Pdhx"    "Dld"     "Dlst"    "Ogdh"    "Acaa1"   "Acox1"   "Cpt1c"   "Cpt1a"   "Cpt1b"   "Alas1"  
# [17] "Alas2"   "Gcat"    "Bbox1"   "Ehhadh"  "Hsd17b4" "Crot"    "Mut"     "Oxct2"   "Oxct1"   "Plod1"   "Plod2"   "Plod3"   "Phyh"    "P4ha1"   "P4hb"    "P4ha2"  
# [33] "Acot8"   "Acot4"   "Tmlhe"   "Slc25a1"
View(df_DifAbundColonAllAges[df_DifAbundColonAllAges$GeneSymbol %in% lst_succinateRelatedGenesFromMM &
                             df_DifAbundColonAllAges$padj <= 0.05,])
View(df_DifAbundBrainAllAges[df_DifAbundBrainAllAges$GeneSymbol %in% lst_succinateRelatedGenesFromMM &
                             df_DifAbundBrainAllAges$padj <= 0.05,])
View(df_DifAbundLiverAllAges[df_DifAbundLiverAllAges$GeneSymbol %in% lst_succinateRelatedGenesFromMM &
                             df_DifAbundLiverAllAges$padj <= 0.05,])
lapply(lst_succinateRelatedGenesFromMM,void_checkGene)

#check for aging and mb overlap
# like Alas1 Alas2 colon liver brain, Pdhx brain, Crot liver, 

######################

#################
# example usage
void_checkGene("Tuba1a")
# [1] "Tuba1a in Colon:"
# RxnName Spearman_rho       p.adj
# 1  2'-Deoxyadenosine 5'-diphosphate:oxidized-thioredoxin 2'-oxidoreductase    0.5731966 0.013170161
# 2                                                      L-Rhamnose Exchange    0.5679894 0.015539226
# 3                              L-Rhamnulose-1-phosphate lactaldehyde-lyase    0.5679894 0.015539226
# 4                                               L-Rhamnose ketol-isomerase    0.5679894 0.015539226
# 5                                    ATP:L-rhamnulose 1-phosphotransferase    0.5679894 0.015539226
# 6                                  L-rhamnose transport via proton symport    0.5679894 0.015539226
# 7  2'-Deoxyguanosine 5'-diphosphate:oxidized-thioredoxin 2'-oxidoreductase    0.5659384 0.016298019
# 8                                       ATP:cytidine 5'-phosphotransferase    0.5646548 0.016789894
# 9      2'-Deoxycytidine diphosphate:oxidized-thioredoxin 2'-oxidoreductase    0.5605868 0.018815806
# 10                               glycinamide ribonucleotide transformylase    0.5577201 0.020000439
# 11                                    Glycerone phosphate phosphohydrolase    0.5519096 0.023281173
# 12                                      ornithine transport via ABC system   -0.5539226 0.022044491
# 13                         Prephenate:NAD+ oxidoreductase(decarboxylating)   -0.5631819 0.017445574
# 14                                               L-Arginine iminohydrolase   -0.5639815 0.017032516
# 15                                                        Urea-e0 Exchange   -0.5737091 0.013016156
# 16                                                    TRANS-RXN4LZ-7041.ce   -0.5737091 0.013016156
# 17                    Carbamoyl-phosphate:L-ornithine carbamoyltransferase   -0.5742024 0.012841597
# 18                                               Chorismate pyruvatemutase   -0.5841870 0.009680960
# 19                                         ATP:GTP 3'-diphosphotransferase   -0.6392972 0.001522053
# 20            guanosine 3'-diphosphate 5'-triphosphate 5'-phosphohydrolase   -0.6392972 0.001522053
# [1] "Tuba1a in Brain:"
#                           RxnName Spearman_rho      p.adj
# 1        Xanthosine ribohydrolase   -0.5779476 0.07050349
# 2         Guanosine ribohydrolase   -0.5831764 0.06021088
# 3 N-D-ribosylpurine ribohydrolase   -0.6046203 0.03574772


