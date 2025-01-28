import gurobipy
import cobra
from cobra.io import load_model, read_sbml_model
from cobra.flux_analysis import flux_variability_analysis
#import cobra.test
from corpse import simpleFastcore
import random
import pickle
import numpy
import pandas
import re
import os.path

## In a first step, build a consistent model through a first pass 
## with FVA and subsequently several passes of fastcc.
## Subsequently, fastcore is being run on the model using the precomputed core reactions.

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def setdiff(lst1, lst2):
    return list(set(lst1) - set(lst2))

###########################################
##   Create consistent model
###########################################
##(step can be skipped if consistent model has already been created)

##Load model
metamodel = read_sbml_model('../models/metamodel_20230510.xml')
metamodel.solver="gurobi"

##First pass (FVA)
core_reaks = random.sample([x.id for x in metamodel.reactions], 5) # just five random reactions
fast_meta = simpleFastcore(model = metamodel, core_set = core_reaks, zero_cutoff=1e-4)
fast_meta.FVA_consistency()

f = open('../models/metamodel_20230510_consistent.obj', 'wb')
pickle.dump(fast_meta.get_model(), f ) 
f.close()

##Next passes with fastcc
for i in range(10):
    metamodel = fast_meta.get_model()
    core_reaks = random.sample([x.id for x in metamodel.reactions], 5) # just five random reactions
    fast_meta = simpleFastcore(model = metamodel, core_set = core_reaks, zero_cutoff=1e-4)
    fast_meta.fastcc()
    f = open('../models/metamodel_20230510_consistent.obj', 'wb')
    pickle.dump(fast_meta.get_model(), f )
    f.close()

#######################################################
##   Run fastcore
#######################################################
##Load Fastcc Output
f = open('../models/metamodel_20230510_consistent.obj', 'rb')
metamodel = pickle.load( f )
metamodel.solver = "gurobi"
origModel = metamodel.copy() ##Store original version of the model

##Adjust lower and upper bounds of reactions to prevent numerical issues
for reak in metamodel.reactions:
    if reak.lower_bound < 0 and reak.lower_bound > -1e-2:
        ##print(reak.id)
        reak.lower_bound = -1e-2

for reak in metamodel.reactions:
    if reak.lower_bound < -10 and reak.lower_bound > -1000:
        ##print(reak.id)
        reak.lower_bound = -10

##Load core reactions
core = pandas.read_csv("../processedData/coreReakMat.tsv",sep="\t")

##Run fastcore on a sample and subsequently flux variability analysis on the output
for i in range(1,53):
    if not(os.path.isfile('../results/mouseModels/'+core.columns.values[i]+'.reaks.csv')):
        ##Write a dummy file to indicate which model is currently been generated
        pandas.DataFrame([x.id for x in metamodel.reactions]).to_csv('../results/mouseModels/'+core.columns.values[i]+'.reaks.csv',index=False,header=False)
        ##Get core reactions for sample
        sub = core[core.iloc[:,i]>0]
        core_reaks = sub['reaction_id'].values
        ##Adjust reaction names
        core_reaks = [re.sub("\\(","_LPAREN_", x) for x in core_reaks]
        core_reaks = [re.sub("\\)","_RPAREN_", x) for x in core_reaks]
        ##Build fastcore object
        fast_meta = simpleFastcore(model = metamodel, core_set = core_reaks, zero_cutoff=1e-4)
        ##Run fastcore
        fast_meta.fastcore()
        csModel = fast_meta.get_model()
        ##Store list of reactions
        pandas.DataFrame([x.id for x in csModel.reactions]).to_csv('../results/mouseModels/'+core.columns.values[i]+'.reaks.csv',index=False,header=False)
        ##Get host and microbiome reactions
        reakList = [x.id for x in csModel.reactions]
        colonReaks = [match for match in reakList if "_Colon" in match]
        brainReaks = [match for match in reakList if "_Brain" in match]
        liverReaks = [match for match in reakList if "_Liver" in match]
        hostReaks = colonReaks + brainReaks + liverReaks
        microbiomeReaks = [match for match in reakList if "_Microbiome" in match]
        ##restore original lower bounds
        exReaks = [match for match in [x.id for x in csModel.reactions] if match.startswith("EX_")]
        for reak in exReaks:
            csModel.reactions.get_by_id(reak).lower_bound=origModel.reactions.get_by_id(reak).lower_bound
        ###################
        #    Host side
        ###################
        ##Run FVA on the WT
        fvaHostWT = flux_variability_analysis(csModel,hostReaks,fraction_of_optimum=0,processes=10)
        ##Run FVA on GF
        csModelGF = csModel.copy()
        for reak in microbiomeReaks:
            csModelGF.reactions.get_by_id(reak).lower_bound=0
            csModelGF.reactions.get_by_id(reak).upper_bound=0
        fvaHostGF = flux_variability_analysis(csModelGF,hostReaks,fraction_of_optimum=0,processes=10)
        ##Result 
        fvaHost = pandas.concat([fvaHostWT,fvaHostGF],axis=1)
        fvaHost.columns = ['minFluxFull', 'maxFluxFull', 'minFluxGF', 'maxFluxGF']
        fvaHost.to_csv('../results/mouseModels/'+core.columns.values[i]+'.fva-mic.csv')
        #########################
        #    Microbiome side
        #########################
        ##Run FVA on the WT
        fvaMicWT = flux_variability_analysis(csModel,microbiomeReaks,fraction_of_optimum=0,processes=10)
        ##Run FVA on GF
        csModelHF = csModel.copy()
        for reak in hostReaks:
            csModelHF.reactions.get_by_id(reak).lower_bound=0
            csModelHF.reactions.get_by_id(reak).upper_bound=0
        fvaMicHF = flux_variability_analysis(csModelHF,microbiomeReaks,fraction_of_optimum=0,processes=10)
        ##Result
        fvaMic = pandas.concat([fvaMicWT,fvaMicHF],axis=1)
        fvaMic.columns = ['minFluxFull', 'maxFluxFull', 'minFluxHF', 'maxFluxHF']
        fvaMic.to_csv('../results/mouseModels/'+core.columns.values[i]+'.fva-host.csv')

