This folder contains analysis scripts for performing the modeling steps described in the paper starting from the mapping of data to the generic metamodel and analysis of the resulting context-specific models for the indidivual mice. Scripts are relatively standalone but some require installation of additional libraries as indicated in the headers of each script. Scripts should be executed in the indicated order.

Content of the individual scripts:

1_createCoreReactions.R
-----------------------
Script for deriving core reactions from transcriptomic and metagenomic data. Preprocesses transcriptomic as well as metagenomic data for StanDep first and maps them to identifiers in the model. Subsequently, 2_runStandep.m needs to be run to execute StanDep and results are collected to create a core reaction matrix afterwards.

2_runStandep.m
--------------
Runs StanDep to derive core reactions based on expressiona and metagenomic data. The script requires installation of StanDep and the CobraToolbox in Matlab (https://github.com/LewisLabUCSD/StanDep).

3_runFastcore.py
----------------
Executes fastcore implemented in the troppo toolbox to reconstruct context-specific models for each mouse based on the core reactions derived from StanDep. Subsequently, models are characterized using flux variability analysis both with and without microbial reactions present.

4_fastcorePostProcessing.R
--------------------------
Further characterization of the derived context-specific models including a correction of FVA results for exchange reactions and determination of maximal production rates for individual metabolites.

5_summarizeResults.R
--------------------
Summarizes results from step 3 and 4.
