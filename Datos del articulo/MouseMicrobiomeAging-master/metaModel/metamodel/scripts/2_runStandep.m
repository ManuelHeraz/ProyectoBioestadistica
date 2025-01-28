%% Script to convert expression and metagenomic data into core reactions
%% Processes host and microbiome data separately
%% Outputs core reaction matrices which have to be combined into a common matrix
%% As input conditions, tissue_age is used

initCobraToolbox(0)
addpath('<Add path to StanDep codes>')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Host data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expression = dlmread('../processedData/standep/expression-host.tsv');
fid = fopen('../processedData/standep/genes-host.tsv');
genes = textscan(fid,'%s');
fclose(fid)
fid = fopen('../processedData/standep/conditions-host.tsv');
tissues = textscan(fid,'%s');
fclose(fid)
load('../models/recon22.mat');

rnaData.gene=genes{1};
rnaData.valuebyTissue=expression;
rnaData.Tissue=tissues{1};

%Create StanDep data structure
modelData = getModelData(rnaData,model);

%Obtain enzyme types
spec = getSpecialistEnzymes(model);  
prom = getPromEnzymes(model);

%calculate enzyme expression
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);

%Further parameters
edgeX = [-2 -1 0 1 2 2.5 3 4]; % bins  
distMethod = '@chi2dist'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering

%Number of clusters (optimized through screening of different cluster numbers)
k=39; %first average distance below 0.05 up to 60 clusters

%calculate clusters of enzyme expression
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,k,distMethod,linkageMethod);

%calculate core reaction sets
coreRxnMat = int8(models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],false,0,[1 1])); 

for i=1:size(expression,2)
    %%Derive unexpressed reactions
    disp(i);
    %%Replace missing genes by arbitrary value since we don't know
    %%their actual expression value
    miss=setdiff(model.genes,genes{1});
    test=struct;
    test.gene=[genes{1}; miss;];
    test.value=[expression(:,i); ones(length(miss),1)*1000;];

    [t1 t2]=mapExpressionToReactions(model,test);
    
    zeroList=find(t1==0);
    coreRxnMat(zeroList,i)=-1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Code for scanning for optimal cluster number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% arr={};

% for i=1:20
%     %Number of clusters
%     k=i*3;

%     %calculate clusters of enzyme expression
%     clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,k,distMethod,linkageMethod);

%     %calculate core reaction sets
%     coreRxnMat = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],false,0,[1 1]); 
    
%     arr{i}=reshape(coreRxnMat,[numel(coreRxnMat),1]);
% end

% dist = zeros(20,20);

% for i=1:19
%     for k=(i+1):20
%         dist(i,k)=pdist([arr{i} arr{k}]','jaccard');
%     end
% end


%store coreRxnMat
%save '../data/coreRxns.mat' coreRxnMat;
csvwrite('../processedData/coreRxns-host.csv',coreRxnMat);
writetable(cell2table(model.rxns),'../processedData/reactions-host.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Microbiome data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expression = dlmread('../processedData/standep/expression-microbiome.tsv');
fid = fopen('../processedData/standep/genes-microbiome.tsv');
genes = textscan(fid,'%s');
fclose(fid)
fid = fopen('../processedData/standep/conditions-microbiome.tsv');
tissues = textscan(fid,'%s');
fclose(fid)
load('../models/metamodel_20230510.mat'); %%Use full metamodel
model.subSystems = cell(size(model.rxns,1),1);

for i=1:size(model.subSystems,1)
    model.subSystems{i,1} = {''};
end

rnaData.gene=genes{1};
%rnaData.valuebyTissue=expression*1e8;
rnaData.valuebyTissue=expression;
rnaData.Tissue=tissues{1};

%Create StanDep data structure
modelData = getModelData(rnaData,model);

%Obtain enzyme types
spec = getSpecialistEnzymes(model);  
prom = getPromEnzymes(model);

%calculate enzyme expression
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);

%Further parameters
%edgeX = [-2 -1 0 1 2 2.5 3 4]; % bins  
distMethod = '@chi2dist'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering

%Redefine bins appropriately
minLog=log10(min(enzymeData.value(enzymeData.value>0)))
maxLog=log10(max(enzymeData.value(enzymeData.value>0)))
edgeX=minLog:((maxLog-minLog)/7):maxLog;

%Number of clusters (result of scanning)
k=15;

%calculate clusters of enzyme expression
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,k,distMethod,linkageMethod);

%calculate core reaction sets
coreRxnMat = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],true,0,[1 1]); 

%store coreRxnMat
csvwrite('../processedData/coreRxns-microbiome.csv',coreRxnMat);
writetable(cell2table(model.rxns),'../processedData/reactions-metamodel.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Code for scanning for optimal cluster number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%arr={};

%%for i=1:20
%%    %Number of clusters
%%    k=i*3;
%%
%%    %calculate clusters of enzyme expression
%%    clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,k,distMethod,linkageMethod);
%%
%%    %calculate core reaction sets
%%    coreRxnMat = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],false,0,[1 1]); 
%%    
%%    arr{i}=reshape(coreRxnMat,[numel(coreRxnMat),1]);
%%end
%%
%%dist = zeros(20,20);
%%
%%for i=1:19
%%    for k=(i+1):20
%%        dist(i,k)=pdist([arr{i} arr{k}]','jaccard');
%%    end
%%end
