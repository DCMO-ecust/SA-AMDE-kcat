c_source = 'r_1714';
model=readCbModel('yeastGEM8.3.4.xml');
model = ravenCobraWrapper(model);
model = minimal_Y6(model);
ecModel_batch=readCbModel('result149onlymodel_5.mat');
posGluc = strcmp(ecModel_batch.rxnNames,'D-glucose exchange (reversible)');
posNH4  = strcmp(ecModel_batch.rxnNames,'ammonium exchange (reversible)');
ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posGluc),inf,'u');
ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posNH4),inf,'u');
[FVA_batch,~,stats_b]     = comparativeFVA(model,ecModel_batch,c_source,false,0);
FVA_batch     = filterDistributions(FVA_batch,1E-10);
distributions = {FVA_batch{1},FVA_batch{2}};
legends       = {'model-batch', 'ecModel-batch'};
titleStr      = 'Flux variability cumulative distribution';
[~, ~]        = plotCumDist(distributions,legends,titleStr);

function newDist = filterDistributions(distribution,treshold)
newDist = [];
for i=1:length(distribution)
    dist = distribution{i};
    dist(dist<treshold) =treshold;
    dist(find(dist==1E-10,1)) = treshold/10;
    newDist = [newDist,{dist}];
end
end

function model = minimal_Y6(model)
exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = 0;
model.ub(exchangeRxns) = 1000;
desiredExchanges = {'r_1654';'r_1992';'r_2005';'r_2060';'r_1861';'r_1832';
    'r_2100'; 'r_4593';'r_4595'; 'r_4596';'r_4597';'r_2049';'r_4594'; 'r_4600'; 'r_2020' };  
blockedExchanges = {'r_1663';'r_4062';'r_4064'}; 
glucoseExchange = {'r_1714'};
uptakeRxnIndexes     = findRxnIDs(model,desiredExchanges);
glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
BlockedRxnIndex      = findRxnIDs(model,blockedExchanges);
if length(find(uptakeRxnIndexes~= 0)) ~= 15
    warning('Not all exchange reactions were found.')
end
model.lb(uptakeRxnIndexes(uptakeRxnIndexes~=0))     = -1000;
model.lb(glucoseExchangeIndex) = -1;
model.lb(BlockedRxnIndex) = 0;
model.ub(BlockedRxnIndex) = 0;
end
