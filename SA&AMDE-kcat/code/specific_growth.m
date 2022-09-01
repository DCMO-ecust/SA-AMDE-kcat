load('growth_data.mat');
data_res =data;
% load('ecmodel.mat');
ecModel=ecModel_batch;
load('yeastGEM8.3.4.mat');
load('C_sources.mat')
load('N_sources.mat')
posGluc = strcmp(ecModel.rxnNames,'D-glucose exchange (reversible)');
posNH4  = strcmp(ecModel.rxnNames,'ammonium exchange (reversible)');
ecModel = setParam(ecModel,'ub',ecModel.rxns(posGluc),0);
ecModel = setParam(ecModel,'ub',ecModel.rxns(posNH4),0);

mu_max     = zeros(size(data));

for i = 1:length(N_sources)
    for j = 1:length(C_sources)
        if strcmp(C_sources{j},'(R)-lactate')
            vC = [+Inf,+Inf];
        else
            vC = +Inf;
        end
        
        ecModel_lim = changeMediaCN(ecModel,C_sources{j},N_sources{i},vC,+Inf);
        sol_lim     = optimizeCbModel(ecModel_lim);
        mu_max(i,j) = sol_lim.f;
        if isempty(sol_lim.x)
             mu_max(i,j)=1;
        end
        disp(['Growing on ' N_sources{i} ' / ' C_sources{j}])
    end
end



