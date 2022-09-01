function error=objective_fun(kcat)
% global ecModel_batch;
kcats_inf=load('kcats_inf.mat');
kcats_inf=kcats_inf.kcats_inf;
% global C_sources;
% global data;
% global N_sources;
C_sources=load('C_sources.mat');
N_sources=load('N_sources.mat');
data=load('data.mat');
ecModel_batch=load('ecYeastGEM_batch.mat');
ecModel_batch=ecModel_batch.ecModel_batch;
data=data.data;
C_sources=C_sources.C_sources;
N_sources=N_sources.N_sources;
C_num=length(C_sources);
N_num=length(N_sources);
data_pool=[];
i=0;
for j=1:C_num
    for k=1:N_num
        i=i+1;
        data_pool(i,1)=data(k,j);
        data_pool(i,2)=j;
        data_pool(i,3)=k; 
    end
end
posGluc = strcmp(ecModel_batch.rxnNames,'D-glucose exchange (reversible)');
posNH4  = strcmp(ecModel_batch.rxnNames,'ammonium exchange (reversible)');
ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posGluc),0,'u');
ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posNH4),0,'u');
for j=1:length(kcat)
%     a=cell2mat(kcats_inf{j}{1});
    a1=kcats_inf{j}{1};
    a2=kcats_inf{j}{2};
    ecModel_batch.S(a1,a2)=kcat(1,j);
end
for k=1:i
    C_s=data_pool(k,2);
    N_s=data_pool(k,3);
    if strcmp(C_sources{C_s},'(R)-lactate')
        vC = [+Inf,+Inf];
    else
        vC = +Inf;
    end
    ecModel_lim=changeMediaCN(ecModel_batch,C_sources{C_s},N_sources{N_s},vC,+Inf);
    solution=solveLP(ecModel_lim);
    result(1,k)=-solution.f;
end

error=nanmean(abs(data_pool(data_pool(:,1)>0,1)'-result(1,(data_pool(:,1)>0)))./data_pool((data_pool(:,1)>0),1)')*100;
end

function model = changeMediaCN(model,C_source,N_source,C_value,N_value)
model = openUptake(model,C_source,C_value);
model = openUptake(model,N_source,N_value);
if strcmp(C_source,'D-fructose')
    model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1134')) = 0;
    model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1134')) = 0;
elseif strcmp(C_source,'D-mannose')
    model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1139')) = 0;
    model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1139')) = 0;
end
end
function model = openUptake(model,source,value)
if strcmp(source,'(R)-lactate')
    model = openUptake(model,'(S)-lactate',value(2));
    value = value(1);
end
posB = strcmp(model.rxnNames,[source ' exchange (reversible)']);
posF = strcmp(model.rxnNames,[source ' exchange']);
if sum(posB) == 0
    model = changeRxnBounds(model,model.rxns(posF),-value,'l');
else
    model = changeRxnBounds(model,model.rxns(posB),value,'u');
end
model = changeRxnBounds(model,model.rxns(posF),0,'u');
end