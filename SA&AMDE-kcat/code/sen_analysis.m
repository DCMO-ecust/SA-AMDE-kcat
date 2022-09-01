clear
clc
load('ecYeastGEM_batch.mat');
p=parpool();
p.IdleTimeout=inf;
load('data.mat');
load('C_sources.mat');
C_num=length(C_sources);
load('N_sources.mat');
N_num=length(N_sources);
data_pool=[];
i=1;
for j=1:C_num
    for k=1:N_num
        data_pool(i,1)=data(k,j);
        data_pool(i,2)=j;
        data_pool(i,3)=k;
        i=i+1;
    end
end
posGluc = strcmp(ecModel_batch.rxnNames,'D-glucose exchange (reversible)');
posNH4  = strcmp(ecModel_batch.rxnNames,'ammonium exchange (reversible)');
ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posGluc),0,'u');
ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posNH4),0,'u');

protposition=[];
for i=1:length(ecModel_batch.mets)
    if ismember('prot',cell2mat(ecModel_batch.mets(i)))
        protposition=[protposition;i];
    end
end

protposition(length(protposition),:)=[];
for i=1:length(ecModel_batch.rxnNames)
    if ismember('draw_prot_',cell2mat(ecModel_batch.rxnNames(i)))
        drawrxn=i;
        break
    end
end

kcats_inf={};
kcat=[];
k=1;
for i=1:length(protposition)
    posi=find(ecModel_batch.S(protposition(i),:)~=0);
    for j=1:length(posi)
        if posi(j)>drawrxn-1
            continue
        end
        kcats_inf{k}={protposition(i),posi(j),ecModel_batch.S(protposition(i),posi(j))};
        k=k+1;
        kcat=[kcat,ecModel_batch.S(protposition(i),posi(j))];
    end
end
kcat=full(kcat);
% train_pool=data_pool;
% train_data=train_pool((train_pool(:,1)>0),1);
for i=1:length(kcat)
    child_kcat=kcat;
    child_kcat(i)=kcat(i)*2;
    for j=1:length(kcat)
        ecModel_batch.S(cell2mat(kcats_inf{j}(1)),cell2mat(kcats_inf{j}(2)))=child_kcat(j);
    end
    parfor k=1:230
        C_s=train_pool(k,2);
        N_s=train_pool(k,3);
        if strcmp(C_sources{C_s},'(R)-lactate')
            vC = [+Inf,+Inf];
        else
            vC = +Inf;
        end
        ecModel_lim=changeMediaCN(ecModel_batch,C_sources{C_s},N_sources{N_s},vC,+Inf);
        solution=solveLP(ecModel_lim);
        result_child(i,k)=-solution.f;
    end
    
    error_child(i)=nanmean(abs((train_data(:,1)'-result_child(i,(train_pool(:,1)>0)))./train_data(:,1)'))*100;
end


ind=cell(1,15);
result=abs((error_child-objective_fun(kcat))/objective_fun(kcat));
for i=1:15
num(i)=length(find(result>0.1^i));
ind(i)={find(result>0.1^i)};
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
    