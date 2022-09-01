posGluc = strcmp(ecModel_batch.rxnNames,'D-glucose exchange (reversible)');
posNH4  = strcmp(ecModel_batch.rxnNames,'ammonium exchange (reversible)');
posO2  = strcmp(ecModel_batch.rxnNames,'oxygen exchange (reversible)');
ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posGluc),inf,'u');
ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posNH4),inf,'u');
ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posO2),inf,'u');
result=[];
for i =1:51
    for j =1:51
        ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posO2),0.2*(i-1),'u');
        ecModel_batch = changeRxnBounds(ecModel_batch,ecModel_batch.rxns(posGluc),0.2*(j-1),'u');
        sol = optimizeCbModel(ecModel_batch);
        result(j,i)=sol.f;
    end
end

s=surf(result)