function elaLimit=getElaLimitFactor(result,elemYield)

numSet=result.numSet;
numStrComp=numSet.numStrComp;
NGK=numSet.NGK;
numElem=numSet.numElem;


EffCalFunc=getConvertFunc(numStrComp);

vertElaStress=result.vertElaStress;
elaStress=mat2cell(vertElaStress,numStrComp*NGK*ones(numElem,1),size(vertElaStress,2));
elaStress=cellfun(@(x)mat2cell(x,numStrComp*ones(NGK,1),ones(1,size(vertElaStress,2))),elaStress,'UniformOutput',0);
effElaStress=cellfun(@(x)cellfun(@(y)EffCalFunc(y),x),elaStress,'UniformOutput',0);
effElaStress=cellfun(@(x)max(x,[],1),effElaStress,'UniformOutput',0);
effElaStress=cell2mat(effElaStress);

[value,PosMax]=max(effElaStress);
elaLimit=elemYield(PosMax)./value(:);
elaLimit=min(elaLimit);

end
