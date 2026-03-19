function VarPos=getSharedKronElemDof(sharedElemID,sharedElemNGK,numStrComp)

sharedElemIDCell=num2cell(sharedElemID);
totalElemStrComp=num2cell(numStrComp.*sharedElemNGK);
VarPosIni=cellfun(@(x,n)kron(x,ones(1,n)),sharedElemIDCell,totalElemStrComp,'UniformOutput',0);
VarPosIni=cellfun(@(x,n)(x-1).*n,VarPosIni,totalElemStrComp,'UniformOutput',0);
VarPosInc=cellfun(@(x)1:x,totalElemStrComp,'UniformOutput',0);
VarPos=cellfun(@(x,y)x+y,VarPosIni,VarPosInc,'UniformOutput',0);
VarPos=cell2mat(transpose(VarPos));

end