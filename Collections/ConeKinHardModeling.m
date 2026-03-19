function [model,ConeHardIndicesCell]=ConeKinHardModeling(model,alphaPos,linConstHardEnrichIncompress,ConeIndicesCell)
global NumSet

% objective
model.modelsense='max';
model.obj=full(sparse(1,alphaPos,1,1,size(linConstHardEnrichIncompress,2)));

% equality
model.A=linConstHardEnrichIncompress;
model.rhs=zeros(size(linConstHardEnrichIncompress,1),1);
model.sense='=';

% inequality constraints
for m=1:numel(ConeIndicesCell)
    model.cones(m).index=ConeIndicesCell{m};
end

% inequality constraints of back stress
numIneq=numel(ConeIndicesCell);
betaHardPos=size(linConstHardEnrichIncompress,2);
ConeHardIndicesCell=cell(NumSet.NumHardGauss,1);
for n=1:NumSet.NumHardGauss
   ConeHardIndices=NumSet.totalVarPerfect+(NumSet.NumU1Comp*(n-1)+1:NumSet.NumU1Comp*n);
   ConeHardIndicesCell{n}=ConeHardIndices;
   model.cones(numIneq+n).index=[betaHardPos,ConeHardIndices];
end

end

