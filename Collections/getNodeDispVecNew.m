function disp=getNodeDispVecNew(globKMat,nodeForceVec,knownDofs,knownValues)


if nargin < 4 || isempty(knownValues)
    knownValues=zeros(numel(knownDofs),1);
end

dimGlobK=size(globKMat,1);
numLoadCases=size(nodeForceVec,2);

allDofs=(1:dimGlobK)';
unknownDofs=setdiff(allDofs,knownDofs);

globKMat_rr=globKMat(unknownDofs,unknownDofs);
globKMat_rk=globKMat(unknownDofs,knownDofs);
force_adj=globKMat_rk*knownValues;

reducedForceMatrix=nodeForceVec(unknownDofs,:)-force_adj;

disp_r=globKMat_rr\reducedForceMatrix;
disp_matrix=zeros(dimGlobK,numLoadCases);
disp_matrix(knownDofs,:)=repmat(knownValues,1,numLoadCases);
disp_matrix(unknownDofs,:)=disp_r;

disp=mat2cell(disp_matrix,dimGlobK,ones(1,numLoadCases));

end

