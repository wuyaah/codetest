function [ABMatColHardNorm,listUHardCell,listXHardCell]=getABMatColHardFunc(indexHardGauss,ngRelHard)
%original eq: ABmatHard*[U_pi; X_pi; beta2], beta2=(sigmaU-sigmaY)/sigmaU
%present  eq: ABmatHard*beta2*[U_pi'; X_pi'; 1], [U_pi'; X_pi']=[U_pi; X_pi]/beta2
%[U_pi; X_pi]=[U_pi'; X_pi']*beta2

global NumSet ABMatSparse

selNGRelaHard=ngRelHard(indexHardGauss);
selNGRelaHardU=kron(selNGRelaHard,ones(NumSet.NumU1Comp,1));

indexHardGauss=num2cell(indexHardGauss);
listUHardCell=cellfun(@(x)NumSet.NumU1Comp*(x-1)+1:NumSet.NumU1Comp*x,indexHardGauss,'UniformOutput',0);
listXHardCell=cellfun(@(x)NumSet.NumU1+(NumSet.NumX1Comp*(x-1)+1:NumSet.NumX1Comp*x),indexHardGauss,'UniformOutput',0);
listUHard=cell2mat(listUHardCell);
listXHard=cell2mat(listXHardCell);

%Extract the columns of the AB matrix to construct the hardening model
BMatColHard=ABMatSparse(:,listXHard);
AMatColHardNorm=ABMatSparse(:,listUHard)*diag(sparse(selNGRelaHardU)); % 为了简化背应力Cone约束的右端项“1”
ABMatColHardNorm=[AMatColHardNorm,BMatColHard];

end

