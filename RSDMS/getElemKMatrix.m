function info=getElemKMatrix(info,elaMatrix)
%%% Compute the element stiffness matrix EleKMat
% Inputs:
%   info.elemBMat   - Cell array of strain-displacement matrices B at Gauss points
%   info.GaussVol   - Weights and volumes for each Gauss integration point
%   elaMatrix       - Constitutive matrix C (elasticity matrix)
% Outputs:
%   info.elemKMat   - Computed element stiffness matrix Ke

%Calculate stiffness contribution at each Gauss point: K_i = B_i' * D * B_i
elemKMat=cellfun(@(x)transpose(x)*elaMatrix*x,info.elemBMat,'UniformOutput',0);

%Multiply each contribution by the corresponding Gauss weight/volume
elemKMat=cellfun(@(k,vol)k*vol,elemKMat,num2cell(info.GaussVol(:)),'UniformOutput',0);

%Sum all Gauss point contributions to form the total element stiffness matrix
elemKMat=sum(cat(3,elemKMat{:}),3);

% Store the result back into info structure
info.elemKMat=elemKMat;

end
