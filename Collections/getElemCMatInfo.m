function info = getElemCMatInfo(info, MatInfo)
% Construct and assign constitutive matrix and related element data.
%
% This function determines the correct approach to compute the constitutive
% (C) matrix and B matrix based on the element type. It also incorporates
% material properties like yield stress, hardening limit, and Poisson's ratio
% into the element info struct.
%
% INPUT:
%   info     : Struct containing element type, node coordinates, materialID, etc.
%   MatInfo  : Cell array containing material information in the form:
%              {E_modulus, YieldStress, HardeningLimit, PoissonRatio}
%
% OUTPUT:
%   info     : Updated struct with:
%              - elemCMat: Constitutive matrix
%              - elemBMat: Strain-displacement matrix
%              - GaussVol: Gauss point volume weights
%              - Yield   : Adjusted yield stress
%              - HardLimit: Adjusted hardening limit
%              - Poisson : Poisson’s ratio
%              - dimension, NGK (Gauss info), etc.

% Extract element type and coordinates
elemType = info.elemType; 
elemNodeCoordSet = info.nodeCoord;

% Extract and replicate material properties
Yield     = MatInfo{2}(1) * ones(2,1);
HardLimit = MatInfo{3}   * ones(2,1);
Poisson   = MatInfo{4};

% Define supported element types for standard processing
normalEleList = {'CPS3','CPE3','CPS4','CPS4R','CPE6','C3D4','C3D6','C3D8R','C3D10'};

% Dispatch to appropriate C-matrix builder function
if ismember(elemType, normalEleList)
    [elemCMat, elemBMat, GaussVol, elemInfo, deriShapeFunc] = getElementCMatNorm(elemNodeCoordSet, elemType);
elseif strcmp(elemType, 'C3D8')
    [elemCMat, elemBMat, GaussVol, elemInfo, deriShapeFunc] = getElementCMatC3D8(elemNodeCoordSet, elemType);
else
    disp("This EleType is not supported");
    return
end

% Assign computed fields to the info struct
info.elemBMat  = elemBMat;
info.elemCMat  = elemCMat;
info.GaussVol  = GaussVol;
info.Yield     = Yield(info.materialID)    .* (1 - info.damageCoe);
info.HardLimit = HardLimit(info.materialID).* (1 - info.damageCoe);
info.Poisson   = Poisson(info.materialID);
info.dimension = elemInfo.Dimension;
info.NGK       = elemInfo.NGK;
% info.DeriShapeFunc = deriShapeFunc;  % Uncomment if needed

% Clean up unused field
info = rmfield(info, 'damageCoe');

end
