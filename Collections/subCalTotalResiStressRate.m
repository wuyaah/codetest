function resiStressRate = subCalTotalResiStressRate(info, resiDispRate, elaMatrix)
% Calculate total residual stress rate from residual displacement rates and elastic matrix
%
% Inputs:
%   info          - Structure containing element information (including elemBMat, IntE)
%   resiDispRate  - Cell array of residual displacement rate vectors per load case
%   elaMatrix     - Elasticity matrix (material stiffness)
%
% Output:
%   resiStressRate - Cell array of residual stress rate tensors per Gauss point and load case

% Scale elastic matrix by material intensity factor IntE
elaMatrix = info.IntE .* elaMatrix;

% Compute residual strain rate at Gauss points by multiplying strain-displacement matrices with displacement rates
GaussResiStrainRate = cellfun(@(x) cellfun(@(b) b * x, info.elemBMat, 'UniformOutput', 0), resiDispRate, 'UniformOutput', 0);
GaussResiStrainRate=cat(2,GaussResiStrainRate{:});

% Calculate residual stress rate by applying elastic matrix to residual strain rate
resiStressRate = cellfun(@(x) elaMatrix * x, GaussResiStrainRate, 'UniformOutput', 0);

end
