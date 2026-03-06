function [coePhi, coePhiTol] = UpdateCoePhionPlaStress(coePhi_old,plaStress,weight)
% UpdateCoePhionPlaStress: Update the scalar coefficients Phi based on effective plastic stress
%
% Inputs:
%   coePhi_old   - Previous iteration's Phi coefficients (row vector)
%   effPlaStress - Cell array of effective plastic stress at Gauss points {elem}{gauss_point}
%   weight       - Integration weight (scalar)
%
% Outputs:
%   coePhi       - Updated Phi coefficients (L2 norm over all elements and Gauss points)
%   coePhiTol    - Relative change in Phi coefficients compared to previous iteration
numElem=numel(plaStress);
NGK=size(plaStress{1},1); % Number of Gauss points per element
plaStress=cat(1,plaStress{:});
effPlaStress=cellfun(@subCalEffStress,plaStress);

coePhi=vecnorm(effPlaStress,2,1); % L2 norm across all Guass points
coePhi=coePhi*weight/sqrt(numElem*NGK); % Normalize by total number of Gauss points

% Compute relative tolerance (used for convergence check)
coePhiTol=abs(coePhi-coePhi_old)./coePhi;

end
