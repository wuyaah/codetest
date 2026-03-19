function [elaStress, elemYield, elemRelHard] = subExtractStress(dataPack)
%% Extract element Gauss point stresses and related material parameters
% Inputs:
%   dataPack  - Cell array of element data structs containing stress and material info
%
% Outputs:
%   elaStress - Element stresses evaluated at Gauss points
%   ngYield   - Yield stress repeated at each Gauss point for all elements
%   ngRelHard - Relative hardening parameter at Gauss points, calculated as (HardLimit - Yield)/Yield

% Get stresses at Gauss points for all elements
elaStress = getElemGaussStress(dataPack);

% Repeat hardening limits and yield values at each Gauss point for each element
% elemHardLimit = cellfun(@(x) repmat(x.HardLimit, x.NGK, 1), dataPack, 'UniformOutput', 0);
% elemYield = cellfun(@(x) repmat(x.Yield, x.NGK, 1), dataPack, 'UniformOutput', 0);
elemHardLimit = cellfun(@(x) x.HardLimit, dataPack);
elemYield = cellfun(@(x) x.Yield, dataPack);

% Calculate relative hardening parameter at Gauss points
elemRelHard = (elemHardLimit - elemYield) ./ elemYield;

% Optional:
% To use maximum hardening parameter over all Gauss points:
% ngRelHard(:) = max(ngRelHard);

end
