function [elemGaussStrainSet, elemNodeDispSet, nodeConstrained] = LoadModelData(DataName)
%% Load element Gauss strains, node displacements, and node constraints from file
% Inputs:
%   DataName - Filename containing saved model data
%
% Outputs:
%   elemGaussStrainSet - Cell array of element strains at Gauss points
%   elemNodeDispSet    - Cell array of element node displacement vectors
%   nodeConstrained    - Matrix or cell array indicating constrained DOFs for nodes

% Load necessary variables from file
loadData = load(DataName, 'NodeConstrainedDof', 'EleGaussStrainSet', 'EleNodeDispSet');

% Extract node constraints and element node displacement data
nodeConstrained = loadData.NodeConstrainedDof;
elemNodeDispSetRaw = loadData.EleNodeDispSet;

% Convert displacement data from 3D array to cell array per element
seqElem = num2cell(transpose(1:size(elemNodeDispSetRaw, 1))); % element indices 1 to N
elemNodeDispSet = cellfun(@(x) squeeze(elemNodeDispSetRaw(x, :, :)), seqElem, 'UniformOutput', 0);
elemNodeDispSet = cellfun(@(x) cell2mat(x), elemNodeDispSet, 'UniformOutput', 0);

% Attempt to load and reshape Gauss strain data if present
try
    elemGaussStrainSetRaw = loadData.EleGaussStrainSet;
    if size(elemGaussStrainSetRaw, 3) > 1 % if multiple Gauss points per element
        elemGaussStrainSet = cellfun(@(x) squeeze(elemGaussStrainSetRaw(x, :, :)), seqElem, 'UniformOutput', 0);
        elemGaussStrainSet = cellfun(@(x) cell2mat(x), elemGaussStrainSet, 'UniformOutput', 0);
    else
        elemGaussStrainSet = elemGaussStrainSetRaw;
    end
catch
    % If EleGaussStrainSet does not exist, return empty
    elemGaussStrainSet = [];
end

% If nodeConstrained has more than 3 elements (rows), convert to cell array per node and transpose
if numel(nodeConstrained) > 3
    [Row, Column] = size(nodeConstrained);
    nodeConstrained = mat2cell(nodeConstrained, ones(Row, 1), Column);
    nodeConstrained = cellfun(@transpose, nodeConstrained, 'UniformOutput', 0);
    NodeConstrainedDof = nodeConstrained;
    % Save the reformatted node constraints back into the file for future use
    save(DataName, 'NodeConstrainedDof', '-append');
end

end
