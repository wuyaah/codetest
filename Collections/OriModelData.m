function [info, elemNodeLabelSet, totalNodeList] = OriModelData(dataName, alternativeLoadCases)
% Restructure original model data into a structured format for analysis.
%
% INPUTS:
%   dataName             : Name of the main .mat file (e.g., 'HoledPlateLoad11_DOF3')
%   alternativeLoadCases : Cell array of alternative load case filenames (without .mat)
%
% OUTPUTS:
%   info                 : Cell array of structured element information
%   elemNodeLabelSet     : Matrix of element-node connectivity (padded with zeros)
%   totalNodeList        : Unique list of all node labels in the model

%% Load element-related data (element labels, nodes, etc.)
loadData = load(dataName, '-regexp', '^(Ele)...');  % Load variables beginning with 'Ele'

% Generate sequential element indices (1, 2, ..., N)
seqElemLable = num2cell(transpose(1:numel(cell2mat(loadData.EleLabelSet))));

% Process and organize the loaded data by element order
loadData = loadDateProcessing(loadData, seqElemLable);

% Convert element-node connectivity data from cell of cells to matrix
elemNodeLabelSet = cellfun(@(x)cell2mat(x), loadData.EleNodeLabelSet, 'UniformOutput', 0);

% Extract all unique node labels used in the model
totalNodeList = cell2mat(elemNodeLabelSet);
totalNodeList = unique(totalNodeList);

% Ensure each element's node label list is a row vector
elemNodeLabelSet = cellfun(@(x)transpose(x), elemNodeLabelSet, 'UniformOutput', 0);
loadData.EleNodeLabelSet = elemNodeLabelSet;

%% Construct structured info for each element
info = cellfun(@(x)getInfoStructGen(loadData, x), seqElemLable, 'UniformOutput', 0);

% If alternative load cases exist, update stress info accordingly
info = getAlternativeLoadGuassStress(info, seqElemLable, alternativeLoadCases);

%% Add alternative nodal displacement results if available
for loadCase = alternativeLoadCases
    fileName = [cell2mat(loadCase) '.mat'];
    if exist(fileName, 'file')
        loadData = load(fileName, 'EleNodeDispSet');
        elemNodeDispSet = cellfun( ...
            @(x)squeeze(loadData.EleNodeDispSet(x,:,:)), ...
            seqElemLable, 'UniformOutput', 0);
        info = cellfun(@(x, disp)getAlternativeNodeDisp(x, disp), ...
            info, elemNodeDispSet, 'UniformOutput', 0);
    end
end

%% Post-process: pad incomplete element-node sets with zeros
elemNodeLabelSet = getModifiedElemNodeLable(elemNodeLabelSet);  % e.g., convert 3-node to 4-node with padding
elemNodeLabelSet = cell2mat(elemNodeLabelSet);  % Final element connectivity matrix

end
