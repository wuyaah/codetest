function parameters=ExtractConvert(dataPack, elemNodeLabelSet, nodeLabelList, nodeConstrainedDoF)
% Determine element type and call corresponding data extraction and conversion subroutine
%
% Inputs:
%   dataPack           - Cell array containing element data structures
%   elemNodeLabelSet   - Set of element node labels
%   nodeLabelList      - List of all node labels in the model
%   nodeConstrainedDoF - Vector or cell array of constrained degrees of freedom per node
%
% Output:
%   parameters         - Extracted and converted parameters specific to the element type

% Attempt to get element type from first element in dataPack
elemType = dataPack{1}.elemType;

% Define lists of supported 2D element types
normalElemListCPS = {'CPS4', 'CPS4R', 'CPS3'};
normalElemListCPE = {'CPE4', 'CPE4R', 'CPE3', 'CPE6'};

% Dispatch processing based on element type category
if sum(ismember(elemType, normalElemListCPS)) > 0
    parameters = SubExtractCPSConvertProcess(dataPack, nodeLabelList, elemNodeLabelSet, nodeConstrainedDoF);
elseif sum(ismember(elemType, normalElemListCPE)) > 0
    parameters = SubExtractCPEConvertProcess(dataPack, nodeLabelList, elemNodeLabelSet, nodeConstrainedDoF);
else
    parameters = SubExtract3DConvertProcess(dataPack, nodeLabelList, elemNodeLabelSet, nodeConstrainedDoF);
end

end
