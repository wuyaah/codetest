function [nodeForceVec, constrainedDoFVec] = OriModelNodeData(dataName, eleNodeLabelSet, alternativeLoadCases)
% Load nodal force vectors and constrained degrees of freedom (DOFs) from model data.
% NOTE: This implementation assumes all constrained DOFs are fixed to zero (Dirichlet BC).
%
% INPUTS:
%   dataName            : Main .mat file name (contains primary model data)
%   eleNodeLabelSet     : Matrix of element-to-node connectivity
%   alternativeLoadCases: Cell array of alternative load case filenames (without .mat)
%
% OUTPUTS:
%   nodeForceVec        : Global nodal force vectors (columns = load cases)
%   constrainedDoFVec   : Vector of constrained DOF indices

% -------------------------------------------------------------------------
% Load nodal information from data file (variables starting with "Node")
loadData = load(dataName, '-regexp', '^(Node)...');
nodeLabelList = loadData.NodeLabelList;  % Node index list
if iscell(nodeLabelList)
    nodeLabelList = cell2mat(nodeLabelList);
end

%% Initialize nodal force vectors (may include multiple load cases)
[nodeForceVec,Dimension]=getNodeForceVec(loadData, eleNodeLabelSet, nodeLabelList);

% Load additional nodal forces from alternative load case files (if any)
for loadCase = alternativeLoadCases
    fileName = [cell2mat(loadCase) '.mat'];
    if exist(fileName, 'file')
        alternativeLoadData = load(fileName, '-regexp', '^(Node)...');
        alternativeNodeForceVec = getNodeForceVec(alternativeLoadData, eleNodeLabelSet, nodeLabelList);
        nodeForceVec = [nodeForceVec, alternativeNodeForceVec];  % Append
    end
end

%% Initialize constrained DOFs
% Generate constraint matrix (size: #nodes ¡Á DOFs per node)
ConstrainedDoFMat = OriNodeConstrained(loadData.NodeConstrainedDof, nodeLabelList);
constrainedDoFVec = reshape(transpose(ConstrainedDoFMat(:, 1:Dimension)), [], 1);
constrainedDoFVec = constrainedDoFVec .* transpose(1:numel(constrainedDoFVec));
constrainedDoFVec = constrainedDoFVec(constrainedDoFVec > 0);

end
