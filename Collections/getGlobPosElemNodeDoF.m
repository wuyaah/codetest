function info = getGlobPosElemNodeDoF(info, nodeLabelList)
% Map local element DOFs (degrees of freedom) to global system DOFs.
%
% This function determines how the degrees of freedom (DOFs) for a single
% finite element correspond to the global system DOFs, based on the node
% labeling and element-node connectivity. It stores the mapping and 
% indexing information in the `info` structure for later use in assembly
% processes like stiffness matrix and force vector assembly.
%
% INPUT:
%   info           : Struct containing:
%                    - info.nodeLabel: Node labels in the element
%                    - info.dimension: Spatial dimension (2 or 3)
%   nodeLabelList  : Global list of all node labels in the model
%
% OUTPUT:
%   info : Updated struct with fields:
%          - posElemNodeInGlobList   : Indices of element nodes in the global node list
%          - posElemDoFInGlobDoF     : Indices of element DOFs in the global DOF vector
%          - posElemMatInGlobMatRow  : Row indices for assembling element matrices into global matrix
%          - posElemMatInGlobMatColumn : Column indices for assembling element matrices into global matrix

% Extract problem dimension (e.g., 2D or 3D)
dimension = info.dimension;

% Get node labels for the current element
elemNodeLabelSet = transpose(info.nodeLabel(:));  % Force column to row vector

% Find position of each element node in the global node label list
[~, posElemNodeInGlobList] = ismember(elemNodeLabelSet, nodeLabelList);

% Compute global DOF indices for this element
% Each node contributes 'dimension' DOFs
% e.g., for 2D: node i => DOFs 2*(i-1)+[1 2]
posElemDoFInGlobDoF = dimension * (posElemNodeInGlobList - 1) + transpose(1:dimension);
posElemDoFInGlobDoF = posElemDoFInGlobDoF(:);  % Flatten to column vector

% Generate row and column index grids for element matrix assembly
[posElemMatInGlobMatColumn, posElemMatInGlobMatRow] = meshgrid(posElemDoFInGlobDoF, posElemDoFInGlobDoF);

% Store computed indices back into the info struct
info.posElemMatInGlobMatRow     = posElemMatInGlobMatRow;
info.posElemMatInGlobMatColumn  = posElemMatInGlobMatColumn;
info.posElemNodeInGlobList      = posElemNodeInGlobList;
info.posElemDoFInGlobDoF        = posElemDoFInGlobDoF;

end
