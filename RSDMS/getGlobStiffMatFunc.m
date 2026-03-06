function globKMat=getGlobStiffMatFunc(dataPack,numNode)
%% Assemble the global stiffness matrix from element stiffness matrices
% Inputs:
%   dataPack       - Cell array of element data structs, each containing:
%                    elemKMat    - Element stiffness matrix (before scaling)
%                    IntE        - Element material property scalar (e.g. interpolated Young's modulus)
%                    posElemMatInGlobMatRow - Row indices for element DOFs in global matrix
%                    posElemMatInGlobMatColumn - Column indices for element DOFs in global matrix
%   numNode        - Total number of nodes in the model
%
% Outputs:
%   globKMat       - Sparse global stiffness matrix assembled from all elements

% Total degrees of freedom in global stiffness matrix
dimGlobK = numNode * dataPack{1}.dimension;

% Extract the global row indices for all element stiffness entries
KRow = cellfun(@(x) x.posElemMatInGlobMatRow(:), dataPack, 'UniformOutput', 0);

% Extract the global column indices for all element stiffness entries
KColumn = cellfun(@(x) x.posElemMatInGlobMatColumn(:), dataPack, 'UniformOutput', 0);

% Scale element stiffness matrices by element material property IntE, then reshape as column vectors
elemKMat = cellfun(@(x) reshape(x.IntE .* x.elemKMat, [], 1), dataPack, 'UniformOutput', 0);

% Assemble the global stiffness matrix as a sparse matrix using all elements' contributions
globKMat = sparse(cell2mat(KRow), cell2mat(KColumn), cell2mat(elemKMat), dimGlobK, dimGlobK);

end
