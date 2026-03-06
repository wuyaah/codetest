function gResiNodeForce = subSortNodeForce2Model_RSDMS(dataPack, resiNodeForceRate, numNode)
% Assemble element residual nodal forces into global residual force vectors
%
% Inputs:
%   dataPack          - Cell array containing element data, including DOF positions:
%                         posElemDoFInGlobDoF (positions of element DOFs in global DOFs)
%   resiNodeForceRate - Cell array of residual nodal force rate vectors for each element
%   numNode           - Total number of nodes in the model
%
% Output:
%   gResiNodeForce    - Cell array of sparse global residual force vectors for each load case

dimGlobF = numNode * dataPack{1}.dimension; % Total global DOFs

% Extract global DOF indices for all elements and concatenate into one vector
FRow = cellfun(@(x) x.posElemDoFInGlobDoF, dataPack, 'UniformOutput', false);
FRow = cell2mat(FRow);

% Convert each element residual force cell array to a numeric matrix and concatenate
resiNodeForceRate = cellfun(@(x) cell2mat(x), resiNodeForceRate, 'UniformOutput', false);
resiNodeForceRate = cell2mat(resiNodeForceRate);

% Split concatenated residual forces column-wise into separate cells (one per load case)
[dim1, dim2] = size(resiNodeForceRate);
resiNodeForceRate = mat2cell(resiNodeForceRate, dim1, ones(1, dim2));
% resiNodeForceRate = transpose(resiNodeForceRate);

% Assemble sparse global residual force vectors for each load case
gResiNodeForce = cellfun(@(x) sparse(FRow, 1, x, dimGlobF, 1), resiNodeForceRate, 'UniformOutput', false);

end
