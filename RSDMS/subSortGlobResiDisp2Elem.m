function resiDispRate = subSortGlobResiDisp2Elem(info, gResiDispRate)
% Sort element-level residual displacement rates from global residual displacement rates
%
% Inputs:
%   info          - Structure containing element data, including DOF mapping
%   gResiDispRate - Cell array of global residual displacement rate vectors (per load case)
%
% Output:
%   resiDispRate  - Cell array of element-level residual displacement rate vectors

% Get the global DOF indices corresponding to the element's DOFs
posElemDoFInGlobDoF = info.posElemDoFInGlobDoF;

% Sort element DOFs from the global residual displacement vectors
resiDispRate = cellfun(@(x) x(posElemDoFInGlobDoF), gResiDispRate, 'UniformOutput', false);

end
