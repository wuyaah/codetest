function resultSet=ShakedownCodeNorm(iskinemHardening,elaStress,elemYield,elemRelHard,parameters,degree,numVert,vertNorm,barQCPConvTol)
% Main driver function for shakedown analysis using Gurobi Solver
%
% Inputs:
%   iskinemHardening - Flag indicating kinematic hardening model (1) or not (0)
%   elaStress        - Elastic stress data at Gauss points
%   GPYield          - Yield stress at Gauss points
%   GPRelHard        - Relative hardening parameters at Gauss points (only used if kinematic hardening)
%   parameters       - Extracted model parameters
%   numVert          - Number of vertices for quadrature/approximation
%   degree           - Polynomial degree for shape functions or approximation
%   vertNorm         - Norm values associated with vertices
%   barQCPConvTol    - (Optional) Convergence tolerance for BarQCP parameter in optimization
% Output:
%   resultSet        - Cell array containing optimization results, including qcpi for convergence

% Add Gurobi path to MATLAB environment
gurobiDir=('C:\Gurobi952\win64\matlab');
addpath(gurobiDir);
%%
degreeSet=deg2rad(degree);

%% SelVert:
vertIndex=getVertIndex(vertNorm,numVert);
%%
elaFactor=1;
eleStressFull=elaFactor.*[elaStress,zeros(size(elaStress,1),3-size(elaStress,2))];
%%
NGK=parameters{3}(1);
NGRelHard=reshape(transpose(repmat(elemRelHard,1,NGK)),[],1);
NGYield=reshape(transpose(repmat(elemYield,1,NGK)),[],1);

if nargin<9
    % If convergence tolerance not provided:
    % Perform shakedown topology optimization with default tolerance
    barQCPConvTol=1e-5;
end

while true
    if iskinemHardening==0  % Call standard shakedown Gurobi cone optimization
        resultSet=SubGurobiCone_EPP(eleStressFull,NGYield,parameters,numVert,vertIndex,degreeSet,barQCPConvTol);
    elseif iskinemHardening==1  % Call kinematic hardening Gurobi cone optimization
        resultSet=SubGurobiCone_KH(eleStressFull,NGYield,NGRelHard,parameters,numVert,vertIndex,degreeSet,barQCPConvTol);
    else
        error('Input error: iskinemHardening must be 0 or 1.');
    end
    barQCPConvTol=barQCPConvTol/1e2;     % Tighten convergence tolerance for next iteration
    if isfield(resultSet{end},'qcpi')    % Iteratively refine tolerance until qcpi parameter is found (converged)
        break;
    end
end

% Remove Gurobi path from MATLAB environment
rmpath(gurobiDir);

end


