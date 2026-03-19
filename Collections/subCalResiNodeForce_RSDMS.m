function nodeForce=subCalResiNodeForce_RSDMS(info,plaStress)
% Calculate residual nodal forces from plastic stresses
%
% Inputs:
%   info       - Struct containing element information including:
%                  elemBMat: Cell array of B-matrices at Gauss points
%                  GaussVol: Vector of Gauss point volumes
%   plaStress  - Plastic stress at Gauss points (cell array)
%
% Output:
%   nodeForce  - Residual nodal forces for each element (cell array)

numTimePoint=size(plaStress,2);
elemBMatCell=repmat(info.elemBMat,1,numTimePoint);
GaussVolCell=repmat(info.GaussVol(:),1,numTimePoint);

% Convert plastic stress at Gauss points to nodal forces for each element and Gauss point
nodeForce=cellfun(@(b,p)transpose(b)*p,elemBMatCell,plaStress,'UniformOutput',0);

% Multiply each nodal force contribution by corresponding Gauss volume (integration weight)
nodeForce=cellfun(@(vol,x)vol.*x,num2cell(GaussVolCell),nodeForce,'UniformOutput',0);

% Sum nodal forces over all Gauss points to get total residual nodal force per element
nodeForce=cellfun(@(x)sum(cat(3,nodeForce{:,x}),3),num2cell(1:numTimePoint),'UniformOutput',0);

end
