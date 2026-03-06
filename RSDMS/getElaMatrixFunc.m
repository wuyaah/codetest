function elaMatrix=getElaMatrixFunc(dimension,v)
% getElaMatrixFunc - Construct elastic (constitutive) matrix for FEM.
%
% INPUT:
%   dimension : 2 or 3 (for 2D plane stress or full 3D)
%   v         : Poisson's ratio
%
% OUTPUT:
%   elaMatrix : Elasticity matrix (C), sparse format

if dimension==2
    % Plane stress condition (e.g., CPS)
    elaMatrix=1/(1-v^2)*[ 1 v 0; v 1 0; 0 0 (1-v)/2 ];
elseif dimension==3
    % 3D isotropic elasticity
    G=1/(2*(1+v));
    elaMatrix=inv([ 1 -v -v  0  0  0; 
                   -v  1 -v  0  0  0;
                   -v -v  1  0  0  0;
                    0  0  0 1/G 0  0;
                    0  0  0  0 1/G 0;
                    0  0  0  0  0 1/G]);
else
    error('Dimension must be 2 or 3.');
end

elaMatrix=sparse(elaMatrix);

end