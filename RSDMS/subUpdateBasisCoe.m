function [a0Coe, akCoe, bkCoe] = subUpdateBasisCoe(a0Coe, akCoe, GPoint, weight, resiStressRate)
% UpdateBasisCoe: Update the Fourier basis coefficients and reconstruct residual stress
%
% Inputs:
%   a0Coe          - Cell array of current constant (a₀) coefficients
%   akCoe          - Cell array of current sine coefficients (aₖ)
%   GPoint         - Gauss-Lobatto points (column vector)
%   weight         - Integration weights (column vector)
%   paramBasis     - {cosComp, sinComp} basis values at Gauss points
%   resiStressRate - Cell array of current residual stress rate (per element)
%
% Outputs:
%   a0Coe     - Updated constant term coefficients
%   akCoe     - Updated sine term coefficients
%   bkCoe     - Updated cosine term coefficients
%   resiStress - Updated residual stress (based on new coefficients)

% Store previous coefficients
a0Coe_old=a0Coe;
akCoe_old=akCoe;

% Convert vectors to cell format for broadcasting
GPoint=num2cell(GPoint);
weight=num2cell(weight);
numTerms=num2cell(transpose(1:numel(akCoe{1})));

% Update sine (ak) and cosine (bk) coefficients using current stress rate
akCoe=cellfun(@(x)SubBasisSinUpdate(x,GPoint,weight,numTerms),resiStressRate,'UniformOutput',0);
bkCoe=cellfun(@(x)SubBasisCosUpdate(x,GPoint,weight,numTerms),resiStressRate,'UniformOutput',0);

% Compute correction to a₀ based on change in ak
a0CoeComp=cellfun(@(r,ak1,ak)SubConstCoeUpdate(r,weight,ak1,ak),resiStressRate,akCoe,akCoe_old,'UniformOutput',0);

% Update a₀
a0Coe = cellfun(@(x,comp)sparse(x+comp),a0Coe_old,a0CoeComp,'UniformOutput',0);

% Reconstruct residual stress from updated Fourier coefficients

end
