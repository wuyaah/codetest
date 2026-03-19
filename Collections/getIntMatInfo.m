function info = getIntMatInfo(info, matInfo)
% Compute interpolated material properties based on element relative density.
%
% This function implements SIMP-based interpolation of material properties,
% commonly used in topology optimization, and stores both the interpolated
% values and their derivatives (gradients) in the info struct.
%
% INPUT:
%   info     : Element data struct, must include field `topoRelDens` (relative density).
%   matInfo  : Material parameter cell array, parsed by getMatInfo function.
%
% OUTPUT:
%   info     : Updated struct with interpolated material properties and gradients.

% Extract base and minimal material parameters + penalization powers
[E, Emin, Yield, YieldMin, Limit, LimitMin, Penal, YPenal, UPenal] = subExtMatInfo(matInfo);

% Element relative density
rho = info.topoRelDens;

% Interpolated material properties
intE     = Emin   + (E     - Emin)   * rho^Penal;
intYield = YieldMin + (Yield - YieldMin) * rho^YPenal;
intLimit = LimitMin + (Limit - LimitMin) * rho^UPenal;

% Gradients of interpolated properties w.r.t. relative density
gradE = Penal  * (E     - Emin)   * rho^(Penal - 1);
gradY = YPenal * (Yield - YieldMin) * rho^(YPenal - 1);
gradU = UPenal * (Limit - LimitMin) * rho^(UPenal - 1);

% Store in struct
info.IntE       = intE;
info.Yield      = intYield;
info.HardLimit  = intLimit;
info.GradMat    = [gradE, gradY, gradU];

end
