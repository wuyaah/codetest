function info = SortElemRelDens(info, ElemRelDens)
% Assign element relative density values to the element information structure.
%
% This function stores the relative density (e.g., from topology optimization)
% for a single element into the field `topoRelDens` of the `info` struct.
%
% INPUT:
%   info        : Struct containing element-specific data.
%   ElemRelDens : Scalar or vector representing the relative density
%                 associated with this element (e.g., between 0 and 1).
%
% OUTPUT:
%   info : Updated struct with the relative density stored in `topoRelDens`.

% Store the relative density
info.topoRelDens = ElemRelDens;

end
