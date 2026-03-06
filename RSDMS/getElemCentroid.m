function [dataPack, elemCentroidSet] = getElemCentroid(dataPack)
% Compute and store element centroid coordinates in dataPack.
%
% INPUT:
%   dataPack         : Cell array of element structures, each with field 'nodeCoord'
%                      (nodeCoord is an N x D matrix of nodal coordinates per element)
%
% OUTPUT:
%   dataPack         : Updated with new field `Centroid` for each element
%   elemCentroidSet  : Matrix (numElems x D) of centroid coordinates

% -------------------------------------------------------------------------
% Compute centroid of each element by averaging its nodal coordinates
elemCentroidSet = cellfun(@(x) mean(x.nodeCoord), dataPack, 'UniformOutput', 0);

% Assign centroid to each element's data struct
for elem = 1:numel(dataPack)
    dataPack{elem}.Centroid = elemCentroidSet{elem};
end

% Convert to numeric matrix for external use
elemCentroidSet = cell2mat(elemCentroidSet);

end
