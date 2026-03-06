function [loadingCoe,GPoint,weight,numSet]=subLoadingCoe(numSet)

numVert = numSet.numVert;
vertNorm = numSet.vertNorm;
% Generate loading paths and Gauss points/weights
[path1, path2, GPoint, weight] = subLoadingPathFunc_new(numVert,vertNorm);

% Evaluate loading coefficients at Gauss points for both load paths
loadingCoe = [path1(GPoint), path2(GPoint)];
loadingCoe = transpose(loadingCoe);

% Save path functions back to params struct
numSet.path1 = path1;
numSet.path2 = path2;

%% Display loading info
% disp(['The Loading Coe (Case1): ', sprintf('%8.4f', loadingCoe(:,1)'), newline, ...
%       'The Loading Coe (Case2): ', sprintf('%8.4f', loadingCoe(:,2)'), newline, ...
%       'The Gauss-Lobatto points: ', sprintf('%8.4f', GPoint'), newline, ...
%       'The weight factors (w/2): ', sprintf('%8.4f', weight')]);

end
