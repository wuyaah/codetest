function resiStressRateEla = subResiDerivativeRate(totalComp,plaComp)
% Compute the elastic part of the residual stress rate by subtracting plastic part
%
% Inputs:
%   totalComp - Cell array containing total residual stress rate components
%   plaComp   - Cell array containing plastic stress rate components
%
% Output:
%   resiStressRateEla - Cell array of elastic residual stress rate components

% Element-wise subtraction: elastic = total - plastic
resiStressRateEla=cellfun(@(x,y)x-y,totalComp,plaComp,'UniformOutput',0);
resiStressRateEla=cellfun(@(x)transpose(x),resiStressRateEla,'UniformOutput',0);

% Convert inner cell arrays to matrices
numTimePoint=size(resiStressRateEla,2);
numTimePointCell=transpose(num2cell(1:numTimePoint));
resiStressRateEla=cellfun(@(x)cell2mat(resiStressRateEla(:,x)),numTimePointCell,'UniformOutput', 0);

end
