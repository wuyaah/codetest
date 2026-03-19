function linConst=getColumnToAlpha(linConst,oriColumnToAlpha,oriVertElaStress,eleYield,numSet)

global ABMatSparse OriFuncHandles
alphaPos=numSet.NumU1+numSet.NumX1+1;
ConvertVariable=OriFuncHandles.ConvertVariable;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The numbering order of indices:
% 2 Vertices:   %4 Vertices:    %8 Vertices:
%       / 1     %  4------ 1    %   6-------1
%     /         %  |       |    %  /|      /|
% 2 /           %  |       |    % 7-|----8/ |
%               %  2 ------3    % | |     | |
                                % |/2 --- | /4
                                % 3------ 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numVert=numSet.numVert;
numU1Comp=numSet.NumU1Comp;
numStrComp=numSet.numStrComp;
numTotalGauss=numSet.NumTotalGauss;
% OriVertElaStress=OriElaStress*Weight;  % Orignal stress at each vertex

% Choose the conversion function depending on the element type
EleYieldCell=num2cell(repmat(eleYield,1,numVert-1));
if numStrComp==6; ConvertFunc=ConvertVariable.Sig2UX3D;      % 3D Model
elseif numStrComp==4; ConvertFunc=ConvertVariable.Sig2UXCPE; % Plane Strain
elseif numStrComp==3; ConvertFunc=ConvertVariable.Sig2UXCPS; % Plane Stress
end
eleConvertMat=cellfun(@(x)ConvertFunc(x),EleYieldCell,'UniformOutput',0);
eleConvertMat=cellfun(@(x)x(1:numU1Comp,:),eleConvertMat,'UniformOutput',0);

diffVertStress=oriVertElaStress(:,2:end)-repmat(oriVertElaStress(:,1),1,numVert-1);
diffVertStressCell=mat2cell(diffVertStress,numStrComp*ones(numTotalGauss,1),ones(1,numVert-1));
diffVertU=cellfun(@(c,y)c*y,eleConvertMat,diffVertStressCell,'UniformOutput',0);
diffVertU=cell2mat(diffVertU);

UX1RHS=sparse(size(ABMatSparse,1),1); % Right hand side corresponds to UX1 are 0s
additionalColumnToAlpha=[UX1RHS;diffVertU(:)];

columnToAlpha=sparse(oriColumnToAlpha+additionalColumnToAlpha);
linConst(:,alphaPos)=columnToAlpha;

end
