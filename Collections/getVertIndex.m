function vertIndex=getVertIndex(vertNorm,numVert)
% NumVert: NV=1:limit; NV=2:Prop shakedown; % NV=4: Shakedown; NV=8: 3D Shakedown
if vertNorm==0
vertIndex=[1 0 1 0 1 0 1 0;  % load1
           1 0 0 1 1 0 0 1;  % load2
           1 0 0 0 0 1 1 1]; % load3
elseif vertNorm==1
    vertIndex=[1 1 0 1 1 0 0 1;  % load1
               1 0 1 0 1 0 1 0;  % load2
               1 0 0 0 0 1 1 1]; % load3
    % vertIndex=[1 0 1 0 1 0 1 0;  % load1
    %            1 1 0 1 1 0 0 1;  % load2
    %            1 0 0 0 0 1 1 1]; % load3
end

numVertSet=[1,2,3,4,8];
if ismember(numVert,numVertSet)
    if numVert==3 && vertNorm==0
        vertIndex=vertIndex(:,[1,3,4]);
    else
        vertIndex=vertIndex(:,1:numVert);
    end
else
    error('Unsupported numVert value. Use 2, 3 or 4.');
end


end 

