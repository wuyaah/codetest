function resultSet=SubGurobiCone_EPP(eleStressFull,NGYield,parameters,numVert,vertIndex,degreeSet,barQCPConvTol)

clear model
global CMatSparse OriFuncHandles
OriFuncHandles=load('Ori.FuncHandles.mat');

%% WeightSet=[1,1,0];  %%% 载荷权重 Const Load-y + varying Load-x
% AlterDegree=deg2rad([0;10;20;30;40;50;60;70;80;90]);
% degreeSet=deg2rad(degree);
weightSet=[cos(degreeSet),sin(degreeSet),zeros(numel(degreeSet),1)];
% WeightSet=[sin(AlterDegree).*cos(degreeSet) sin(AlterDegree).*sin(degreeSet) cos(AlterDegree)];         

% Critical Positions
numSet=getNumSetFunc(parameters,numVert);
alphaPos=numSet.NumU1+numSet.NumX1+1;
uInitPos=numSet.NumU1Comp*(0:numSet.NumTotalGauss-1)+1; % Initial position of every U matrix
if numVert>1
    unCellStart=numSet.NumU1+numSet.NumX1+2+numSet.NumU1Comp*(0:(numVert-1)*numSet.NumTotalGauss-1);
    uInitPos=[uInitPos,unCellStart]; 
end

%% Gurobi Paras
params.QCPDual=1;
params.OutputFlag=0; %% 控制计算输出日志
params.BarQCPConvTol=barQCPConvTol; %% 对偶变量求解精度
% params.FeasibilityTol=1e-4;
% params.OptimalityTol=1e-4;
model.ub=[ inf.*ones(numSet.totalVar-1,1);1];  % Set the upper bound to minus infinite
model.lb=[-inf.*ones(numSet.totalVar-1,1);1];  % Set the lower bound to minus infinite   
model.lb(alphaPos)=0;                          % Load factor should be a positive number
resultSet=cell(numel(degreeSet),1);

%% Start looping through all load combination
for iter=1:numel(degreeSet)

    weight=weightSet(iter,:);
    vertWeight=diag(weight)*vertIndex;
    vertElaStress=eleStressFull*vertWeight;  
    w1Combined=CMatSparse*vertElaStress(:,1); % Vert1ElaStress=OriVertElaStress(:,1);

    %% get the LinConst
    linConst=getEqConstrBesideAlphaCol(w1Combined,numSet);
    columnToAlpha=full(linConst(:,alphaPos));
    linConst=getColumnToAlpha(linConst,columnToAlpha,vertElaStress,NGYield,numSet);
    LinConstEnrich=[linConst,sparse(size(linConst,1),1)]; % Extend the matrix with the position of Beta
    
    %% get the ConeConst
    betaPos=size(LinConstEnrich,2);
    QuadIndicesCell=cellfun(@(x)x:x+numSet.NumU1Comp-1,num2cell(uInitPos),'UniformOutput',0); 
    ConeIndicesCell=cellfun(@(x)[betaPos,x],transpose(QuadIndicesCell),'UniformOutput',0);
    model=ConeModeling(model,alphaPos,LinConstEnrich,ConeIndicesCell);
    
    %% Solve the Rroblem
    result=gurobi(model,params); % gurobi_write(model,'ModelCheck.lp'); %% mps、rew、lp
    result=GurobiPostProcess(result,LinConstEnrich,NGYield,vertElaStress,numSet);
    result.vertWeight=vertWeight;%(1:size(oriElaStress,2),:);
    result.vertElaStress=vertElaStress;
    result.ColumnToAlpha=columnToAlpha;
    result.degree=degreeSet(iter);
    result.numSet=numSet;
    kktErr=GurobiKKTvalid(model,alphaPos,LinConstEnrich,QuadIndicesCell,result);
    result.kktErr=kktErr;
    
    resultSet{iter}=result;
    
    %%
    if numel(degreeSet)>1
        plotShakedownDomain(resultSet,NGYield(1));  % Optional: plot shakedown limit curve
    end
    % disp({'LoadFactor:' result.alpha;
    %       'MaxYield:'       result.ProcessedRes.MaxYield;
    %       'EquiliViolation:'    result.ProcessedRes.EquiViolateAB;
    %       'UAdditionConvertTolerance:' result.ProcessedRes.UAdditionConvertTolerance;}); 
end
%%
% rmpath(gurobiPath);

end
