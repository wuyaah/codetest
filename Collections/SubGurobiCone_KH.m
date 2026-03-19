function resultSet=SubGurobiCone_KH(elaStressFull,ngYield,ngRelHard,parameters,numVert,vertIndex,degreeSet,barQCPConvTol)
% toltalVar=[U1;X1;-alpha;U2;...;U8;beta1;U_pi;X_pi;beta2]; beta(=1) is used to construct the cone constraint
% problem statement: maximize alpha
clear model
global NumSet CMatSparse OriFuncHandles
OriFuncHandles=load('Ori.FuncHandles.mat');
% gurobiPath=('C:\Gurobi9\win64\matlab');
% addpath(gurobiPath); 
%%
% degreeSet=deg2rad([0;30;45;60;90]);
% degreeSet=deg2rad([0;10;20;30;40;50;60;70;80;90]);
% degreeSet=deg2rad([0;30;45;60;90;120;135;150;180;210;225;240;270;300;315;330;360]);
% degreeSet=deg2rad(90);
% degreeSet=deg2rad(0);

% degreeSet=deg2rad(degreeSet);
weightSet=[cos(degreeSet),sin(degreeSet),zeros(numel(degreeSet),1)];

NumSet=getNumSetFunc(parameters,numVert);
NumSet.NumHardGauss=NumSet.NumTotalGauss;
NumSet.totalVarPerfect=NumSet.NumVert*(NumSet.NumU1)+NumSet.NumX1+2; %Number of total variables for elastic-perfect plastic material
NumSet.hardVar=NumSet.NumHardGauss*(NumSet.NumU1Comp+NumSet.NumX1Comp)+1; %Number of variables to construct the hardening (U_Pi, X_Pi, Beta_Pi)
NumSet.totalVar=NumSet.totalVarPerfect+NumSet.hardVar;
NumSet.NumIneqs=NumSet.NumIneqs+NumSet.NumTotalGauss; %Number of total Ineqs for kinecmatic hardening material contains the back stress constraint.

%% two-surface model for bounded kinematic hardening
% Binder.sigmaY=280;  % Yield strength of the binder phase
% Binder.sigmaU=341;  % Ultimate strength of the binder phase
% Binder.relHard=(Binder.sigmaU-Binder.sigmaY)./Binder.sigmaY;
% elemYield=kron(elemYield,ones(NumSet.NGK,1)); % Note: eleYield is already the yield strength of every GAUSSIAN point.
% elemRelHard=kron(elemRelHard,ones(NumSet.NGK,1));

%% Critical Positions
alphaPos=NumSet.NumU1+NumSet.NumX1+1;
UInitPos=NumSet.NumU1Comp*(0:NumSet.NumTotalGauss-1)+1; % Initial position of every U matrix
if numVert>1
    UNCellStart=NumSet.NumU1+NumSet.NumX1+2+NumSet.NumU1Comp*(0:(numVert-1)*NumSet.NumTotalGauss-1);
    UInitPos=[UInitPos,UNCellStart]; 
end

%% Start looping through all load combination
params.QCPDual=1;
params.OutputFlag=0; %控制计算输出日志
params.BarQCPConvTol=barQCPConvTol; %对偶变量求解精度
model.ub=[ inf.*ones(NumSet.totalVarPerfect-1,1);1; inf.*ones(NumSet.hardVar-1,1);1]; % Set the upper bound to minus infinite
model.lb=[-inf.*ones(NumSet.totalVarPerfect-1,1);1;-inf.*ones(NumSet.hardVar-1,1);1]; % Set the lower bound to minus infinite
model.lb(alphaPos)=0;  % load factor should be a positive number
resultSet=cell(numel(degreeSet),1);
vertIndex=vertIndex(:,1:numVert);

%% Introducing Hardening
indexHardGauss=1:NumSet.NumHardGauss; %Indices of Gaussian that should be enriched with hardening model
[ABMatColHard,listUHardCell,listXHardCell]=getABMatColHardFunc(indexHardGauss,ngRelHard); 

%% Start looping through all load combination
for iter=1:numel(degreeSet)

    weight=weightSet(iter,:);
    vertWeight=diag(weight)*vertIndex;
    vertElaStress=elaStressFull*vertWeight;  
    w1Combined=CMatSparse*vertElaStress(:,1); % Vert1ElaStress=OriVertElaStress(:,1);

    %% get the LinConst
    linConst=getEqConstrBesideAlphaCol(w1Combined);
    columnToAlpha=full(linConst(:,alphaPos));
    linConst=getColumnToAlpha(linConst,columnToAlpha,vertElaStress,ngYield); %Update the Column to Alpha
    linConstEnrich=[linConst,sparse(size(linConst,1),1)]; %Extend the matrix with the position of Beta (=1)
    zeroRowsHard=sparse(size(linConstEnrich,1)-size(ABMatColHard,1),size(ABMatColHard,2)); %Zero rows appear during the introduction of hardening
    linConstHard=[linConstEnrich,[ABMatColHard;zeroRowsHard]];
    linConstHardEnrich=[linConstHard,sparse(size(linConstHard,1),1)]; % Extend the matrix with the position of Beta2
    
    %% Imcompressible condition of the Back Stress
    ngYieldHard=ngYield(indexHardGauss);  % Yield Strength of elements with hardening behavior
    incompressCell=cellfun(@(x)getIncompressCoe(x,NumSet.numStrComp),num2cell(ngYieldHard),'UniformOutput',0);    % values used to construct the incompressible condition
    incompressCellNorm=cellfun(@(x,i)[x.*i(1:2),i(3:end)],num2cell(ngRelHard),incompressCell,'UniformOutput',0);     % 为了简化背应力Cone约束的右端项“1”
    incompressValue=transpose(cell2mat(incompressCellNorm));
    rowIndexIncompress=(kron(1:numel(indexHardGauss),ones(1,size(incompressValue,1)))).';              % list of row number (used to construct the incompressible condition)
    columIndexIncompress=cellfun(@(u,x)[u(1:2),x],listUHardCell,listXHardCell,'UniformOutput',0);      % [U_pi1,U_pi2,X_pi]
    columIndexIncompress=NumSet.totalVarPerfect+cell2mat(columIndexIncompress);    
    incompressMat=sparse(rowIndexIncompress,columIndexIncompress(:),incompressValue(:),numel(indexHardGauss),NumSet.totalVar);
    linConstHardEnrich=[linConstHardEnrich;incompressMat]; %#ok<AGROW> %Introducing the incompressible condition
    
    clear rowIndexIncompress columIndexIncompress incompressValue incompressCell
    clear linConst zeroRowsHard linConstHard
    
    %%
    % Inequality constraints
    betaPos=size(linConstEnrich,2);
    QuadIndicesCell=cellfun(@(x)x:x+NumSet.NumU1Comp-1,num2cell(UInitPos),'UniformOutput',0);
    ConeIndicesCell=cellfun(@(x)[betaPos,x],transpose(QuadIndicesCell),'UniformOutput',0);
    
%     model=ConeKinHardModeling(model,alphaPos,linConstHardEnrichIncompress,ConeIndicesCell);
    [model,ConeHardIndicesCell]=ConeKinHardModeling(model,alphaPos,linConstHardEnrich,ConeIndicesCell);

    %% Solve the Rroblem
    result=gurobi(model,params);
    result=GurobiPostProcessHard(result,linConstHardEnrich,ngYield,ngRelHard,vertElaStress,indexHardGauss);  
    result.ColumnToAlpha=linConstHardEnrich(:,alphaPos);
    result.VertWeight=vertWeight;
    result.vertElaStress=vertElaStress;
    result.Degree=degreeSet(iter);
%     result.IncompressValue=transpose(cell2mat(incompressCell));
    resultSet{iter}=result;

    %%
    % plotShakedownDomain(resultSet);
    % disp({'LoadFactor:' Result.alpha;
    %       'MaxYield:'       Result.ProcessedRes.MaxYield;
    %       'EquiliViolation:'    Result.ProcessedRes.EquiViolateAB;
    %       'UAdditionConvertTolerance:' Result.ProcessedRes.UAdditionConvertTolerance;});
    
  
end
%%
% save(DataName,'resultSet','-append');
% err=-1+transpose(result.pi(:))*linConstHardEnrich(:,alphaPos); %gradL2alpha
% GurobiKKTvalidKinHard(model,linConstHardEnrich,QuadIndicesCell,ConeHardIndicesCell,result);
% rmpath(gurobiPath);


end
