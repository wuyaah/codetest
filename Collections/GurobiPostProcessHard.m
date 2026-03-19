function result=GurobiPostProcessHard(result,linConstHardEnrich,ngYield,ngRelaHard,vertElaStress,indexHardGauss,ConstElaStress)
global NumSet CMatSparse
% ResultSum=postProcessingHard(Result,U1X,GaussStrength,ConvertHandleSet)
% this function interprets the result
% Result: Result of Gurobi calculation
% GaussStrength: [SY@Int1 SY@Int2 ... SY@Int_m]
% OriElaStress: Original Elastic Stress Vector
% IndexHardGauss: indices of Gaussian points have hardening behaviors
% ConvertFunc: The function handle used for the conversion.
% OutPut "TotalStressAdditional" is obatained by \alpha*SigmaEla+ResidualStress 
% OutPut "UAdditionConvertTolerance" measures the discrepency between
% S2U(\alpha*ElaStress@Vert_n+\rho)-U2  which is the error caused during the conversion between variables
%-------------------------------------------------------------------------------
% The numbering order of indices:
% 2 Vertices:   %4 Vertices:    %8 Vertices:
%       / 1     %  4------ 1    %   6-------1
%     /         %  |       |    %  /|      /|
% 2 /           %  |       |    % 7-|----8/ |
%               %  2 ------3    % | |     | |
                                % |/2 --- | /4
                                % 3------ 5
%-------------------------------------------------------------------------------
if nargin<7
    ConstElaStress=sparse(size(vertElaStress,1),1);
end
numStrComp=NumSet.numStrComp;
numTotalGauss=NumSet.NumTotalGauss;
numHardGauss=NumSet.NumHardGauss;
numU1Comp=NumSet.NumU1Comp;
numVert=NumSet.NumVert;
numU1=NumSet.NumU1;
numX1=NumSet.NumX1;

result.Time=clock; 
result.NumVert=numVert; 
result.alpha=abs(result.objval);
result=rmfield(result,'objval');
result=rmfield(result,'versioninfo');
result.NumVert=numVert;
loadFactor=result.alpha;
UXList=result.x;

indexHardGauss=indexHardGauss';
ngYieldHard=ngYield(indexHardGauss);

% Choose the conversion function depending on the element type
[effCalFunc,~,UX2Sigma,sigma2UX]=getConvertFunc(numStrComp); 

U1=UXList(1:numU1);
U1Cell=mat2cell(U1,numU1Comp*ones(numTotalGauss,1),1);
U1Regular=reshape(U1,numU1Comp,numTotalGauss); 
X1=UXList(numU1+1:numU1+numX1);
U1XRegular=[U1Regular;X1'];
U1XRegular=mat2cell(U1XRegular,size(U1XRegular,1),ones(size(U1XRegular,2),1));
U1XRegular=transpose(U1XRegular);
VertUlist=U1;

%Process Back Stress
UPi=UXList(NumSet.totalVarPerfect+1:NumSet.totalVarPerfect+(numU1Comp*numHardGauss));
XPi=UXList(NumSet.totalVarPerfect+1+(numU1Comp*numHardGauss):end-1);
UPiCell=mat2cell(UPi,numU1Comp*ones(numHardGauss,1),1);
UPiRegular=reshape(UPi,numU1Comp,numHardGauss);
UPiRegular=UPiRegular*diag(sparse(ngRelaHard)); %恢复由于简化背应力Cone约束的右端项为“1”，导致的BackUX造成的缩放
BackUXRegular=[UPiRegular;XPi'];
BackUXRegular=mat2cell(BackUXRegular,size(BackUXRegular,1),ones(size(BackUXRegular,2),1));
BackUXRegular=transpose(BackUXRegular);

%Difference of stress
DiffStressVert1Cell=cellfun(@(x,y)UX2Sigma(x)*y,num2cell(ngYield),U1XRegular,'UniformOutput',0);    % DiffStress: Total Stress- Back Stress 

%Back stress
BackStressCell=cellfun(@(x,y)UX2Sigma(x)*y,num2cell(ngYieldHard),BackUXRegular,'UniformOutput',0);
MinusBackStressCell=cellfun(@(x)-1.*x,BackStressCell,'UniformOutput',0); %-1*BackStress           
EffBackStress=cellfun(@(x)effCalFunc(x),BackStressCell);

%Total stress
totalStressVert1Cell=getTotalStress(DiffStressVert1Cell,BackStressCell,indexHardGauss);
EffTotalStressVert1=cellfun(@(x)effCalFunc(x),totalStressVert1Cell);

%Fill back stress back to the global list of gaussian points
BackStressBlock=mat2cell(sparse(numStrComp,1),numStrComp,1); 
GlobalBackStressCell=repmat(BackStressBlock,numTotalGauss,1); %list of back stresses in the entire model, includes also non-hardening elements
GlobalBackStressCell(indexHardGauss)=BackStressCell;
GlobalEffBackStress=zeros(numTotalGauss,1);
GlobalEffBackStress(indexHardGauss)=EffBackStress;

totalStressVert1=cell2mat(totalStressVert1Cell);
ResidualStress=totalStressVert1-loadFactor.*vertElaStress(:,1)-ConstElaStress;
ResidualStressCell=mat2cell(ResidualStress,numStrComp*ones(numTotalGauss,1),1);
EffectiveResidual=cellfun(@(x)effCalFunc(x),ResidualStressCell);

rhsInConst=CMatSparse*ConstElaStress;
EquilViolateAB=max(abs(linConstHardEnrich*result.x-[rhsInConst;sparse(size(linConstHardEnrich,1)-size(rhsInConst,1),1)]));
EquiViolateC=max(abs(CMatSparse*ResidualStress));
%
BackYield=cellfun(@(x)norm(x),UPiCell);
Yield=cellfun(@(x)norm(x),U1Cell);
% MaxEffBack=max(cellfun(@(x)x'*x,UPiCell));
% MaxYield=max(cellfun(@(x)x'*x,U1Cell));
%
totalStressAdditional=[];
UAdditionConvertTolerance=0;

if numVert>1
    %Var=[U1; X; -alpha; U2; U3; U4; .... Beta; U_Pi(Binder); X_Pi(); Beta2] 
    UAdditional=UXList(numU1+numX1+2:NumSet.totalVarPerfect-1);  % get U for all other load verteces (2...8)
    UAdditional=reshape(UAdditional,numU1,numVert-1);
    UAdditionalCell=mat2cell(UAdditional,numU1Comp*ones(numTotalGauss,1),ones(numVert-1,1));
    YieldAdditional=cellfun(@(x)norm(x),UAdditionalCell);
    Yield=[Yield,YieldAdditional];
    totalStressAdditional=loadFactor.*vertElaStress(:,2:numVert)+repmat(ResidualStress+ConstElaStress,1,numVert-1);
    totalStressAdditionalCell=mat2cell(totalStressAdditional,numStrComp*ones(numTotalGauss,1),ones(numVert-1,1));   %alpha*ElaStress+Residual
    diffStressAdditionalCell=totalStressAdditionalCell;
    for m=1:(numVert-1)
        diffStressAdditionalCell(:,m)=getTotalStress(totalStressAdditionalCell(:,m),MinusBackStressCell,indexHardGauss);          %Obtain the difference of stress by extracting back stress from the total stress 
    end
    additionalVertYield=num2cell(repmat(ngYield,1,numVert-1));
    ConvertedUXAdditionalCell=cellfun(@(x,y)sigma2UX(x)*y,additionalVertYield,diffStressAdditionalCell,'UniformOutput',0);
    if numStrComp==3 % CPS case
        ConvertedUAdditional=cell2mat(ConvertedUXAdditionalCell);
    else
        ConvertedUAdditional=cell2mat(cellfun(@(x)x(1:end-1),ConvertedUXAdditionalCell,'UniformOutput',0));
    end
     UAdditionConvertTolerance=max(max(abs(ConvertedUAdditional-UAdditional)));
     VertUlist=[VertUlist,UAdditional];
end

totalStress=[totalStressVert1,totalStressAdditional];
VertUlist=[VertUlist,UPi];

%Switch between load cases
ProcessedRes.TotalStress=totalStress;
ProcessedRes.ResidualStress=ResidualStress;
ProcessedRes.BackStress=cell2mat(GlobalBackStressCell);

ProcessedRes.UAdditionConvertTolerance=UAdditionConvertTolerance;
ProcessedRes.EquilViolateAB=EquilViolateAB;
ProcessedRes.EquiViolateC=EquiViolateC;

ProcessedRes.BackYiled=BackYield;
ProcessedRes.Yield=Yield;
ProcessedRes.MaxBackYiled=max(BackYield);
ProcessedRes.MaxYield=max(max(Yield));

ProcessedRes.EffTotalStress=EffTotalStressVert1;
ProcessedRes.EffResidualStress=EffectiveResidual;
ProcessedRes.EffBackStress=GlobalEffBackStress;

% result.BackUXRegular=oriBackUXRegular;
result.ProcessedRes=ProcessedRes;
result.VertUlist=VertUlist;

end


