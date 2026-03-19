function result=GurobiPostProcess(result,LinConst,EleYield,VertElaStress,numSet)

global CMatSparse
% this function interprets the result
% Result: Result of Gurobi calculation
% GaussStrength: [SY@Int1 SY@Int2 ... SY@Int_m]
% OriElaStress: Original Elastic Stress Vector
% ConvertFunc: The function handle used for the conversion.
% OutPut "UAdditionConvertTolerance" measures the discrepency between
% S2U(\alpha*ElaStress@Vert_n+\rho)-U2  which is the error caused during
% the conversion between variables
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
% Var=[U1; X; alpha; U2; U3; U4; .... Beta]
%%
numStrComp=numSet.numStrComp;
numTotalGauss=numSet.NumTotalGauss;
numU1Comp=numSet.NumU1Comp;
numVert=numSet.numVert;
numU1=numSet.NumU1;
numX=numSet.NumX1;

result.Time=datetime("now"); 
result.numVert=numVert; 
% result.x=result.x(1:size(LinConst,2));
result.alpha=abs(result.objval);
result=rmfield(result,'objval');
result=rmfield(result,'versioninfo');

UXList=result.x;
U1=UXList(1:numU1);
X=reshape(UXList(numU1+1:numU1+numX),1,[]);
LoadFactor=result.alpha;
vertUlist=U1;

% Choose the conversion function depending on the element type
% [UConvertSigma,SigmaConvertU,EffCalFunc]=getConvertFunc(numStrComp);
[EffCalFunc,~,UConvertSigma,SigmaConvertU]=getConvertFunc(numStrComp);

U1Cell=mat2cell(U1,numU1Comp*ones(numTotalGauss,1),1);
U1Regular=reshape(U1,numU1Comp,numTotalGauss); 
U1XRegular=[U1Regular;X];
U1XRegular=mat2cell(U1XRegular,size(U1XRegular,1),ones(size(U1XRegular,2),1));
TotalStressVert1Cell=cellfun(@(x,y)UConvertSigma(x)*y,num2cell(EleYield),transpose(U1XRegular),'UniformOutput',0);
EffTotalStressVert1=cellfun(@(x)EffCalFunc(x),TotalStressVert1Cell);
StressVert1=cell2mat(TotalStressVert1Cell);
ResidualStress=StressVert1-LoadFactor.*VertElaStress(:,1);
ResidualStressCell=mat2cell(ResidualStress,numStrComp*ones(numTotalGauss,1),1);
EffectiveResidual=cellfun(@(x)EffCalFunc(x),ResidualStressCell);
EquiViolateAB=max(abs(LinConst*result.x));
EquiViolateC=max(abs(CMatSparse*ResidualStress));
Yield=cellfun(@(x)norm(x),U1Cell);
UAdditionConvertTolerance=0;
TotalStressAdditional=[];

if numVert>1
    UAdditional=UXList(numU1+numX+2:end-1);          % Get U for all other load verteces (2...8)
    UAdditional=reshape(UAdditional,numU1,numVert-1);
    UAdditionalCell=mat2cell(UAdditional,numU1Comp*ones(numTotalGauss,1),ones(numVert-1,1));
    YieldAdditional=cellfun(@(x)norm(x),UAdditionalCell);
    Yield=[Yield,YieldAdditional];
    TotalStressAdditional=LoadFactor*VertElaStress(:,2:end)+repmat(ResidualStress,1,size(2:numVert,2));
    TotalStressAdditionalCell=mat2cell(TotalStressAdditional,numStrComp*ones(numTotalGauss,1),ones(numVert-1,1));
    ConvertedUXAdditionalCell=cellfun(@(x,y)SigmaConvertU(x)*y,...
         num2cell(repmat(EleYield,1,size(TotalStressAdditionalCell,2))),TotalStressAdditionalCell,'UniformOutput',0);
    if numStrComp~=3
       ConvertedUAdditional=cell2mat(cellfun(@(x)x(1:end-1),ConvertedUXAdditionalCell,'UniformOutput',0));
    else
        ConvertedUAdditional=cell2mat(ConvertedUXAdditionalCell);
    end
    UAdditionConvertTolerance=max(max(abs(ConvertedUAdditional-UAdditional)));
    vertUlist=[vertUlist,UAdditional];
end

TotalStress=[StressVert1,TotalStressAdditional];

% Switch between load cases
ProcessedRes.TotalStress=TotalStress;
ProcessedRes.ResidualStress=ResidualStress;
ProcessedRes.EquiViolateAB=EquiViolateAB;
ProcessedRes.EquiViolateC=EquiViolateC;
ProcessedRes.EffTotalStress=EffTotalStressVert1;
ProcessedRes.EffResidualStress=EffectiveResidual;
ProcessedRes.UAdditionConvertTolerance=UAdditionConvertTolerance;
ProcessedRes.MaxYield=max(max(Yield));
ProcessedRes.Yield=Yield;

result.ProcessedRes=ProcessedRes;
result.vertUlist=vertUlist;

end

