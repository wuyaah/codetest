%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main Func %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;        % 建议不要用 clear all（效率低并且会清除 breakpoint）
close all;
% warningInfo = warning('off', 'all');
outputPath = 'D:\outputTopoRSDMS';
if ~isfolder(outputPath)       % 若文件夹不存在则创建
    mkdir(outputPath);
end
addpath(genpath(pwd));
addpath(outputPath);
tolTime=tic;
%%
dataDir=('D:\dateDir\abqDir\Cantilever2D');
dataName=('Cantilever2DLoad22');
alternativeLoadCases={};

addpath(dataDir);
eval(['mat_' dataName]);

%% matInfo and Optimization Params
% oriRelDens=VolumeFrac;
oriRelDens=1;
% oriRelDens=1;
degree=45;
if isempty(alternativeLoadCases); degree=0; end
delta=0.01; maxIter=50; minRelDens=1e-3; Penal=3; YPenal=1; UPenal=1; %opt params (default: YPenal=1)
matInfo={E.*[1,1e-8]; Yield.*[1,1e-8]; Limit.*[1,1e-8]; mu; Penal; YPenal; UPenal}; %material
iskinemHardening=0; % "0": EPP; "1": KH;  (perfectly elastoplastic)
updateScheme=1;     % "1": MMA "2": OC
vertNorm=0;
numVert=2;
optType=1;          % "1": SD  "2": SE "3": VM

%% Initialize model
[elemGaussStrainSet,elemNodeDispSet,nodeConstrained]=LoadModelData(dataName);
[dataPack,elemLabelList,elemNodeLabelSet,nodeLabelList,frozenElem]=OriModelDataDamageFreeNew(dataName,alternativeLoadCases);
[nodeForceVec,constrainedDoFVec]=OriModelNodeData(dataName,elemNodeLabelSet,alternativeLoadCases);
dataPack=cellfun(@(x)getElemCMatInfo(x,matInfo),dataPack,'UniformOutput',0);
dataPack=cellfun(@(x)getGlobPosElemNodeDoF(x,nodeLabelList),dataPack,'UniformOutput',0);
[relDensSet,paramsIter]=OriElemRelDens(dataPack,frozenElem,oriRelDens,VolumeFrac,minRelDens);
[dataPack,elemCentroidSet]=getElemCentroid(dataPack);
% relDensSet=LoadDensHistory; %% 检验结果
elaMatrix=getElaMatrixFunc(dataPack{1}.dimension,dataPack{1}.Poisson);
dataPack=cellfun(@(x)getElemsINGivenRadius(x,elemCentroidSet,FilterRadius),dataPack,'UniformOutput',0);
dataPack=cellfun(@(x,d)SortElemRelDens(x,d),dataPack,num2cell(relDensSet),'UniformOutput',0);
dataPack=cellfun(@(x)getElemKMatrix(x,elaMatrix),dataPack,'UniformOutput',0);

error_Matrix=ValidForMatrix(dataPack,elemGaussStrainSet); %B矩阵验证

clear elemGaussStrainSet VolumeFrac oriRelDens minRelDens setFrozen
clear E mu Yield Limit Penal YPenal

%% Ori Simp Optimization
objSet=[]; volumeFracAttainedSet=[]; relDensHistory=[]; elaLimitSet=[];
alphaScaling=3; %alpha=alphaFactor*Yield/Sigma_max
numTerms=3;     %number of the Fourier series
numLoadCase=size(nodeForceVec,2);
params=getNumSetRSDMS(dataPack,alphaScaling,numLoadCase,numVert,vertNorm,numTerms);
tolSet.innerTol=1e-3; %1e-4~1e-3
tolSet.outerTol=1e-2; %1e-3~1e-2

%% iterate
disp("Performing shakedown oriented TO via RSDM ...");
getRunTime(tolTime);
iter=0;

%%
while true
    
    t=tic; %clc
    iter=iter+1;% iter=iter-1
    oldRelDensSet=relDensSet;

    %% sensitivity analysis
    dataPack=cellfun(@(x)getIntMatInfo(x,matInfo),dataPack,'UniformOutput',0);
    globKMat=getGlobStiffMatFunc(dataPack,numel(nodeLabelList));
    % ValidateFEARes(dataName,dataPack,globKMat,elaMatrix,nodeForceVec,elemNodeLabelSet,nodeLabelList,constrainedDoFVec);
    nodeDispVecCell=getNodeDispVecNew(globKMat,nodeForceVec,constrainedDoFVec);
    dataPack=cellfun(@(x)SortGlobDisp2Ele(x,nodeDispVecCell),dataPack,'UniformOutput',0);
    dataPack=cellfun(@(x)CalGaussStress(x,elaMatrix),dataPack,'UniformOutput',0);

    %% 计算结果对比: Gurobi
    %%1：result_Quad
    % [elaStress,NGYield,NGRelHard]=subExtractStress(dataPack);
    % parameters=ExtractConvert(dataPack,elemNodeLabelSet,nodeLabelList,nodeConstrained);
    % result_Quad=SubSAGurobiCSCone(elaStress,NGYield,parameters,numVert,degree);
    % result2=result_Quad{end};
    % dataPack2=SortLambda2Ele(dataPack,result2);
    % PlaDissErr2=Paper_GetValidForErr_CS(dataPack2,result2);
    % lambda2=cellfun(@(x)x.lambda,dataPack2,'UniformOutput',0);
    % lambda2=cell2mat(lambda2);

    %%2：result_Cone
    % [elaStress,NGYield,NGRelHard]=subExtractStress(dataPack);
    % parameters=ExtractConvert(dataPack,elemNodeLabelSet,nodeLabelList,nodeConstrained);
    % result_Cone=ShakedownCodeNorm(iskinemHardening,elaStress,NGYield,NGRelHard,parameters,degree,numVert,vertNorm);
    % result3=result_Cone{end};
    % dataPack3=SortLambda2Ele(dataPack,result3);
    % PlaDissErr3=Paper_GetValidForErr(dataPack3,result3);
    % lambda3=cellfun(@(x)x.lambda,dataPack3,'UniformOutput',0);
    % lambda3=cell2mat(lambda3);

    %%
    % elemYield=cellfun(@(x)x.Yield,dataPack);
    % elaLimitSet2(iter,1)=getElaLimitFactor(result2,elemYield); %#ok<SAGROW> 
    % SD_Elaratio=result2.alpha/elaLimitSet2;
    % sensitivity2=sensAnalysisNorm(iskinemHardening,dataPack,result2,elemYield);
    % sensitivity2=subAlternativeSensNorm(sensitivity2); 
    % stressTrial=CalStressTrial(dataPack,result2,elaMatrix);
    
    %% ============================== RSDMS ===============================
    result_RSDMS=SubCodeNormRSDMS(dataPack,globKMat,nodeLabelList,constrainedDoFVec,elaMatrix,params,degree,tolSet);
    result=result_RSDMS{end};
    [lambda,deltaLambda]=Paper_SortLambda2Ele(result);
    % deltaLambdaMat=cell2mat(deltaLambda);
    % PlaDissErr=validateKKTCondition(result,deltaLambda);

    %%
    elemYield=cellfun(@(x)x.Yield,dataPack);
    elaLimitSet(iter,1)=getElaLimitFactor(result,elemYield); %#ok<SAGROW>
    %% 屈服强度项敏感度+弹性项
    % dataPack=CalGradYieldFunc2Sigma(dataPack,lambda,result);
    % elemArtiDisp2=ArtiClobFSens_Paper(dataPack,globKMat,elaMatrix);
    % sensEla2Dens2=cellfun(@(x,d)Paper_SubArtiClobFSens(x,d,elaMatrix,result.vertWeight),dataPack,elemArtiDisp2);
    % sensEla2Dens2=result.alpha.*sensEla2Dens2;
    % sensYield2Dens=cellfun(@(x,l)sensAnalysisBasedOnRSDMS(x,l),dataPack,lambda);
    % sensitivity=sensEla2Dens2+sensYield2Dens;
    % sensitivity=sensitivity./max(abs(sensitivity));
    
    %% 屈服程度+屈服强度项敏感度
    YieldDeg=-mean(result.ProcessedRes.Yield,2);
    YieldDeg=mat2cell(YieldDeg,dataPack{1}.NGK*ones(numel(dataPack),1),1);
    YieldDeg=cellfun(@(x,y)x.GaussVol*y,dataPack,YieldDeg);
    sensitivity1=YieldDeg./max(abs(YieldDeg));
    % sensitivity=subAlternativeSensNorm(sensitivity);
    sensitivity2=cellfun(@(x,l)sensAnalysisBasedOnRSDMS(x,l),dataPack,lambda);
    sensitivity2=sensitivity2./max(abs(sensitivity2));
    sensitivity=0.01*sensitivity1+sensitivity2;

    %% 屈服程度+屈服强度项敏感度+弹性项
    % YieldDeg=-mean(result.ProcessedRes.Yield,2);
    % YieldDeg=mat2cell(YieldDeg,dataPack{1}.NGK*ones(numel(dataPack),1),1);
    % YieldDeg=cellfun(@(x,y)x.GaussVol*y,dataPack,YieldDeg);
    % sensitivity1=YieldDeg./max(abs(YieldDeg));
    % dataPack=CalGradYieldFunc2Sigma(dataPack,lambda,result);
    % % elemArtiDisp2=ArtiClobFSens(dataPack,globKMat,elaMatrix);
    % elemArtiDisp2=ArtiClobFSens_Paper(dataPack,globKMat,elaMatrix);
    % sensEla2Dens2=cellfun(@(x,d)Paper_SubArtiClobFSens(x,d,elaMatrix,result.vertWeight),dataPack,elemArtiDisp2);
    % sensEla2Dens2=result.alpha.*sensEla2Dens2;
    % sensYield2Dens=cellfun(@(x,l)sensAnalysisBasedOnRSDMS(x,l),dataPack,lambda);
    % sensitivity2=sensEla2Dens2+sensYield2Dens;
    % sensitivity2=sensitivity2./max(abs(sensitivity2));
    % sensitivity=0.01*sensitivity1+sensitivity2;

    %% 屈服强度项敏感度
    % sensitivity=cellfun(@(x,l)sensAnalysisBasedOnRSDMS(x,l),dataPack,lambda);
    % sensitivity=sensitivity./max(abs(sensitivity));

    %%
    obj=result.alpha;

    %% Update Design Variables
    % sensitivity=subAlternativeSensNorm(sensitivity);
    sensitivity=SetFrozenSens(sensitivity,frozenElem,paramsIter);
    sensitivity=GetFilteredSens(dataPack,sensitivity,FilterRadius);
    [dataPack,volumeFracAttained,relDensSet,paramsIter]=updateSchemeNorm(dataPack,frozenElem,sensitivity,paramsIter,updateScheme,iter);
    % [dataPack,volumeFracAttained,relDensSet,paramsIter]=updateDensVal_RSDM(dataPack,frozenElem,sensitivity,paramsIter,updateScheme,iter);
    volumeFracAttainedSet(iter,1)=volumeFracAttained;
    relDensHistory(:,iter)=relDensSet;
    objSet(iter,1)=obj;
    
    %% PostProcessing/Termination
    change=max(abs(relDensSet-oldRelDensSet)); %GetRelDensNormfit(relDensSet,oldRelDensSet);
    iterationInfo(volumeFracAttainedSet,elaLimitSet,objSet,t,change);
    visualizeTopoRes(elemCentroidSet,relDensSet,dataPack{1}.dimension);
    subMakeGif(fullfile(outputPath,['gifTO_',dataName,'_',num2str(iskinemHardening),'.gif']),iter);
    clear obj oldRelDensSet volumeFracAttained
    %=======================================================
    if (change<delta)||(iter>=maxIter)
        break;
    end
end
%%
print('-dtiffn','-r300',fullfile(outputPath,['TO_',dataName,'_',num2str(iskinemHardening),'.tiff']));
plotTopoCurve(volumeFracAttainedSet,optType,objSet,elaLimitSet); % ObjSet./ElaLimitSet
% print('-dtiffn','-r300',fullfile(outputPath,['obj_',dataName,'_',num2str(optType),num2str(iskinemHardening),'.tiff']));
getRunTime(tolTime);
topoResultSave;
