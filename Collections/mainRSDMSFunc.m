%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main Func %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;        % 建议不要用 clear all（效率低并且会清除 breakpoint）
close all;
% warningInfo = warning('off', 'all');
outputPath = 'D:\outputTopoRSDMS';
if ~isfolder(outputPath)       % 若文件夹不存在则创建
    mkdir(outputPath);
end
addpath(outputPath);
addpath(genpath(pwd));
tolTime=tic;
%%
%%% Set data directory
% dataDir=('D:\dateDir\abqDir\Cantilever2D');
% dataName=('Cantilever2DLoad22');
% alternativeLoadCases={};

dataDir=('D:\dateDir\abqDir\HoledPlate\HoledPlate_2Layer');
dataName=('HoledPlateLoad11');
alternativeLoadCases={};
% alternativeLoadCases={'HoledPlateLoad22'};
addpath(dataDir);

%%% load angle in degrees (scalar or vector)
numVert=2;  % defines loading domain or path
vertNorm=0;
% degree=[0;30;45;60;90];
% degree=(0:30:360)';
degree=45;
if isempty(alternativeLoadCases); degree=0; end

%%% Material parameter settings: [E, Yield, Limit, μ, Penal, YPenal, UPenal]
eval(['mat_' dataName]);%known material parameters
Penal=3; YPenal=1; UPenal=1;
matInfo={E.*[1,1e-8]; Yield.*[1,1e-8]; Limit.*[1,1e-8]; mu; Penal; YPenal; UPenal}; %Material Info

clear E Yield Limit mu Penal YPenal UPenal FilterRadius VolumeFrac

%% Initialize model
[dataPack,elemNodeLabelSet,nodeLabelList]=OriModelData(dataName,alternativeLoadCases);
[nodeForceVec,constrainedDoFVec]=OriModelNodeData(dataName,elemNodeLabelSet,alternativeLoadCases);
[dataPack,elemCentroidSet]=getElemCentroid(dataPack);
dataPack=cellfun(@(x)getElemCMatInfo(x,matInfo),dataPack,'UniformOutput',0); %/BMat/CMat
dataPack=cellfun(@(x)getGlobPosElemNodeDoF(x,nodeLabelList),dataPack,'UniformOutput',0);
dataPack=cellfun(@(x)SortElemRelDens(x,1),dataPack,'UniformOutput',0); %original density: 1
dataPack=cellfun(@(x)getIntMatInfo(x,matInfo),dataPack,'UniformOutput',0);
elaMatrix=getElaMatrixFunc(dataPack{1}.dimension,dataPack{1}.Poisson);
dataPack=cellfun(@(x)getElemKMatrix(x,elaMatrix),dataPack,'UniformOutput',0);
globKMat=getGlobStiffMatFunc(dataPack,numel(nodeLabelList));
% globKMat=readAbaqusStiffnessMatrix(dataName,'Load11_STIF1.mtx');

%% Valid the results from FE calculattion between Matlab and abaqus
% ValidateFEARes(dataName,dataPack,globKMat,elaMatrix,nodeForceVec,elemNodeLabelSet,nodeLabelList,constrainedDoFVec);
% nodeDispVecCell=getNodeDispVecNew(globKMat,nodeForceVec,constrainedDoFVec);
% dataPack=cellfun(@(x)SortGlobDisp2Ele(x,nodeDispVecCell),dataPack,'UniformOutput',0);
% dataPack=cellfun(@(x)CalGaussStress(x,elaMatrix),dataPack,'UniformOutput',0);% ValidateFEARes(dataName,dataPack,globKMat,nodeForceVec,elemNodeLabelSet,nodeLabelList,constrainedDoFVec);
% nodeDispVecCell2=getNodeDispVecNew(globKMat,nodeForceVec,constrainedDoFVec);
% plot(abs(nodeDispVecCell{1}-nodeDispVecCell2{1}));

%% SOCP-Gurobi
[elaStress,elemYield,elemRelHard]=subExtractStress(dataPack);
[~,~,nodeConstrained]=LoadModelData(dataName);
parameters=ExtractConvert(dataPack,elemNodeLabelSet,nodeLabelList,nodeConstrained);
result_Gurobi=ShakedownCodeNorm(0,elaStress,elemYield,elemRelHard,parameters,degree,numVert,vertNorm);
result_Gurobi{1}.alpha;
% plotShakedownDomain(result_Gurobi,matInfo{2}(1));
% result=result_Gurobi{end}; 

%% RSDM-S
% numVert=4;
pause(0.5); clc; close
alphaScaling=3; %alpha=alphaFactor*Yield/Sigma_max
numTerms=3; %number of the Fourier series
numLoadCase=size(nodeForceVec,2);
numSet=getNumSetRSDMS(dataPack,alphaScaling,numLoadCase,numVert,vertNorm,numTerms);
tolSet.innerTol=1e-3; %1e-4~1e-3
tolSet.outerTol=1e-2; %1e-3~1e-2
% ==========================
% result=SubCodeNormRSDM(dataPack,globKMat,nodeLabelList,constrainedDoFVec,elaMatrix,params,degree,Tol);
result_RSDMS=SubCodeNormRSDMS(dataPack,globKMat,nodeLabelList,constrainedDoFVec,elaMatrix,numSet,degree,tolSet);
% result_RSDMS=SubCodeNormRSDMS_Vert(dataPack,globKMat,nodeLabelList,constrainedDoFVec,elaMatrix,numSet,degree,tolSet);
% plotShakedownDomain(result_RSDMS,matInfo{2}(1));
% plotRSDMSDescending(result_RSDMS);
% ==========================
% save(dataName,'result_Gurobi','-append');
% save(dataName,'result_RSDMS','-append');
result=result_RSDMS{end};
getRunTime(tolTime);
rmpath(dataDir);
