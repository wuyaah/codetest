function parameters=SubExtract3DConvertProcess(dataPack,NodeLabelList,EleNodeLabelSet,NodeConstrainedDof)

global ABMatSparse CMatSparse

numElem=numel(dataPack);
numNode=numel(NodeLabelList);
dimension=dataPack{1}.dimension;
numStrComp=dataPack{1}.numStrComp; % Stress=[S11;S22;S12;S12;S23;S13];
elemCMat=cellfun(@(x)x.elemCMat,dataPack,'UniformOutput',0);
elemNGK=cellfun(@(x)x.NGK,dataPack);

% Specified Parameters
% L_MatTrans=transpose([sqrt(2) 0 0 0 0; 1/sqrt(2) sqrt(3)/sqrt(2) 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1;]);
% T_MatInv=[1 -1 0 0 0 0; 0  1 -1 0 0 0; 1 0 1 0 0 0; 0 0 0 sqrt(6) 0 0; 0 0 0 0 sqrt(6) 0; 0 0 0 0 0 sqrt(6)];   
% S2UMat=(1/sqrt(2))*L_MatTrans*T_MatInv([1,2,4:6],:);

%% DataProcessing
elemYieldCell=cellfun(@(x)x.Yield,dataPack,'UniformOutput',0); 
% EleYield=kron(cell2mat(elemYieldCell),ones(NGK,1));
% EleHardLimit=cellfun(@(x) x.HardLimit,dataPack,'UniformOutput',0); 
% EleHardLimit=cell2mat(EleHardLimit);
% EleHardLimit=kron(EleHardLimit,ones(NGK,1));
% EleRelaHard=(EleHardLimit-EleYield)./EleYield;
% ElePoisson=cellfun(@(x) x.Poisson,dataPack,'UniformOutput',0); 
% ElePoisson=cell2mat(ElePoisson);
% ElePoisson=kron(ElePoisson,ones(NGK,1));

[elemA_Mat,elemB_Mat]=cellfun(@(C,Y)VarConvert3D(C,Y),elemCMat,elemYieldCell,'UniformOutput',0);
dataPack=cellfun(@(Info,A,B)SortABMat2Elem(Info,A,B),dataPack,elemA_Mat,elemB_Mat,'UniformOutput',0);

[AMatSparseFull, VarU1EntryList]=GetConstrainedGlobalMatSparseList(dataPack,'elemA_Mat',dimension,5,EleNodeLabelSet,NodeLabelList,NodeConstrainedDof);
[BMatSparseFull, VarX1EntryList]=GetConstrainedGlobalMatSparseList(dataPack,'elemB_Mat',dimension,1,EleNodeLabelSet,NodeLabelList,NodeConstrainedDof);
[CMatSparseFull,StressEntryList]=GetConstrainedGlobalMatSparseList(dataPack,'elemCMat',dimension,6,EleNodeLabelSet,NodeLabelList,NodeConstrainedDof);
ABMatSparseFull=[AMatSparseFull,BMatSparseFull];

CMatSparseFull=mat2cell(CMatSparseFull,ones(size(CMatSparseFull,1),1),size(CMatSparseFull,2));
VarUX1EntryList=[VarU1EntryList,max(max(VarU1EntryList))+VarX1EntryList];

[rowABSparse,colABSparse]=size(ABMatSparseFull); 
parameters={numElem;numNode;elemNGK;numStrComp;rowABSparse;colABSparse;dimension};
% EleLabelList=num2str(num2str(EleLabelList));
CMatSparseFull=cell2mat(CMatSparseFull);

% % Sorting Displacement data if it has been extracted 
% if exist('EleNodeDispSet')
%    SortedNodeDispSet=permute(EleNodeDispSet,[1 3 2]);
%    SortedNodeDispSet=cell2mat(reshape(SortedNodeDispSet,[size(EleNodeDispSet,1) size(EleNodeDispSet,2)*size(EleNodeDispSet,3)]));
%    SortedNodeDispSet=mat2cell(SortedNodeDispSet,ones(ElemNum,1),size(EleNodeDispSet,3)*ones(size(EleNodeDispSet,2),1));
%    EleNodeLabelSet=cell2mat(EleNodeLabelSet);
%    NodeLabelListCell=num2cell(NodeLabelList);
%    PosList=cellfun(@(x)find(EleNodeLabelSet==x),NodeLabelListCell,'UniformOutput',0);
%    SortedDisp=cellfun(@(x)transpose(SortedNodeDispSet{x(1)}),PosList,'UniformOutput',0);
%    SortedDisp=cell2mat(SortedDisp);
%    clear PosList NodeLabelListCell SortedNodeDispSet
% else
%     SortedDisp=[];
% end

% % Sorting Strain data if it has been extracted
% if exist('EleGaussStrainSet')
%     SortedStrain=permute(EleGaussStrainSet,[1 3 2]);
%     SortedStrain=cell2mat(reshape(SortedStrain,[size(SortedStrain,1) size(SortedStrain,2)*size(SortedStrain,3)]));
%     SortedStrain=mat2cell(SortedStrain,ones(ElemNum,1),size(SortedStrain,3)*size(SortedStrain,2));
%     SortedStrain=cellfun(@(x) reshape(transpose(x),6,size(x,2)/6),SortedStrain,'UniformOutput',0);
%     SortedStrain=cellfun(@(x) [x(1,:); x(2,:); x(3,:); x(4,:); x(6,:); x(5,:)],SortedStrain,'UniformOutput',0);
%     SortedStrain=cellfun(@(x) reshape(x,size(x,1)*size(x,2),1),SortedStrain,'UniformOutput',0);
%     SortedStrain=cell2mat(SortedStrain);
% else
%     SortedStrain=[];
% end

% Convert AB Matrix to real sparse form
[RowNum,MaxRelevantColumnNum]=size(ABMatSparseFull);
RowNumVec=kron((1:RowNum)',ones(MaxRelevantColumnNum,1));
ABMatEntryReform=mat2cell(VarUX1EntryList,ones(RowNum,1),MaxRelevantColumnNum);
ABMatEntryReform=cellfun(@(x)transpose(x),ABMatEntryReform,'UniformOutput',0);
ColumnNumVec=cell2mat(ABMatEntryReform);

ABMatReform=mat2cell(ABMatSparseFull,ones(RowNum,1),MaxRelevantColumnNum);
ABMatReform=cellfun(@(x)transpose(x),ABMatReform,'UniformOutput',0);
ValueVec=cell2mat(ABMatReform);

ABMatSparse=sparse(RowNumVec,ColumnNumVec,ValueVec);
ABMatSparse=[ABMatSparse(:,1:(numStrComp-1)*sum(elemNGK)),ABMatSparse(:,(numStrComp-1)*sum(elemNGK)+2:end-1)];  %The idle columns are already removed
clear RowNumVec ColumnNumVec ValueVec

% Convert C Matrix to real sparse form
[RowNum,MaxRelevantColumnNum]=size(CMatSparseFull);
RowNumVec=kron((1:1:RowNum)',ones(MaxRelevantColumnNum,1));
CMatEntryReform=mat2cell(StressEntryList,ones(RowNum,1),MaxRelevantColumnNum);
CMatEntryReform=cellfun(@(x)transpose(x),CMatEntryReform,'UniformOutput',0);
ColumnNumVec=cell2mat(CMatEntryReform);

CMatReform=mat2cell(CMatSparseFull,ones(RowNum,1),MaxRelevantColumnNum);
CMatReform=cellfun(@(x)transpose(x),CMatReform,'UniformOutput',0);
ValueVec=cell2mat(CMatReform);

CMatSparse=sparse(RowNumVec,ColumnNumVec,ValueVec);
CMatSparse=CMatSparse(:,1:end-1);
clear NumStrComp LineNumVec ColumnNumVec ValueVec 
%%% ori
% save(['C:\CodeForWork\Output\Extracted' DataName],'dataPack','EleStressData','EleYield','NodeLabelList','EleBMat', ...
%       'NodeConstrainedDof1','NodeConstrainedDof2','NodeConstrainedDof3','GaussVol','Parameters', ...
%       'ABMatSparse','S2UMat','CMatSparse','AlternativeLoadCaseStressSet','SortedDisp','SortedStrain','-v7.3');

% save(['C:\CodeForWork\Output\Extracted' DataName],'dataPack','EleYield','NodeLabelList','EleBMat', ...
%       'NodeConstrainedDof','GaussVol','Parameters', ...
%       'ABMatSparse','S2UMat','CMatSparse','AlternativeLoadCaseStressSet','SortedDisp','SortedStrain','-v7.3');
% 

end
