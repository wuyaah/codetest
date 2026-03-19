function parameters=SubExtractCPSConvertProcess(dataPack,NodeLabelList,EleNodeLabelSet,NodeConstrained)

global ABMatSparse CMatSparse

numElem=numel(dataPack);
numNode=numel(NodeLabelList);
dimension=dataPack{1}.dimension;
numStrComp=dataPack{1}.numStrComp; % Stress=[S11;S22;S12];
elemCMat=cellfun(@(x)x.elemCMat,dataPack,'UniformOutput',0);
elemNGK=cellfun(@(x)x.NGK,dataPack);

%%% Specified Parameters
% InvL_MatTrans=[1, sqrt(3)/3, 0; 0, sqrt(3)*2/3, 0; 0, 0, sqrt(3)/3];
% S2UMat=inv(InvL_MatTrans);

%% DataProcessing
elemYieldCell=cellfun(@(x)x.Yield,dataPack,'UniformOutput',0);
elemAB_Mat=cellfun(@(C,Y)VarConvertCPS(C,Y),elemCMat,elemYieldCell,'UniformOutput',0);
dataPack=cellfun(@(x,AB)SortABMat2Elem(x,AB),dataPack,elemAB_Mat,'UniformOutput',0);

[ABMatSparseFull,VarUX1EntryList]=GetConstrainedGlobalMatSparseList(dataPack,'elemAB_Mat',dimension,3,EleNodeLabelSet,NodeLabelList,NodeConstrained);
[CMatSparseFull,StressEntryList]=GetConstrainedGlobalMatSparseList(dataPack,'elemCMat',dimension,3,EleNodeLabelSet,NodeLabelList,NodeConstrained);

CMatSparseFull=mat2cell(CMatSparseFull,ones(size(CMatSparseFull,1),1),size(CMatSparseFull,2));
CMatSparseFull=cell2mat(CMatSparseFull);
[rowABSparse,colABSparse]=size(ABMatSparseFull); 
parameters={numElem;numNode;elemNGK;numStrComp;rowABSparse;colABSparse;dimension};
% EleLabelList=num2str(num2str(EleLabelList));

% Convert AB Matrix to real sparse form
[rowNum,maxRelevantColumnNum]=size(ABMatSparseFull);
RowNumVec=kron((1:rowNum)',ones(maxRelevantColumnNum,1));
ABMatEntryReform=mat2cell(VarUX1EntryList,ones(rowNum,1),maxRelevantColumnNum);
ABMatEntryReform=cellfun(@(x)transpose(x),ABMatEntryReform,'UniformOutput',0);
ColumnNumVec=cell2mat(ABMatEntryReform);

ABMatReform=mat2cell(ABMatSparseFull,ones(rowNum,1),maxRelevantColumnNum);
ABMatReform=cellfun(@(x)transpose(x),ABMatReform,'UniformOutput',0);
ValueVec=cell2mat(ABMatReform);

ABMatSparse=sparse(RowNumVec,ColumnNumVec,ValueVec);
ABMatSparse=ABMatSparse(:,1:end-1);
clear RowNumVec ColumnNumVec ValueVec

% Convert C Matrix to real sparse form
[rowNum,maxRelevantColumnNum]=size(CMatSparseFull);
RowNumVec=kron((1:1:rowNum)',ones(maxRelevantColumnNum,1));
CMatEntryReform=mat2cell(StressEntryList,ones(rowNum,1),maxRelevantColumnNum);
CMatEntryReform=cellfun(@(x)transpose(x),CMatEntryReform,'UniformOutput',0);
ColumnNumVec=cell2mat(CMatEntryReform);

CMatReform=mat2cell(CMatSparseFull,ones(rowNum,1),maxRelevantColumnNum);
CMatReform=cellfun(@(x) transpose(x),CMatReform,'UniformOutput',0);
ValueVec=cell2mat(CMatReform);

CMatSparse=sparse(RowNumVec,ColumnNumVec,ValueVec);
CMatSparse=CMatSparse(:,1:end-1);
clear NumStrComp LineNumVec ColumnNumVec ValueVec 

end
