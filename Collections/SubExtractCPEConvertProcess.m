function parameters=SubExtractCPEConvertProcess(dataPack,NodeLabelList,EleNodeLabelSet,NodeConstrained)

global ABMatSparse CMatSparse

numElem=numel(dataPack);
numNode=numel(NodeLabelList);
dimension=dataPack{1}.Dimension;
numStrComp=dataPack{1}.NumStrComp; %Stress=[S11;S22;S33;S12];
EleCMat=cellfun(@(x) x.EleCMat,dataPack,'UniformOutput',0); 
elemNGK=cellfun(@(x)x.NGK,dataPack);

% Specified Parameters
% T_MatInv=[1 -1 0 0; 0 1 -1 0; 1 0 1 0; 0 0 0 sqrt(6)];   
% L_MatTrans=transpose([sqrt(2) 0 0; 1/sqrt(2) sqrt(3)/sqrt(2) 0; 0 0 1]);
% S2UMat=(1/sqrt(2))*L_MatTrans*T_MatInv([1:2,4],:);
%% DataProcessing
elemYieldCell=cellfun(@(x) x.Yield,dataPack,'UniformOutput',0); 
% EleYield=kron(cell2mat(elemYieldCell),ones(NGK,1));
% EleHardLimit=cellfun(@(x) x.HardLimit,dataPack,'UniformOutput',0); 
% EleHardLimit=cell2mat(EleHardLimit);
% EleHardLimit=kron(EleHardLimit,ones(NGK,1));
% EleRelaHard=(EleHardLimit-EleYield)./EleYield;
% ElePoisson=cellfun(@(x) x.Poisson,dataPack,'UniformOutput',0); 
% ElePoisson=cell2mat(ElePoisson);
% ElePoisson=kron(ElePoisson,ones(NGK,1));
[elemA_Mat,elemB_Mat]=cellfun(@(C,Y)VarConvertCPE(C,Y),EleCMat,elemYieldCell,'UniformOutput',0);
dataPack=cellfun(@(Info,A,B)SortABMat2Elem(Info,A,B),dataPack,elemA_Mat,elemB_Mat,'UniformOutput',0);

[AMatSparseFull,VarU1EntryList]=GetConstrainedGlobalMatSparseList(dataPack,'elemA_Mat',dimension,3,EleNodeLabelSet,NodeLabelList,NodeConstrained);
[BMatSparseFull,VarX1EntryList]=GetConstrainedGlobalMatSparseList(dataPack,'elemB_Mat',dimension,1,EleNodeLabelSet,NodeLabelList,NodeConstrained);
[CMatSparseFull,StressEntryList]=GetConstrainedGlobalMatSparseList(dataPack,'elemCMat',dimension,3,EleNodeLabelSet,NodeLabelList,NodeConstrained);
ABMatSparseFull=[AMatSparseFull,BMatSparseFull];

CMatSparseFull=mat2cell(CMatSparseFull,ones(size(CMatSparseFull,1),1),size(CMatSparseFull,2));
VarUX1EntryList=[VarU1EntryList,max(max(VarU1EntryList))+VarX1EntryList];

[rowABSparse,colABSparse]=size(ABMatSparseFull); 
parameters={numElem;numNode;elemNGK;numStrComp;rowABSparse;colABSparse;dimension};
% EleLabelList=num2str(num2str(EleLabelList));
CMatSparseFull=cell2mat(CMatSparseFull);

% Convert AB Matrix to real sparse form
[LineNum,MaxRelevantColumnNum]=size(ABMatSparseFull);
LineNumVec=kron((1:LineNum)',ones(MaxRelevantColumnNum,1));
ABMatEntryReform=mat2cell(VarUX1EntryList,ones(LineNum,1),MaxRelevantColumnNum);
ABMatEntryReform=cellfun(@(x) transpose(x),ABMatEntryReform,'UniformOutput',0);
ColumnNumVec=cell2mat(ABMatEntryReform);

ABMatReform=mat2cell(ABMatSparseFull,ones(LineNum,1),MaxRelevantColumnNum);
ABMatReform=cellfun(@(x) transpose(x),ABMatReform,'UniformOutput',0);
ValueVec=cell2mat(ABMatReform);

ABMatSparse=sparse(LineNumVec,ColumnNumVec,ValueVec);
ABMatSparse=[ABMatSparse(:,1:(numStrComp-1)*sum(elemNGK)),ABMatSparse(:,(numStrComp-1)*sum(elemNGK)+2:end-1)];  %The idle columns are already removed
clear LineNumVec ColumnNumVec ValueVec

% Convert C Matrix to real sparse form
[LineNum,MaxRelevantColumnNum]=size(CMatSparseFull);
LineNumVec=kron((1:1:LineNum)',ones(MaxRelevantColumnNum,1));
CMatEntryReform=mat2cell(StressEntryList,ones(LineNum,1),MaxRelevantColumnNum);
CMatEntryReform=cellfun(@(x) transpose(x),CMatEntryReform,'UniformOutput',0);
ColumnNumVec=cell2mat(CMatEntryReform);

CMatReform=mat2cell(CMatSparseFull,ones(LineNum,1),MaxRelevantColumnNum);
CMatReform=cellfun(@(x) transpose(x),CMatReform,'UniformOutput',0);
ValueVec=cell2mat(CMatReform);

CMatSparse=sparse(LineNumVec,ColumnNumVec,ValueVec);
CMatSparse=CMatSparse(:,1:end-1);
clear NumStrComp LineNumVec ColumnNumVec ValueVec 

end
