function  [info,totalElemList,elemNodeLabelSet,totalNodeList,frozenElem]=OriModelDataDamageFreeNew(dataName,alternativeLoadCases)
%%% Resort the Original model data to Structured Data

loadData=load(dataName,'-regexp','^(Ele)...');
totalElemList=cell2mat(loadData.EleLabelSet);
seqElem=num2cell(transpose(1:numel(totalElemList))); % all elements are sorted from 1 to N

loadData=loadDateProcessing(loadData,seqElem);
elemNodeLabelSet=cellfun(@(x)cell2mat(x),loadData.EleNodeLabelSet,'UniformOutput',0);
totalNodeList=cell2mat(elemNodeLabelSet);
totalNodeList=unique(totalNodeList);

elemNodeLabelSet=cellfun(@(x)transpose(x),elemNodeLabelSet,'UniformOutput',0);
loadData.EleNodeLabelSet=elemNodeLabelSet;
info=cellfun(@(x)getInfoStructGen(loadData,x),seqElem,'UniformOutput',0);

for loadCase=alternativeLoadCases
    if exist([cell2mat(loadCase) '.mat'],'file')
        loadData=load([cell2mat(loadCase) '.mat'],'EleNodeDispSet');
        elemNodeDispSet=cellfun(@(x)squeeze(loadData.EleNodeDispSet(x,:,:)),seqElem,'UniformOutput',0);
        info=cellfun(@(x,disp)getAlternativeNodeDisp(x,disp),info,elemNodeDispSet,'UniformOutput',0);
    end
end

elemNodeLabelSet=getModifiedElemNodeLable(elemNodeLabelSet); %不足最大单元节点数的补足"0"
elemNodeLabelSet=cell2mat(elemNodeLabelSet);
frozenElem=getFrozenElem(dataName,totalElemList);

end

