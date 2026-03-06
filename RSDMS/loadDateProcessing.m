function loadData=loadDateProcessing(loadData,seqElem)

if numel(loadData.EleNodeLabelSet{1})==1 % the number of element types used in the model equal to 1
    elemNodeLabelSet=cellfun(@(x)reshape(loadData.EleNodeLabelSet(x,:),[],1),seqElem,'UniformOutput',0);
    elemNodeCoordSet=cellfun(@(x)squeeze(loadData.EleNodeCoordSet(x,:,:)),seqElem,'UniformOutput',0);
    elemNodeDispSet=cellfun(@(x)squeeze(loadData.EleNodeDispSet(x,:,:)),seqElem,'UniformOutput',0);
    elemGaussStressSet=cellfun(@(x)squeeze(loadData.EleGaussStressSet(x,:,:)),seqElem,'UniformOutput',0);
    try
        elemGaussStrainSet=cellfun(@(x)squeeze(loadData.EleGaussStrainSet(x,:,:)),seqElem,'UniformOutput',0);
    catch
        elemGaussStrainSet=[];
    end
else
    return
end

%elemNodeLabelSet=cellfun(@(x)cell2mat(x),elemNodeLabelSet,'UniformOutput',0);
loadData.EleGaussStressSet=elemGaussStressSet;
loadData.EleGaussStrainSet=elemGaussStrainSet;
loadData.EleNodeLabelSet=elemNodeLabelSet;
loadData.EleNodeCoordSet=elemNodeCoordSet;
loadData.EleNodeDispSet=elemNodeDispSet;

end
