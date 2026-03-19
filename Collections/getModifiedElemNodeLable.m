function elemNodeLabelModified=getModifiedElemNodeLable(elemNodeLabelSet)

elemNodeLabelModified=elemNodeLabelSet;
numElemNodeLabelCell=cellfun(@(x)size(x,2),elemNodeLabelSet,'UniformOutput',0);
numElemNodeLabel=cell2mat(numElemNodeLabelCell);
maxNumElemNodeLabel=max(numElemNodeLabel);
if mean(numElemNodeLabel)==maxNumElemNodeLabel
    return
else
    diffNumElemNode=cellfun(@(x)maxNumElemNodeLabel-x,numElemNodeLabelCell,'UniformOutput',0);
    PosModified=(cell2mat(diffNumElemNode)>0);
    SelDiffNumElemNode=diffNumElemNode(PosModified);
    SelElemNodeLabel=elemNodeLabelModified(PosModified);
    SelElemNodeLabel=cellfun(@(x,n)[x,zeros(size(x,1),n)],SelElemNodeLabel,SelDiffNumElemNode,'UniformOutput',0);
    elemNodeLabelModified(PosModified)=SelElemNodeLabel;
end

end

