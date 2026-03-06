function ConstrainedDoFMat=OriNodeConstrained(NodeConstrainedDof,NodeLabelList)

NodeConstrained=cellfun(@cell2mat,NodeConstrainedDof,'UniformOutput',0);

ConstrainedDoFMat=[ismember(NodeLabelList,NodeConstrained{1}), ...
                        ismember(NodeLabelList,NodeConstrained{2}), ...
                             ismember(NodeLabelList,NodeConstrained{3})];


end