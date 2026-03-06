function Info=SortGlobDisp2Ele(Info,NodeDispVec)

EleNodeDisp=cellfun(@(x)x(Info.posElemDoFInGlobDoF),NodeDispVec,'UniformOutput',0);
Info.nodeDisp=cell2mat(EleNodeDisp);

end
