function info=getAlternativeNodeDisp(info,NodeDisp)

AlternativeNodeDisp=reshape(transpose(cell2mat(NodeDisp)),[],1);
info.nodeDisp=[info.nodeDisp,AlternativeNodeDisp];
    
end