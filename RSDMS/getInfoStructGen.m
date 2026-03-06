function info=getInfoStructGen(loadData,x)

info.elemType=loadData.EleTypeSet{x};
info.elemLabel=loadData.EleLabelSet{x};
info.materialID=loadData.EleMaterialIDSet{x};
info.damageCoe=0; %DamageFree

info.nodeLabel=loadData.EleNodeLabelSet{x};
info.nodeCoord=cell2mat(loadData.EleNodeCoordSet{x});
info.nodeDisp=reshape(transpose(cell2mat(loadData.EleNodeDispSet{x})),[],1);

info.GaussStress=getInfoGuassStress(loadData.EleGaussStressSet{x});

try
    info.numStrComp=size(cell2mat(info.GaussStress),2);
catch
    info.numStrComp=size(info.GaussStress,2);
end

end
