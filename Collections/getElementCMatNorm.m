function [ElementCMat,OriBMat,GaussVolMat,EleInfo,DeriShapeFuncGauss]=getElementCMatNorm(EleNodeCoordSet,EleType)

[DeriShapeFunc,GaussVolMat,EleInfo]=GetEleGlobDeriShapeFunc(EleNodeCoordSet,EleType);

EleNK=EleInfo.NK;
EleNGK=EleInfo.NGK;
EleDim=EleInfo.Dimension;
if EleDim==3; EleStrainCom=6; else; EleStrainCom=3; end

% Original B Matrix
OriBMat=DeriShapeFunc;
OriBMat=mat2cell(OriBMat,EleDim,ones(EleNK*EleNGK,1));
OriBMat=cellfun(@(x)GlobDeriShapeFunc2BMatCom(x,EleDim),OriBMat,'UniformOutput',false);  
OriBMat=transpose(mat2cell(cell2mat(OriBMat),EleStrainCom,(EleDim*EleNK)*ones(1,EleNGK)));

GaussVolCell=num2cell(GaussVolMat);
DeriShapeFuncGauss=mat2cell(DeriShapeFunc,EleDim,EleNK*ones(EleNGK,1));
DiscDeriShapeFuncByVol=cellfun(@(x,y)mtimes(x,y),GaussVolCell,DeriShapeFuncGauss,'UniformOutput',0); 
DiscDeriShapeFuncByVol=cell2mat(DiscDeriShapeFuncByVol);
DiscDeriShapeFuncByVol=mat2cell(DiscDeriShapeFuncByVol,EleDim,ones(EleNK*EleNGK,1));
OriBMatGaussByVol=cellfun(@(x)GlobDeriShapeFunc2BMatCom(x,EleDim),DiscDeriShapeFuncByVol,'UniformOutput',0);  
OriBMatGaussByVol=cell2mat(OriBMatGaussByVol);

SortedBMatByVolFull=transpose(mat2cell(OriBMatGaussByVol,EleStrainCom,(EleDim*EleNK)*ones(1,EleNGK)));
SortedBMatByVolFull=cell2mat(SortedBMatByVolFull);

if EleDim==3
    ElementCMat=mat2cell(SortedBMatByVolFull,6*EleNGK,3*ones(1,EleNK));
elseif EleDim==2
    ElementCMat=mat2cell(SortedBMatByVolFull,3*EleNGK,2*ones(1,EleNK));
end

ElementCMat=cellfun(@(x)transpose(x),ElementCMat,'UniformOutput',0);
ElementCMat=transpose(ElementCMat);

end
