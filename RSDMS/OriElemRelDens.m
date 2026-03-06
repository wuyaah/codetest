function [OriRelDensSet,Param]=OriElemRelDens(dataPack,PosFrozenEle,OriRelDens,VolumeFrac,MinRelDens)

DesignalElems=~PosFrozenEle;
ElemVolSet=cellfun(@(x)sum(x.GaussVol),dataPack);

DesiredVolume=sum(ElemVolSet)*VolumeFrac;
FrozenEleVolume=sum(ElemVolSet(PosFrozenEle));
DesignalEleVolume=sum(ElemVolSet(DesignalElems));

OPTVolumeFrac=(DesiredVolume-FrozenEleVolume)/DesignalEleVolume;

OriRelDensSet=ones(size(dataPack));

if sum(PosFrozenEle)>0 && OriRelDens==VolumeFrac %% 如果存在冻结区域，OriRelDens暂不能自定义。
    OriRelDensSet(DesignalElems)=OPTVolumeFrac;
else
    OriRelDensSet=OriRelDens.*OriRelDensSet;
end

NumOPTEle=numel(find(DesignalElems));
Param=SetmmaParams(NumOPTEle,OriRelDensSet,MinRelDens);
Param.TolEleVolume=sum(ElemVolSet);
Param.OptEleVolume=DesignalEleVolume;
Param.FroEleVolume=FrozenEleVolume;
Param.OptVolumeFrac=OPTVolumeFrac;
Param.VolumeFrac=VolumeFrac;

end