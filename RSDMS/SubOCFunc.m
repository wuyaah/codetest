% function [VolumeFracAttained,RelDensSet]=SubOCFunc(DataPack,Sensitivity,PosFrozenEle,VolumeFrac,MinRelDens)
function RelDensSet=SubOCFunc(DataPack,Sensitivity,PosFrozenEle,VolumeFrac,MinRelDens)

Move=0.2; % Default Move=0.2;
ElemVolSet=cellfun(@(x)sum(x.GaussVol),DataPack);
DesiredVolume=VolumeFrac*sum(ElemVolSet);
OriDens=cellfun(@(x)x.topoRelDens,DataPack);
Sensitivity=abs(Sensitivity); %%% Sensitivity=-Sensitivity;

Lb=0; Ub=1e9;
while Ub-Lb>1e-5
    Threshold=(Ub+Lb)/2;
    RelDensSet=max(MinRelDens,max(OriDens-Move,min(1., ...
               min(OriDens+Move,OriDens.*sqrt(Sensitivity./Threshold)))));
    RelDensSet(PosFrozenEle)=1;
    if sum(RelDensSet.*ElemVolSet)>DesiredVolume
        Lb=Threshold;
    else
        Ub=Threshold;
    end
end

% DataPack=cellfun(@(x,Dens)GetEleRelDens(x,Dens),DataPack,num2cell(RelDensSet),'UniformOutput',0);
% VolumeFracAttained=sum(RelDensSet.*ElemVolSet)/sum(ElemVolSet);

end

