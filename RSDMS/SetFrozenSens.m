function Sensitivity=SetFrozenSens(OriSensitivity,PosFrozenEle,Params)

Sensitivity=OriSensitivity;
VolumeFrac=Params.VolumeFrac;

NumElem=numel(Sensitivity);
SelPosSens=fix(NumElem*VolumeFrac);
rankingSens=sort(Sensitivity);
SelSens=rankingSens(SelPosSens);

SensFreeze=Sensitivity(PosFrozenEle);
SensFreeze(SensFreeze>SelSens)=SelSens;
Sensitivity(PosFrozenEle)=SensFreeze;

end