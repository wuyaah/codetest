function Sensitivity=GetFilteredSens(DataPack,OriSens,Radius)

GaussVol=cellfun(@(x)sum(x.GaussVol),DataPack);
Sensitivity=cellfun(@(x)SubFilteredSens(x,OriSens,GaussVol,Radius),DataPack); 

end

%% SubFunc
function Sens=SubFilteredSens(info,OriSens,GaussVol,Radius)

AdjacentEle=info.EntriesAdjacentEle;
GaussVol=GaussVol(AdjacentEle);
Weight=Radius-info.DistEleWithinRadius; 
Weight(Weight<0)=0;

Sens=sum(Weight.*GaussVol.*OriSens(AdjacentEle));
Sens=Sens/sum(Weight.*GaussVol);

end


