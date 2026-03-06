function info=getElemsINGivenRadius(info,EleCentroidSet,Radius)
% DataPack=GetElemsWithinGivenRadius(DataPack,Radius)
% Date: 2022.11.6
% This function find elements whose centroid lies within the given filter radius for each element belonging to DataPack
% Note that for convinience, instead of labels of adjacent elements, their entries in EleLabelList are recorder

SelEleCentroidMat=repmat(info.Centroid,size(EleCentroidSet,1),1);
Dist2SelEle=vecnorm(EleCentroidSet-SelEleCentroidMat,2,2);
ValidCriterion=(Dist2SelEle<Radius)&(Dist2SelEle~=0);
DistEleWithinRadius=Dist2SelEle(ValidCriterion);
EntriesAdjacentEle=find(ValidCriterion);

info.DistEleWithinRadius=DistEleWithinRadius;
info.EntriesAdjacentEle=EntriesAdjacentEle;

end

