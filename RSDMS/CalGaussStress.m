function info=CalGaussStress(info,elaMatrix)

NGK=info.NGK;
elemBMat=info.elemBMat;
elaMatrix=info.IntE.*elaMatrix;
numStrComp=info.numStrComp;
nodeDisp=info.nodeDisp;

GaussStrain=cellfun(@(x)x*nodeDisp,elemBMat,'UniformOutput',0);
GaussStress=cellfun(@(x)elaMatrix*x,GaussStrain,'UniformOutput',0);
GaussStress=mat2cell(cell2mat(GaussStress),numStrComp*NGK,ones(1,size(nodeDisp,2)));
GaussStress=cellfun(@(x)reshape(x,numStrComp,[]),GaussStress,'UniformOutput',0);
GaussStress=cellfun(@(x)transpose(x),GaussStress,'UniformOutput',0);

info.GaussStrain=transpose(cell2mat(GaussStrain'));
info.GaussStress=GaussStress(:);

end
