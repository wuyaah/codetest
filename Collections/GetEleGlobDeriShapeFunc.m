function [GaussGlobalDiffShapes,GaussVol,EleInfo]=GetEleGlobDeriShapeFunc(EleNodeCoordSets,ElementType)
% Given ElementType as a String, and EleNodeCoordSets Coordinates of Each Nodes in the Format [x1 y1 z1;...; x_n y_n z_n]; 
% Function returns a [3,NK(n)*NGS(m)] matrix, whose element is
% [N1,x(Gauss1) N2,x(Gauss1)... N_n,x(Gauss_1) N1,x(Gauss_2)...N_n,x(Gauss_m);
%  N1,y(Gauss1) N2,y(Gauss1)... N_n,y(Gauss_1) N1,y(Gauss_2)...N_n,y(Gauss_m);
%  N1,z(Gauss1) N2,z(Gauss1)... N_n,z(Gauss_1) N1,z(Gauss_2)...N_n,z(Gauss_m)];
% e.g.: For a C3D8 Element with 8 Gauss Point, Result is a (3x64) matrix
% Note: GaussVol is already associated with the weight factor

ElePara=eval([ElementType '_Para']);
EleInfo=ElePara.EleInfo;
GaussPoint=ElePara.GaussCoord;
ELeGKWeight=ElePara.EleInfo.Weight;
EleDim=ElePara.EleInfo.Dimension;

if iscell(EleNodeCoordSets); EleNodeCoordSets=cell2mat(EleNodeCoordSets); end
if EleDim==2; EleNodeCoordSets=EleNodeCoordSets(:,1:2); end

GaussCell=mat2cell(GaussPoint,ones(size(GaussPoint,1),1),size(GaussPoint,2)); % Local CS of each integral poin
GaussShapeFuncs=eval(['cellfun(@' ElementType ',GaussCell)']); % eval([cellfun(@CPS4R,GaussCell)])
GaussShapeFuncs=struct2cell(GaussShapeFuncs);
GaussLocalDiffShapes=GaussShapeFuncs(3,:); 
GaussJacobian=cellfun(@(x)transpose(x)*EleNodeCoordSets,GaussLocalDiffShapes,'UniformOutput',0);
GaussVolUnweight=cellfun(@(x)det(x),GaussJacobian,'UniformOutput',false);         
GaussVol=ELeGKWeight.*cell2mat(GaussVolUnweight);
GaussInvJacobian=cellfun(@(x)inv(x),GaussJacobian,'UniformOutput',false);
GaussGlobalDiffShapes=cellfun(@(x,y)x*transpose(y),GaussInvJacobian,GaussLocalDiffShapes,'UniformOutput',false);
GaussGlobalDiffShapes=cell2mat(GaussGlobalDiffShapes);

end

