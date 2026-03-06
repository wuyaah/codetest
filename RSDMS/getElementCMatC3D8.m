function [ElementCMat,SortedBMatFull,GaussVolMat,EleInfo,DiscDeriShapeFunc]=getElementCMatC3D8(EleNodeCoordSet,EleType)
% ElementCMat is a cell whose elements are 3-by-6 matrices with the order
% N_m,@Int_n=[N_m,x    0      0     N_m,y      0      N_m,z;
%              0    N_m,y    0     N_m,x    N_m,z       0 ;
%              0      0    N_m,z     0      N_m,y    N_m,x;]@Int_n;
% 
% ElementCMat(1)=[N_1,@Int_1 N_1,@Int_2   ...   N_1,@Int_n;] (3-by-(6*NGK))
% Where:
% N_1,@Int_1= [N_1,x(Int_1)      0          0        N_1,y(Int_1)      0         N_1,z(Int_1);
%                  0      N_1,y(Int_1)     0        N_1,x(Int_1) N_1,z(Int_1)        0;                   
%                  0            0      N_1,z(Int_1)     0        N_1,y(Int_1)   N_1,x(Int_1);]                

% GetElementCMatC3D8 construct a special case C matrix because in ABAQUS, the geometry function of C3D8 is 
% e=sym(du/dx)+f*I*[trace(du_avg/dt)-trace(du-dt)] 
% where: sym(du/dx)->original;f*I*trace(du_avg/dt)->Teil1;f*I*-trace(du-dt)->Teil2
% Rank of ElementCMat is EleNK-by-EleNGK. with each element 3-by-6

[DeriShapeFunc,GaussVol,EleInfo]=GetEleGlobDeriShapeFunc(EleNodeCoordSet,EleType);

% GaussVol=transpose(GaussVol);
% GaussVolMat=cell2mat(GaussVol);
% SumVol=sum(cell2mat(GaussVol));
% GaussVolShare=cellfun(@(x)x./SumVol,GaussVol);

GaussVolMat=transpose(GaussVol);
SumVol=sum(GaussVolMat);
GaussVolShare=GaussVolMat./SumVol;
GaussVolShare=num2cell(GaussVolShare);

EleNK=EleInfo.NK;
EleNGK=EleInfo.NGK;
EleDim=EleInfo.Dimension;
EleNStrainCom=6;  %Number of stress components;

% Original B Matrix
DiscDeriShapeFunc=mat2cell(DeriShapeFunc,EleDim,ones(EleNK*EleNGK,1));

OriBMatGauss=cellfun(@(x)GlobDeriShapeFunc2BMatCom(x,EleDim),DiscDeriShapeFunc,'UniformOutput',false);  
OriBMatGauss=cell2mat(OriBMatGauss);  
SortedOriBMatGauss=transpose(mat2cell(OriBMatGauss,EleNStrainCom,(EleDim*EleNK)*ones(1,EleNGK)));
SortedBMatUni=cellfun(@(x)[x(1,1:3:end);x(2,2:3:end);x(3,3:3:end)],SortedOriBMatGauss,'UniformOutput',0); 
% B Matrix Extended Part1
VolShareMulDiffShape=cellfun(@(x,y)mtimes(x,y),SortedBMatUni,GaussVolShare,'UniformOutput',0); 
VolShareMulDiffShape=cell2mat(VolShareMulDiffShape);
AveVolMulDiffShapeX=sum(VolShareMulDiffShape(1:3:end,:),1);
AveVolMulDiffShapeY=sum(VolShareMulDiffShape(2:3:end,:),1);
AveVolMulDiffShapeZ=sum(VolShareMulDiffShape(3:3:end,:),1);
AveVolMulDiffShape=[AveVolMulDiffShapeX;AveVolMulDiffShapeY;AveVolMulDiffShapeZ];
AveVolMulDiffShape=mat2cell(AveVolMulDiffShape,3,ones(1,EleNK));
AveVolMulDiffShape=transpose(cell2mat(transpose(AveVolMulDiffShape)));
SortedBMatExt1Sec=[AveVolMulDiffShape;AveVolMulDiffShape;AveVolMulDiffShape;zeros(3,24)];
SortedBMatExt1={SortedBMatExt1Sec;SortedBMatExt1Sec;SortedBMatExt1Sec;SortedBMatExt1Sec;...
                SortedBMatExt1Sec;SortedBMatExt1Sec;SortedBMatExt1Sec;SortedBMatExt1Sec};

% B Matrix Extended Part2            
SortedBMatExt2=cellfun(@(x)kron(x,[1 1 1]),SortedBMatUni,'UniformOutput',false); %Pay Attention to the Node Number!!
SortedBMatExt2=cell2mat(SortedBMatExt2);
SortedBMatExt2=mat2cell(SortedBMatExt2,3*ones(1,8),3*ones(1,8));
SortedBMatExt2=cellfun(@transpose,SortedBMatExt2,'UniformOutput',false);
SortedBMatExt2=cell2mat(SortedBMatExt2);
SortedBMatExt2=mat2cell(SortedBMatExt2,3*ones(1,8),24);
SortedBMatExt2=cellfun(@(x)[x;zeros(3,24)],SortedBMatExt2,'UniformOutput',false);

% Sum to Total B Matrix and Calculate Strain Component at Integral Point.
SortedBMatFull=cellfun(@(x,y,z)x+1./3.*(y-z),SortedOriBMatGauss,SortedBMatExt1,SortedBMatExt2,'UniformOutput',false);
% SortedBMatFull=cell2mat(SortedBMatFull);
% SelDispResData=cell2mat(SelDispResData);
% SelDispResData=mat2cell(SelDispResData,ones(8,1),3);
% SelDispResData=cellfun(@transpose,SelDispResData,'UniformOutput',false); 
% SelDispResDataSort=cell2mat(SelDispResData);
% StrainGauss=SortedBMatFull*SelDispResDataSort;
% StrainGaussSorted=mat2cell(StrainGauss,6*ones(1,8),1);
% StrainGaussSorted=transpose(cell2mat(transpose(StrainGaussSorted)));
% 
% Should Be Modified,
% B matrix component should be scaled with Gauss volume
GaussVolCell=num2cell(GaussVolMat);
SortedBMatByVolFull=cellfun(@(x,y)x.*y,SortedBMatFull,GaussVolCell,'UniformOutput',0);
SortedBMatByVolFull=cell2mat(SortedBMatByVolFull);
ElementCMat=mat2cell(SortedBMatByVolFull,6*EleNGK,3*ones(1,EleNK));
ElementCMat=cellfun(@transpose,ElementCMat,'UniformOutput',0);
ElementCMat=transpose(ElementCMat);

end
