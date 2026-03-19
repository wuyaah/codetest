function [EleAMat,EleBMat]=VarConvertCPE(EleCMat,EleYield)
%This function converts the primitive stress variable to corresponding variable easier for calculating equivalent Mises stress

EleCMat=cell2mat(EleCMat);
EleCMat=mat2cell(EleCMat,size(EleCMat,1),3*ones(size(EleCMat,2)/3,1));
%Inserting a zero column to the equilibrium matrix. 
EleCMat=cellfun(@(C)[C(:,1:2),zeros(size(C,1),1),C(:,3)],EleCMat,'UniformOutput',0);

T_MatInv=[1 -1 0 0; 0 1 -1 0; 1 0 1 0; 0 0 0 sqrt(6)];
L_MatTrans=transpose([sqrt(2) 0 0; 1/sqrt(2) sqrt(3)/sqrt(2) 0; 0 0 1]);

FullCTMat=cellfun(@(C)C/T_MatInv,EleCMat,'UniformOutput',0);

EleDMat=cellfun(@(C)[C(:,1:2),C(:,4)],FullCTMat,'UniformOutput',0);
EleAMat=cellfun(@(D)D/L_MatTrans*EleYield*sqrt(2),EleDMat,'UniformOutput',0);
EleAMat=cell2mat(EleAMat);
EleAMat=mat2cell(EleAMat,2*ones(size(EleAMat,1)/2,1),size(EleAMat,2));

EleBMat=cellfun(@(C)C(:,3),FullCTMat,'UniformOutput',0);
EleBMat=cell2mat(EleBMat);
EleBMat=mat2cell(EleBMat,2*ones(size(EleBMat,1)/2,1),size(EleBMat,2));

end
