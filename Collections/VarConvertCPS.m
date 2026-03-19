function  EleABMat=VarConvertCPS(EleCMat,EleYield)
%This function converts the primitive stress variable to corresponding variable easier for calculating equivalent Mises stress
InvL_MatTrans=[ 1, 1/sqrt(3),   0;
                0, 2/sqrt(3),   0; 
                0,    0,    1/sqrt(3)];

EleCMat=cell2mat(EleCMat);
EleCMat=mat2cell(EleCMat,size(EleCMat,1),3*ones(size(EleCMat,2)/3,1));
EleABMat=cellfun(@(C)C*InvL_MatTrans*EleYield,EleCMat,'UniformOutput',0);
EleABMat=cell2mat(EleABMat);
EleABMat=mat2cell(EleABMat,2*ones(size(EleABMat,1)/2,1),size(EleABMat,2));

end

