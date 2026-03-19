function deltaLambda1_new=subNormUpdate(deltaLambda1,stress,elemElaStress)

gradFlowDir=cellfun(@(x)cellfun(@(y)subCalEffStressDiff(y),x,'UniformOutput',0),stress,'UniformOutput',0);

strainEngy=cellfun(@(g,e)cellfun(@(x,y)transpose(x)*y,g,e),gradFlowDir,elemElaStress,'UniformOutput',0);
strainEngy=cellfun(@(x,y)sum(sum(x.*y)),strainEngy,deltaLambda1);

weight=1/sum(strainEngy);

deltaLambda1_new=cellfun(@(x)weight.*x,deltaLambda1,'UniformOutput',0);

% %%验证塑形功归一化结果为1
% strainEngy=cellfun(@(g,e)cellfun(@(x,y)transpose(x)*y,g,e),gradFlowDir,elemElaStress,'UniformOutput',0);
% strainEngy=cellfun(@(x,y)sum(sum(x.*y)),strainEngy,deltaLambda1_new);
% sum(strainEngy);

end

% function deltaLambda1_new=subNormUpdate(deltaLambda1,stress,elemElaStress)
% 
% gradFlowDir=cellfun(@(x)cellfun(@(y)subCalEffStressDiff(y),x,'UniformOutput',0),stress,'UniformOutput',0);
% 
% strainEngy=cellfun(@(g,e)cellfun(@(x,y)transpose(x)*y,g,e),gradFlowDir,elemElaStress,'UniformOutput',0);
% strainEngy=cellfun(@(x,y)sum(sum(x.*y)),strainEngy,deltaLambda1);
% 
% pos=(strainEngy>0);
% weight=zeros(size(strainEngy));
% weight(pos)=1./strainEngy(pos);
% 
% deltaLambda1_new=cellfun(@(w,x)w.*x,num2cell(weight),deltaLambda1,'UniformOutput',0);
% 
% % %%验证塑形功归一化结果为1
% % strainEngy=cellfun(@(g,e)cellfun(@(x,y)transpose(x)*y,g,e),gradFlowDir,elemElaStress,'UniformOutput',0);
% % strainEngy=cellfun(@(x,y)sum(sum(x.*y)),strainEngy,deltaLambda1_new);
% % sum(strainEngy);
% 
% end

