%% 大工论文敏感度表达式+高斯权重
function sens = sensAnalysisBasedOnRSDMS(info, lambda)

GradY=info.GradMat(2);
Weight=repmat(info.GaussVol(:),1,size(lambda, 2));% 获取高斯权重 [NGK x numVert]
sens=-GradY*sum(sum(lambda.*Weight));

end

%% 大工论文敏感度表达式
% function sens=sensAnalysisBasedOnRSDMS(info,lambda)
% 
% 
% elemYield=info.Yield;
% GradY=info.GradMat(2);
% sens=-GradY*sum(sum(lambda))*elemYield;
% 
% end

%% 
% function sens=sensAnalysisBasedOnRSDMS(info,lambda)
% 
% elemYield=info.Yield;
% GradY=info.GradMat(2);
% Weight=repmat(info.GaussVol(:),1,size(lambda, 2));% 获取高斯权重 [NGK x numVert]
% sens=-GradY*sum(sum(lambda.*Weight))*elemYield;
% 
% end

%%
% function sens=sensAnalysisBasedOnRSDMS(dataPack,result)
% 
% numElem=num2cell(1:numel(dataPack))';
% plaCoeSet=result.plaCoeSet;
% plaCoeSet=cellfun(@(x)cellfun(@(y)cell2mat(y'),x,'UniformOutput',0),plaCoeSet,'UniformOutput',0);
% 
% plaCoeMat=cellfun(@(x)cellfun(@(y)y{x},plaCoeSet,'UniformOutput',0),numElem,'UniformOutput',0);
% plaCoeMat=cellfun(@(x)sum(cat(3,x{:}),3),plaCoeMat,'UniformOutput',0);
% plaCoeMat=cell2mat(plaCoeMat);
% 
% NGK=result.params.NGK;
% numElem=result.params.numElem;
% elemPlaCoe=mat2cell(plaCoeMat,NGK*ones(numElem,1),size(plaCoeMat,2));
% % elemPlaCoe=cellfun(@(x)mat2cell(x,size(x,1),ones(1,size(x,2))),elemPlaCoe,'UniformOutput',0);
% 
% sens=cellfun(@(x,p)subCalSensFunc(x,p),dataPack,elemPlaCoe,'UniformOutput',1);
% 
% % sens2=-log10(-sens);
% 
% 
% end
