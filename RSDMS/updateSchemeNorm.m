function [dataPack,volumeFracAttained,relDensSet,params]=updateSchemeNorm(dataPack,posFrozenEle,sensitivity,params,updateScheme,iter)

minRelDens=params.MinRelDens;

if updateScheme==1
    
    posOptEle=~posFrozenEle;
    data=dataPack(posOptEle);
    sensitivity=sensitivity(posOptEle);

    xVal       = cellfun(@(x) x.topoRelDens, data);
    elemVolSet = cellfun(@(x) sum(x.GaussVol), data); 
    fVal = sum(xVal .* elemVolSet) / sum(elemVolSet) - params.OptVolumeFrac;     % 当前体积分数约束
    dfdx = transpose(elemVolSet) / sum(elemVolSet);                    % 约束的一阶导数
    df0dx = sensitivity;    % 目标函数梯度

    [optRelDensSet,params]=SubMMAFunc(xVal,fVal,dfdx,df0dx,params,iter);
    relDensSet=ones(size(dataPack));
    relDensSet(posOptEle)=optRelDensSet;

elseif updateScheme==2
    volumeFrac=params.VolumeFrac;
    relDensSet=SubOCFunc(dataPack,sensitivity,posFrozenEle,volumeFrac,minRelDens);
end

% OldRelDensSet=cellfun(@(x)x.topoRelDens,DataPack);
% RelDensSet=(RelDensSet+OldRelDensSet)/2;
% RelDensSet=GetFilteredSens(DataPack,RelDensSet,2);
% if Iter>30
%     RelDensSet(RelDensSet>0.5)=Params.MaxRelDens;
%     RelDensSet(RelDensSet<0.5)=Params.MinRelDens;
% end
% RelDensSet(RelDensSet>=1-MinRelDens)=1; % AlternativeUpdateRelDens
% oriRelDensSet=cellfun(@(x)x.topoRelDens,dataPack);

for elem=1:numel(dataPack)
    dataPack{elem}.topoRelDens=relDensSet(elem); 
end

elemVolSet=cellfun(@(x)sum(x.GaussVol),dataPack);
volumeFracAttained=sum(relDensSet.*elemVolSet)/sum(elemVolSet);

end
