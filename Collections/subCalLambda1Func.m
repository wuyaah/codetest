function [deltaLambda1,plaStress]=subCalLambda1Func(info,deltaLambda1,stress,plaStress)
%% matInfo
E  = info.IntE;       % Young's modulus
nu = info.Poisson;    % Poisson's ratio
% Lame constant value
mu=E/(2*(1+nu));      % shear modulus
lambda=(E*nu)/((1+nu)*(1-2*nu));
elaMatrix_Prime=2*mu*eye(3)+lambda*ones(3); % 主应力空间下的等效弹性矩阵 Eij

%%
Yield=info.Yield;
% Yield=(1-0.01)*Yield; %软化屈服面: 定义有效集合
% plaStress=subCalElemPlaStressNew(stress,Yield);

pos=cellfun(@(x)norm(x)>0,plaStress,'UniformOutput',1);
[stress,dir]=cellfun(@(x)getEigenStress(x),stress(pos),'UniformOutput',0);% dir{1}*diag(stress{1})*dir{1}'
gradFlowDir=cellfun(@(x)subFlowDirect(x),stress,'UniformOutput',0);
deltaStress=cellfun(@(x,d)subCalDeltaStress(x,Yield),stress,'UniformOutput',0);
eigenStress=cellfun(@(x)elaMatrix_Prime*x,gradFlowDir,'UniformOutput',0);

deltaLambda1_new=zeros(size(deltaLambda1));
deltaLambda1_new(pos)=cellfun(@(d,x)mean(d./x),deltaStress,eigenStress,'UniformOutput',1);
% deltaLambda1_new(pos)=cellfun(@(d,x)d./x,deltaStress,eigenStress,'UniformOutput',0);

deltaLambda1=deltaLambda1+deltaLambda1_new; 

plaStress_pos_tensor=cellfun(@(x,d)d*diag(x)*d',deltaStress,dir,'UniformOutput',0);
switch info.numStrComp
    case 6
        plaStress_pos=cellfun(@(x)[x(1,1);x(2,2);x(3,3);x(1,2);x(1,3);x(2,3)],plaStress_pos_tensor,'UniformOutput',0);
    case 4
        plaStress_pos=cellfun(@(x)[x(1,1);x(2,2);x(3,3);x(1,2)],plaStress_pos_tensor,'UniformOutput',0);
    case 3
        plaStress_pos=cellfun(@(x)[x(1,1);x(2,2);x(1,2)],plaStress_pos_tensor,'UniformOutput',0);
end
plaStress(pos)=plaStress_pos;

end


function plaStress=subCalElemPlaStressNew(stress,elemYield)

effStress = cellfun(@subCalEffStress, stress);
plaCoe = effStress > elemYield;
plaCoe = plaCoe.*(effStress-elemYield)./(effStress+eps);
plaStress = cellfun(@(x, Coe)Coe.*x, stress, num2cell(plaCoe), 'UniformOutput', false);

end


% function [deltaLambda1, plaStress] = subCalLambda1Func(info, deltaLambda1, stress, plaStress, plaCoe)
%     %% matInfo 初始化
%     E  = info.IntE;       
%     nu = info.Poisson;    
%     mu = E / (2 * (1 + nu));      
%     lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
%     elaMatrix_Prime = 2 * mu * eye(3) + lambda * ones(3); 
% 
%     Yield = info.Yield;
% 
%     % ================= 修复核心区 =================
%     % 1. 重新计算真实的等效应力
%     effStress = cellfun(@subCalEffStress, stress);
% 
%     % 2. 弃用 norm(plaStress)，改用严格的物理截断准则
%     % 引入 0.9999 系数是为了容忍浮点数舍入误差，防止临界单元漏判
%     pos = effStress >= (Yield * 0.9999); 
%     % ==============================================
% 
%     deltaLambda1_new = zeros(size(deltaLambda1));
% 
%     % --- 分支 A：处理真正发生塑性屈服的激活单元 (精确矩阵计算) ---
%     if any(pos)
%         [stress_active, dir] = cellfun(@(x)getEigenStress(x), stress(pos), 'UniformOutput', 0);
%         gradFlowDir = cellfun(@(x)subFlowDirect(x), stress_active, 'UniformOutput', 0);
%         deltaStress = cellfun(@(x,d)subCalDeltaStress(x, Yield(pos)), stress(pos), 'UniformOutput', 0);
%         eigenStress = cellfun(@(x)elaMatrix_Prime*x, gradFlowDir, 'UniformOutput', 0);
% 
%         % 这里不再会有 0/0，因为进入此分支的必然是 deltaStress > 0 的单元
%         deltaLambda1_new(pos) = cellfun(@(d,x) mean(d./x), deltaStress, eigenStress, 'UniformOutput', 1);
% 
%         plaStress_pos_tensor = cellfun(@(x,d) d*diag(x)*d', deltaStress, dir, 'UniformOutput', 0);
%         switch info.numStrComp
%             case 6
%                 plaStress_pos = cellfun(@(x)[x(1,1);x(2,2);x(3,3);x(1,2);x(1,3);x(2,3)], plaStress_pos_tensor, 'UniformOutput', 0);
%             case 4
%                 plaStress_pos = cellfun(@(x)[x(1,1);x(2,2);x(3,3);x(1,2)], plaStress_pos_tensor, 'UniformOutput', 0);
%             case 3
%                 plaStress_pos = cellfun(@(x)[x(1,1);x(2,2);x(1,2)], plaStress_pos_tensor, 'UniformOutput', 0);
%         end
%         plaStress(pos) = plaStress_pos;
%     end
% 
%     % --- 分支 B：处理弹性安全区的未屈服单元 (注入连续平滑伪梯度) ---
%     not_pos = ~pos;
%     if any(not_pos)
%         % 直接将外层 Softplus 算出的平滑系数赋给弹性单元
%         deltaLambda1_new(not_pos) = plaCoe(not_pos);
%     end
% 
%     % 累加迭代量
%     deltaLambda1 = deltaLambda1 + deltaLambda1_new; 
% end

