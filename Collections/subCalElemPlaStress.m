function plaStress=subCalElemPlaStress(stress,elemYield)

effStress=cellfun(@subCalEffStress,stress);
plaCoe=effStress>elemYield;
plaCoe=plaCoe.*(effStress-elemYield)./(effStress+eps);
plaStress=cellfun(@(x,Coe)Coe.*x,stress,num2cell(plaCoe),'UniformOutput',0);

end


% function [plaStress, plaCoe] = subCalElemPlaStress(stress, elemYield, beta)
% % -------------------------------------------------------------------------
% % 输入:
% %   stress    - 单元应力元胞数组
% %   elemYield - 单元屈服强度
% %   beta      - 平滑控制参数 (建议外围拓扑迭代中从 50 逐渐增大至 1e4)
% % 输出:
% %   plaStress - 虚拟塑性应力 (用于RSDMS内循环残余应力迭代)
% %   plaCoe    - 平滑塑性系数 (即伪对偶乘子 \xi_smooth，用于外循环敏感度计算)
% % -------------------------------------------------------------------------
% 
% % 如果未输入beta，赋予一个初始中等平滑度作为默认值
% if nargin < 3
%     beta = 50; 
% end
% 
% % 计算等效应力 (调用你原有的外置函数)
% effStress = cellfun(@subCalEffStress, stress);
% 
% % 计算无量纲自变量 x = \beta * (\bar{\sigma} - \sigma_Y)
% x = beta * (effStress - elemYield);
% 
% % --- 核心：数值稳定的 Softplus 实现 ---
% % softplus(x) = log(1 + exp(x))
% % 为了防止 x 很大时 exp(x) 溢出，采用等价数学形式：max(x,0) + log(1 + exp(-abs(x)))
% softplus_x = max(x, 0) + log(1 + exp(-abs(x)));
% 
% % 计算平滑塑性系数 \xi_{smooth}
% % 注意：分母仍需保留 +eps 防止低应力/零应力区域除零错
% plaCoe = softplus_x ./ (beta * (effStress + eps));
% 
% % 计算虚拟塑性应力分量 (利用 cellfun 将标量乘子作用于张量)
% plaStress = cellfun(@(s, Coe) Coe .* s, stress, num2cell(plaCoe), 'UniformOutput', false);
% 
% end

