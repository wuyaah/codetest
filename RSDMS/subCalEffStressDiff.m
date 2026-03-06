function diffVal=subCalEffStressDiff(stressVec)
% SUBCALEFFSTRESSDIFF 计算 von Mises 应力对各分量的导数
%
% 输入:
%   stressVec - 应力向量或矩阵
%
% 输出:
%   diffVal   - 导数矩阵/向量，维度与输入 stressVec 完全一致
%
% 逻辑:
%   1. 自动识别 3/4/6 分量 (CPS/CPE/3D)。
%   2. 自动识别单列向量并处理。
%   3. 处理 0 应力奇异点 (返回 0)。

% --- 1. 维度与方向判断 ---
[rows, cols] = size(stressVec);
isColVector = false;

% 线性判断：如果是单列且行数为 3,4,6，标记为列向量并转置处理
if cols == 1 && (rows == 3 || rows == 4 || rows == 6)
    x = stressVec';   % 转置为 [1 x M] 计算
    numStrComp = rows;
    isColVector = true;
else
    x = stressVec;    % 保持 [N x M]
    numStrComp = cols;
end

% --- 2. 提取分量 & 计算 VM 应力 (作为分母) ---
sx = x(:, 1);
sy = x(:, 2);

% 预先计算 von Mises 应力 (为了后面做分母)
switch numStrComp
    case 6 % 3D
        sz  = x(:, 3); txy = x(:, 4); tyz = x(:, 5); tzx = x(:, 6);
        vm = sqrt(0.5 * ((sx - sy).^2 + (sy - sz).^2 + (sz - sx).^2 + 6*(txy.^2 + tyz.^2 + tzx.^2)));
    case 4 % CPE
        sz  = x(:, 3); txy = x(:, 4);
        vm = sqrt(0.5 * ((sx - sy).^2 + (sy - sz).^2 + (sz - sx).^2 + 6*txy.^2));
    case 3 % CPS
        txy = x(:, 3);
        vm = sqrt(sx.^2 - sx.*sy + sy.^2 + 3*txy.^2);
    otherwise
        error('输入维度错误 (需为3, 4, 6)');
end

% --- 3. 准备公共系数 (1 / (2*vm)) 并处理奇异点 ---
% 技巧：为了避免除以0，先计算系数，然后将 vm=0 的位置强制置0
% 公式为: d(vm)/d(comp) = (d(vm^2)/d(comp)) * (1 / (2*vm))

coef = zeros(size(vm));
mask_nonzero = vm > 1e-12; % 设置一个极小的容差

% 仅对非零应力点计算系数，零应力点保持为 0
coef(mask_nonzero) = 0.5 ./ vm(mask_nonzero);

% --- 4. 计算导数 ---
% 初始化输出矩阵
d_out = zeros(size(x)); 

switch numStrComp
    case 6 % 3D: [Sx, Sy, Sz, Txy, Tyz, Tzx]
        % d(vm^2)/dSx = 2Sx - Sy - Sz
        d_out(:,1) = coef .* (2*sx - sy - sz); 
        d_out(:,2) = coef .* (2*sy - sx - sz);
        d_out(:,3) = coef .* (2*sz - sx - sy);
        % d(vm^2)/dTxy = 6*Txy
        d_out(:,4) = coef .* (6*txy); % = 3*Txy/vm
        d_out(:,5) = coef .* (6*tyz);
        d_out(:,6) = coef .* (6*tzx);

    case 4 % CPE: [Sx, Sy, Sz, Txy] (包含 Sz 的导数)
        d_out(:,1) = coef .* (2*sx - sy - sz);
        d_out(:,2) = coef .* (2*sy - sx - sz);
        d_out(:,3) = coef .* (2*sz - sx - sy);
        d_out(:,4) = coef .* (6*txy);

    case 3 % CPS: [Sx, Sy, Txy] (Sz=0, 不参与求导)
        % d(vm^2)/dSx = 2Sx - Sy
        d_out(:,1) = coef .* (2*sx - sy);
        d_out(:,2) = coef .* (2*sy - sx);
        d_out(:,3) = coef .* (6*txy);
end

% --- 5. 恢复原始维度 ---
if isColVector
    diffVal = d_out'; % 如果输入是列向量，输出转回列向量
else
    diffVal = d_out;
end

end
