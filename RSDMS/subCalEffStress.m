function effStress = subCalEffStress(stressVec)
% SUBCALEFFSTRESS 计算 von Mises 等效应力:
%   1. 支持 [N x M] 矩阵批量计算 (N个点，M个分量)
%   2. 支持 [M x 1] 单列向量输入 (自动转置为 1 x M 计算)
%
% 输入维度 M 说明:
%   3 -> 平面应力 (CPS): [Sx; Sy; Txy]
%   4 -> 平面应变 (CPE): [Sx; Sy; Sz; Txy]
%   6 -> 三维情况 (3D) : [Sx; Sy; Sz; Txy; Tyz; Tzx]

% --- 1. 行列预处理 ---
[rows, cols] = size(stressVec);

% 如果是“单列”且行数符合应力分量特征 (3, 4, 6)，则判定为单个向量
if cols == 1 && (rows == 3 || rows == 4 || rows == 6)
    x = stressVec';   % 强制转置为行向量 [1 x M]
    numStrComp = rows;
else
    x = stressVec;    % 保持原样 [N x M]
    numStrComp = cols;
end

% --- 2. 提取分量 (使用点运算 . 兼容矩阵和单点) ---
sx = x(:, 1);
sy = x(:, 2);

switch numStrComp
    case 6 % 3D: [Sx, Sy, Sz, Txy, Tyz, Tzx]
        sz  = x(:, 3);
        txy = x(:, 4);
        tyz = x(:, 5);
        tzx = x(:, 6);
        % 3D 完整公式
        effStress = sqrt(0.5 * ((sx - sy).^2 + (sy - sz).^2 + (sz - sx).^2 + 6 * (txy.^2 + tyz.^2 + tzx.^2)));
    case 4 % CPE: [Sx, Sy, Sz, Txy]
        sz  = x(:, 3);
        txy = x(:, 4);
        % 3D 公式退化 (无 Tyz, Tzx)
        effStress = sqrt(0.5 * ((sx - sy).^2 + (sy - sz).^2 + (sz - sx).^2 + 6 * txy.^2));
    case 3 % CPS: [Sx, Sy, Txy]
        txy = x(:, 3);
        % 平面应力简化公式 (Sz=0)
        effStress = sqrt(sx.^2 - sx.*sy + sy.^2 + 3 * txy.^2);
    otherwise
        error('输入错误: 应力分量必须为 3, 4 或 6 (列或行)。当前识别为 %d 列。', numStrComp);
end

end
