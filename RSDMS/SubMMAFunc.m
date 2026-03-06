function [RelDensSet, params] = SubMMAFunc(xVal, fVal, dfdx, df0dx, params, iter)
% SubMMAFunc: 求解 MMA 子问题，更新相对密度变量
%
% 输入参数：
%   dataPack    : 包含每个单元 GaussVol 和 topoRelDens 的结构体数组
%   sensitivity : 目标函数关于相对密度的导数
%   params      : 参数结构体，包含 MMA 所需的变量和常量
%   iter        : 当前迭代次数
%
% 输出参数：
%   RelDensSet  : 更新后的相对密度值
%   params      : 更新后的参数结构体，包括历史变量和对偶变量信息

%% 提取体积、变量和灵敏度信息
% elemVolSet = cellfun(@(x) sum(x.GaussVol), dataPack);            % 每个单元的体积
% xVal       = cellfun(@(x) x.topoRelDens, dataPack);              % 当前设计变量值（列向量）
% volumeFrac = params.OptVolumeFrac;
% df0dx = sensitivity;    % 目标函数梯度
% fVal = sum(xVal .* elemVolSet) / sum(elemVolSet) - volumeFrac;     % 当前体积分数约束
% dfdx = transpose(elemVolSet) / sum(elemVolSet);                    % 约束的一阶导数

numVar     = params.NumVar;
numConstr  = params.NumConstr;

maxRho = params.MaxRelDens * ones(numVar, 1);
minRho = params.MinRelDens * ones(numVar, 1);

xold1 = params.rho1;    % 上一轮变量
xold2 = params.rho2;    % 上上轮变量

%% MMA 参数设置
epsimin = 1e-8;         % 解精度下限
RAA0    = 1e-5;         % 正则化小扰动项
move    = 1.0;          % 最大变量移动量
albefa  = 0.1;          % alpha/beta界限权重
asyinit = 0.5;          % 初始渐近值
asyincr = 1.2;          % 渐近值增长因子
asydecr = 0.7;          % 渐近值缩小因子

NEye = ones(numVar, 1);
MEye = ones(numConstr, 1);

%% 设置设计变量的上下限（渐近边界）
IntRho = maxRho - minRho;

if iter <= 2
    Low = xVal - asyinit * IntRho;
    Upp = xVal + asyinit * IntRho;
else
    [Low, Upp] = getLowUppValue(xVal, params, NEye, asyincr, asydecr);
end

%% 计算 alpha 和 beta
alpha = max(max(Low + albefa * (xVal - Low), xVal - move * IntRho), minRho);
beta  = min(min(Upp - albefa * (Upp - xVal), xVal + move * IntRho), maxRho);

%% 构造 MMA 子问题的目标函数部分
xmamiinv = NEye ./ max(IntRho, 1e-5 * NEye);  % 避免除零
UX2 = (Upp - xVal).^2;
XL2 = (xVal - Low).^2;

uxinv  = NEye ./ (Upp - xVal);
xlinv  = NEye ./ (xVal - Low);

P0 = max(df0dx, 0);
Q0 = max(-df0dx, 0);
PQ0 = 0.001 * (P0 + Q0) + RAA0 * xmamiinv;
P0 = (P0 + PQ0) .* UX2;
Q0 = (Q0 + PQ0) .* XL2;

%% 构造约束部分
P = max(dfdx, 0);
Q = max(-dfdx, 0);
PQ = 0.001 * (P + Q) + RAA0 * MEye * xmamiinv';
P = (P + PQ) * spdiags(UX2, 0, numVar, numVar);
Q = (Q + PQ) * spdiags(XL2, 0, numVar, numVar);

params.b = P * uxinv + Q * xlinv - fVal;  % 右侧常数向量 b

%% 求解子问题
params.Upp = Upp;
params.Low = Low;

[RelDensSet, Ymma, Zmma, Lam, Xsi, Eta, Mu, Zet, S] = getSubSolv(epsimin, alpha, beta, P0, Q0, P, Q, params);

%% 存储求解信息和历史变量
info.Ymma = Ymma;
info.Zmma = Zmma;
info.Lam  = Lam;
info.Xsi  = Xsi;
info.Eta  = Eta;
info.Mu   = Mu;
info.Zet  = Zet;
info.S    = S;

params.info = info;
params.rho2 = xold1;
params.rho1 = xVal;

end
