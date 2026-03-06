function [path1, path2, GP, weight] = subLoadingPathFunc_new(numVert, fixLoad2)
% subLoadingPathFunc_new:
%   当 numVert = 2 且 fixLoad2 = 1:
%       path2 = 常值
%   当 numVert = 2 且 fixLoad2 = 0:
%       path2 = 镜像 path1
%   当 numVert = 3 或 4:
%       fixLoad2 必须为 0，否则报错
%       path2 = 镜像 path1

if nargin < 2
    fixLoad2 = 0;
end

if ~ismember(numVert, [2, 3, 4])
    error('numVert must be 2, 3 or 4');
end

%% ★ 新增规则检查：fixLoad2 仅允许在 numVert=2 时为 1
if fixLoad2 == 1 && numVert ~= 2
    error('fixLoad2 = 1 is only allowed when numVert = 2.');
end

%% Step 1: 构造基础函数 func1
func1 = @(x) 0;
segmentLength = 1 / numVert;
for i = 1:numVert
    a = (i - 1) * segmentLength;
    b = i * segmentLength;
    h = b - a;
    switch numVert
        case 4
            if i == 1
                fseg = @(x) sinpi(0.5*(x-a)/h);
            elseif i == 2
                fseg = @(x) 1;
            elseif i == 3
                fseg = @(x) sinpi(0.5*(b-x)/h);
            else
                fseg = @(x) 0;
            end
        case 3
            if i <= 2
                fseg = @(x) sinpi(0.5*x/h);
            else
                fseg = @(x) 0;
            end
        case 2
            fseg = @(x) sinpi(x);
    end
    func1 = @(x) func1(x) + (x >= a & x < b).*fseg(x);
end


%% Step 2: 构造 func2（载荷2）
if numVert == 2 && fixLoad2 == 1
    func2 = @(x) 1;              % ★ 载荷2恒定
else
    func2 = @(x) func1(1 - x);   % 原镜像
end


%% Step 3: Gauss-Legendre 映射
[GP, weight] = GaussLegendreXW(numVert + 1, 0, 1);
move = GP(1);

GP_shifted = GP - move;
GP_shifted(end) = 1;

SP = linspace(0, 1, numVert+1);

func1_scaled = @(x) 0;
func2_scaled = @(x) 0;

for i = 1:numVert
    leftG = GP_shifted(i);
    rightG = GP_shifted(i+1);
    lenG = rightG - leftG;

    leftS = SP(i);
    rightS = SP(i+1);
    lenS = rightS - leftS;

    rescale = @(x) leftS + (x-leftG)/lenG * lenS;

    func1_scaled = @(x) func1_scaled(x) + (x>=leftG & x<rightG).*func1(rescale(x));
    func2_scaled = @(x) func2_scaled(x) + (x>=leftG & x<rightG).*func2(rescale(x));
end

func1_scaled = @(x) func1_scaled(x) + (x==1).*func1(1);
func2_scaled = @(x) func2_scaled(x) + (x==1).*func2(1);


%% Step 4: 周期化 + 平移
periodic1 = @(x) func1_scaled(mod(x,1));
periodic2 = @(x) func2_scaled(mod(x,1));

path1 = @(x) periodic1(x - move);
path2 = @(x) periodic2(x - move);

end
