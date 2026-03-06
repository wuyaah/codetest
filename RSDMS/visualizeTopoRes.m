function visualizeTopoRes(elemCentroidSet, relDensSet, dimension, colorScheme)

if nargin < 4
    colorScheme = 'gray'; % 默认灰度
end

% 清除当前坐标轴内容
cla;

% 根据 colorScheme 设置 colormap
n = 256;
relDensSetGrey = relDensSet;
switch lower(colorScheme)
    case 'gray'
        cmap = gray(n);
        relDensSetGrey = mat2gray(1-relDensSet,[0,1]);
    case 'red'
        cmap = [linspace(1,0.5,n)' linspace(1,0,n)' linspace(1,0,n)']; % 深红 -> 浅红
    case 'blue'
        cmap = [linspace(1,0,n)' linspace(1,0,n)' linspace(1,0.5,n)']; % 深蓝 -> 浅蓝
    case 'green'
        cmap = [linspace(1,0,n)' linspace(1,0.5,n)' linspace(1,0,n)']; % 深绿 -> 浅绿
    otherwise
        warning('未知 colorScheme，使用灰度');
        cmap = gray(n);
end
colormap(cmap);

% 绘制
if dimension==2
    scatter(elemCentroidSet(:,1), elemCentroidSet(:,2), 50, relDensSetGrey, 'filled', 's');
elseif dimension==3
    visualEle = find(relDensSet>0.5);
    relDensSetGrey = relDensSetGrey(visualEle,:);
    elemCentroidSet = elemCentroidSet(visualEle,:);
    scatter3(elemCentroidSet(:,1), elemCentroidSet(:,2), elemCentroidSet(:,3), 50, relDensSetGrey, 'filled', 's');
end

% 坐标轴属性
axis equal;
axis off;
drawnow;

end


