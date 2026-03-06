function plotShakedownDomain(resultSet,Yield,appliedLoadSet)

if nargin<=2
    appliedLoadSet=[1,0.1];
end
numLoadCase=size(appliedLoadSet,2);
emptyIndex=cellfun(@isempty,resultSet);
iter=sum(~emptyIndex);

resultSet=resultSet(1:iter);
degree=cellfun(@(x)x.degree,resultSet);
% Degree=deg2rad(Degree);
% AlterDegree=cellfun(@(x)x.AlterDegree,resultSet);
weightSet=[cos(degree),sin(degree),ones(numel(degree),1)];
% WeightSet=[sin(AlterDegree).*cos(degreeSet) sin(AlterDegree).*sin(degreeSet) cos(AlterDegree)]; 

loadFactorSet=cellfun(@(x)x.alpha,resultSet,'UniformOutput',0);
appliedLoadCell=cellfun(@(x)x.*appliedLoadSet,loadFactorSet,'UniformOutput',0);
appliedLoadSet=cell2mat(appliedLoadCell);
appliedLoadSet=[appliedLoadSet,zeros(size(appliedLoadSet,1),3-size(appliedLoadSet,2))];
appliedLoadSet=weightSet.*appliedLoadSet;
limitLoadSet=appliedLoadSet./Yield;

%%
% 绘图，根据加载工况维度决定使用 2D 还是 3D
if numLoadCase == 3
    plot3(limitLoadSet(:,1), limitLoadSet(:,2), limitLoadSet(:,3), '--rs', 'LineWidth', 2);
elseif numLoadCase == 2
    plot(limitLoadSet(:,1), limitLoadSet(:,2), '--rs', 'LineWidth', 2);
end

% 标签和坐标轴设置
xlabel('$\it{Px}/\it{\sigma}_Y$', 'Interpreter', 'latex');
ylabel('$\it{Py}/\it{\sigma}_Y$', 'Interpreter', 'latex');

% 图像显示设置
axis equal
grid on

% 坐标轴范围
xlim([-1,1]);
ylim([-1,1]);

% 设置坐标轴刻度
xticks(-1:0.1:1);
yticks(-1:0.1:1);

drawnow

end
