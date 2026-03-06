function plotTopoCurve(volumeFracAttainedSet, optType, objSet, elaLimitSet)

iter=numel(objSet);
figure;
% 设置当前图窗和坐标轴的背景颜色为白色
set(gcf, 'PaperPositionMode', 'auto');  % 设置为自动调整纸张尺寸
set(gcf,'Color','w');
set(gca,'Color','w');

%%%%%%%%%%%%%%%% yyaxis left
yyaxis left
if optType==1
%     plot(1:iter,objSet,1:iter,elaLimitSet,'LineWidth',1.5);
%     % legend('Shakedown','Elastic','Interpreter','latex');
%     legend('${\it{\alpha}}_{\rm{SD}}$', '${\it{\alpha}}_{\rm{Ela}}$', 'Interpreter', 'latex');

    plot(1:iter, objSet, 1:iter, objSet./elaLimitSet, 'LineWidth', 1.5);
    ylabel('Load multiplier','Interpreter','latex');
    % legend('Shakedown','Shakedown/Elastic','VolumeFrac','Interpreter','latex');
    legend('${\it{\alpha}}_{\it{SD}}$', '${\it{\alpha}}_{\it{SD}}/{\it{\alpha}}_{\it{Ela}}$', 'Interpreter', 'latex');

elseif optType==2
    objSet = objSet/1000;
    plot(1:iter, objSet, 'LineWidth', 1.5);
    ylabel('Strain Energy [J]', 'Interpreter', 'latex');
elseif optType==3
    plot(1:iter, objSet, 'LineWidth', 1.5);
    ylabel('Von Mises [MPa]', 'Interpreter', 'latex');
end

xlabel('Iteration', 'Interpreter', 'latex');
xlim([0,numel(objSet)]);
%%%%%%%%%%%%%%%% yyaxis right
yyaxis right
plot(1:iter, volumeFracAttainedSet, 'LineWidth', 1.5, 'HandleVisibility','off');
ylabel('$\it{\Omega}/\it{\Omega}_{\rm{0}}$','Interpreter','latex');
set(gca, 'LineWidth', 1, 'FontName', 'TeXGyrePagella', 'FontSize', 15);
ylim([0, 1]);
hold on

end

