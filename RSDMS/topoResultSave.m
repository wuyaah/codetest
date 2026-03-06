%%%%%%%%%%%%%%%%%%%% ——保存计算结果—— %%%%%%%%%%%%%%%%%%%%
% OutputPath=('D:\OutputTopo');
% mainPath=('C:\2023.CodeforWork\GurobiNestedTopo'); cd(mainPath);
PathSavedFile=fullfile(outputPath,[dataName,'_',num2str(optType),num2str(iskinemHardening)]);

%%
if optType==1     % shakedown based optimization
    try
        save([PathSavedFile '_SDRes.mat'], ...
            'dataPack','matInfo','optType','result', ...
            'relDensHistory','volumeFracAttainedSet','objSet','elaLimitSet', ...
            'nodeLabelList','nodeConstrained','elemCentroidSet','alternativeLoadCases','-v7.3');
        disp([dataName '_SDRes was saved']);
    catch
        disp([dataName '_SDRes was not saved']);
    end
elseif optType==2 % compliance based optimization
    try
        save([PathSavedFile '_SERes.mat'], ...
            'dataPack','matInfo','optType',...
            'relDensHistory','volumeFracAttainedSet','objSet', ...
            'nodeLabelList','nodeConstrained','elemCentroidSet','-v7.3');
        disp([dataName '_SERes was saved']);
    catch
        disp([dataName '_SERes was not saved']);
    end
elseif optType==3 % stress based optimization
    try
        save([PathSavedFile '_VMRes.mat'], ...
            'dataPack','matInfo','optType',...
            'relDensHistory','volumeFracAttainedSet','objSet', ...
            'nodeLabelList','nodeConstrained','elemCentroidSet','-v7.3');
        disp([dataName '_VMRes was saved']);
    catch
        disp([dataName '_VMRes was not saved']);
    end
end
%%%%%%%%%%%%%%%%%%%% ——保存计算结果—— %%%%%%%%%%%%%%%%%%%%
