function elaStress = getElemGaussStress(dataPack)
% Extract elastic stress values at Gauss points from dataPack
% Supports both single and multiple load cases
%
% Inputs:
%   dataPack - Cell array containing element data with GaussStress field
%
% Outputs:
%   elaStress - Concatenated elastic stress vectors for all elements and load cases

caseStress = dataPack{1}.GaussStress;

if iscell(caseStress)  % Multiple load cases
    % Extract and reshape stresses from each element and load case
    EleStressData = cellfun(@(x) transpose(x.GaussStress(:)), dataPack, 'UniformOutput', false);
    EleStressData = cellfun(@(x) cellfun(@(y) reshape(transpose(y), [], 1), x, 'UniformOutput', false), EleStressData, 'UniformOutput', false);
    EleStressData = cellfun(@(x) cell2mat(x), EleStressData, 'UniformOutput', false);
else
    % Single load case, reshape stresses accordingly
    EleStressData = cellfun(@(x) x.GaussStress, dataPack, 'UniformOutput', false);
    EleStressData = cellfun(@(x) reshape(transpose(x), [], 1), EleStressData, 'UniformOutput', false);
end

% Combine all element stress data into one matrix
elaStress = cell2mat(EleStressData);

end

%% 
% function OriElaStress=getElemGaussStress(dataPack,AlternativeLoadCases)
% % Load the stress of alternative load cases (e.g. Ela22, Ela33, ...)
% 
% CaseStress=dataPack{1}.GaussStress;
% 
% if iscell(CaseStress) %如果是cell结构，则代表是由matlab代码计算的弹性场
%     if numel(CaseStress)>1
%         % if size(CaseStress{1},2)==6 %numStrComp==6;
%         % dataPack=cellfun(@(x)rankAlternativeCaseGaussStress(x),dataPack,'UniformOutput',0);
%         % end
%         EleStressData=cellfun(@(x)x.GaussStress,dataPack,'UniformOutput',0);
%         EleStressData=cellfun(@(x)cellfun(@(y)reshape(transpose(y),[],1),x,'UniformOutput',0),EleStressData,'UniformOutput',0);
%         EleStressData=cellfun(@(x)cell2mat(x),EleStressData,'UniformOutput',0);
%     else
%         EleStressData=cellfun(@(x)cell2mat(x.GaussStress),dataPack,'UniformOutput',0);
%         EleStressData=cellfun(@(x)reshape(transpose(x),[],1),EleStressData,'UniformOutput',0);
%     end
%     OriElaStress=cell2mat(EleStressData);
% else %对应从abaqus导出的弹性场
%     EleStressData=cellfun(@(x)x.GaussStress,dataPack,'UniformOutput',0);
%     EleStressData=cellfun(@(x)reshape(transpose(x),[],1),EleStressData,'UniformOutput',0);
%     AlternativeStress=GetAlternativeLoadCaseStress(AlternativeLoadCases);
%     %OriEleStress=[cell2mat(EleStressData),10.*AlternativeStress];
%     OriElaStress=[cell2mat(EleStressData),AlternativeStress];
% end
% 
% end
% 
% %%%-----------------------------------------------------------------
% function res=rankAlternativeCaseGaussStress(Info)
% 
% AlternativeCase=cellfun(@(g)g(:,[1:4,6,5]),Info.GaussStress(:,2:end),'UniformOutput',0);
% res.GaussStress=[Info.GaussStress(1),AlternativeCase];
% 
% end
