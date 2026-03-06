function frozenElems = getFrozenElem(dataName, totalEleList)

if isfile([dataName '.mat']) % 额外检查文件是否存在，增加鲁棒性
    loadData = load(dataName, 'FrozenElems');
else
    loadData = struct();
end

%% 判断变量是否存在 (Check Existence)
if isfield(loadData, 'FrozenElems')
    frozenList = cell2mat(loadData.FrozenElems);
else
    frozenList = [];
end

%% Logical Mask
frozenElems = ismember(totalEleList, frozenList);
    
end