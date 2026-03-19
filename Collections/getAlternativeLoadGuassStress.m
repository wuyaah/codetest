function info=getAlternativeLoadGuassStress(info,seqElemLable,alternativeLoadCases)

for LoadCase=alternativeLoadCases
    if exist([cell2mat(LoadCase) '.mat'],'file')
        loadData=load([cell2mat(LoadCase) '.mat'],'EleGaussStressSet');
        elemGaussStressSet=cellfun(@(x)squeeze(loadData.EleGaussStressSet(x,:,:)),seqElemLable,'UniformOutput',0);
        caseStress=cellfun(@(x)getInfoGuassStress(elemGaussStressSet{x}),seqElemLable,'UniformOutput',0);
        caseStress=cellfun(@(x)cell2mat(x),caseStress,'UniformOutput',0);
    for elem=1:numel(info)
        info{elem}.GaussStress=[info{elem}.GaussStress;caseStress(elem)]; 
    end
    end
end

end