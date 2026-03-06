function error=ValidForMatrix(dataPack,elemGaussStrainSet)
% 通过与ABAQUS几何线性分析得到的进行比较验证本程序中所涉及的矩阵
if ~isempty(elemGaussStrainSet)
    % NumElem=numel(DataPack);
    elemBMatCell=cellfun(@(x)x.elemBMat,dataPack,'UniformOutput',0);
    elemDispCell=cellfun(@(x)x.nodeDisp(:,1),dataPack,'UniformOutput',0);
    if size(elemGaussStrainSet{1},2)==1 %缩减积分？？ or 不同单元类型？？
        elemBMatCell=cellfun(@(x)sum(cat(3,x{:}),3)./numel(x),elemBMatCell,'UniformOutput',0);
        EleStrainSetMLab=cellfun(@(BMat,Disp)BMat*Disp,elemBMatCell,elemDispCell,'UniformOutput',0);
    else
        EleStrainSetMLab=cellfun(@(BMat,Disp)cellfun(@(B)B*Disp,BMat,'UniformOutput',0),elemBMatCell,elemDispCell,'UniformOutput',0);
        %EleStrainSetMLab=cellfun(@(x)cellfun(@(y)cell2mat(y),x,'UniformOutput',0),EleStrainSetMLab,'UniformOutput',0);
        EleStrainSetMLab=cellfun(@(x)cell2mat(x),EleStrainSetMLab,'UniformOutput',0);
    end
    if iscell(elemGaussStrainSet{1})
        elemGaussStrainSet=cellfun(@(x)cell2mat(x),elemGaussStrainSet,'UniformOutput',0);
    end
    EleStrainSetAba=cellfun(@(x)reshape(transpose(x),[],1),elemGaussStrainSet,'UniformOutput',0);
    %EleStrainSetMLab=cellfun(@(x)reshape(transpose(x),[],1),EleStrainSetMLab,'UniformOutput',0);
    error=abs(cell2mat(EleStrainSetMLab)-cell2mat(EleStrainSetAba));
    maxErr=max(error);
    fprintf('The error of strain between MATLAB and ABAQUS codes is %.6f\n', maxErr);
else
    fprintf('The elemGaussStrainSet is empty\n');
end

end