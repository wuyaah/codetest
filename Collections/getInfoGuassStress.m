function GuassStress=getInfoGuassStress(GuassStress)

GuassStress=cell2mat(GuassStress);

%针对缩减积分
if iscolumn(GuassStress) 
    GuassStress=transpose(GuassStress);
end

numStrComp=size(GuassStress,2);
if numStrComp==6
    GuassStress=GuassStress(:,[1:4,6,5]);
end
GuassStress=mat2cell(GuassStress,size(GuassStress,1),size(GuassStress,2));

end
