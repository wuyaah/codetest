function resiStress=subUpdateResiStress_RSDMS(a0Coe,akCoe,bkCoe,cosComp,sinComp)

% cosComp=paramBasis{1}; 
% sinComp=paramBasis{2};

akCoeComp=cellfun(@(x)cellfun(@(y,ak)y.*ak,num2cell(x),akCoe,'UniformOutput',0),cosComp,'UniformOutput',0);
akCoeComp=cellfun(@(x)sum(cat(3,x{:}),3),akCoeComp,'UniformOutput',0);

bkCoeComp=cellfun(@(x)cellfun(@(y,bk)y.*bk,num2cell(x),bkCoe,'UniformOutput',0),sinComp,'UniformOutput',0);
bkCoeComp=cellfun(@(x)sum(cat(3,x{:}),3),bkCoeComp,'UniformOutput',0);

resiStress=cellfun(@(x,y)a0Coe+x+y,akCoeComp,bkCoeComp,'UniformOutput',0);
resiStress=cellfun(@(x)mat2cell(x,ones(size(x,1),1),size(x,2)),resiStress,'UniformOutput',0);
resiStress=cellfun(@(x)cellfun(@transpose,x,'UniformOutput',0),resiStress,'UniformOutput',0);
resiStress=transpose(resiStress);
resiStress=cat(2,resiStress{:});

% info.resiStress=resiStress;

end 
