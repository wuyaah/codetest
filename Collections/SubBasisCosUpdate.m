function bkCoe=SubBasisCosUpdate(resiStressRate,GPoint,weight,numTerms)

% NGK=numel(resiStressRate{1});
% numStrComp=numel(resiStressRate{1}{1});

bkCoe=cellfun(@(k) ...
    cellfun(@(r,t,w)(1)/(k*pi)*(w*cospi(2*k*t).*r),resiStressRate,GPoint,weight,'UniformOutput',0), ...
    numTerms,'UniformOutput',0);
bkCoe=cellfun(@(x)sum(cat(3,x{:}),3),bkCoe,'UniformOutput',0); %时间积分

end 

