function akCoe=SubBasisSinUpdate(resiStressRate,GPoint,weight,numTerms)

% NGK=numel(resiStressRate{1});
% numStrComp=numel(resiStressRate{1}{1});

akCoe=cellfun(@(k) ...
    cellfun(@(r,t,w)(-1)/(k*pi)*(w*sinpi(2*k*t).*r),resiStressRate,GPoint,weight,'UniformOutput',0), ...
    numTerms,'UniformOutput',0);
akCoe=cellfun(@(x)sum(cat(3,x{:}),3),akCoe,'UniformOutput',0); %时间积分

end 