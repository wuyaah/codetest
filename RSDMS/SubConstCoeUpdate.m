function a0CoeComp=SubConstCoeUpdate(resiStressRate,weight,akCoe,akCoe_old)

coeBasisCosComp=sum(cat(3,akCoe{:}),3); %级数求和: a(K+1)

% oldCoeBasisCosComp=cellfun(@(x)full(x),akCoe_old,'UniformOutput',0);
oldCoeBasisCosComp=sum(cat(3,akCoe_old{:}),3); %级数求和: a(K)

resiStressRate=cellfun(@(x,w)w.*x,resiStressRate,weight,'UniformOutput',0);
resiStressRateComp=sum(cat(3,resiStressRate{:}),3); %时间积分

a0CoeComp=-coeBasisCosComp+oldCoeBasisCosComp+resiStressRateComp;

end 
