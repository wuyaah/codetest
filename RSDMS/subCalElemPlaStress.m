function plaStress=subCalElemPlaStress(stress,elemYield)


effStress = cellfun(@subCalEffStress, stress);

% Determine plastic coefficients: how much effective stress exceeds yield
plaCoe = effStress > elemYield;
plaCoe = plaCoe.*(effStress-elemYield)./(effStress+eps);

% Calculate plastic stress components by scaling total stress with plastic coefficients
plaStress = cellfun(@(x, Coe)Coe.*x, stress, num2cell(plaCoe), 'UniformOutput', false);

end
