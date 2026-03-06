function [E,Emin,Yield,YieldMin,Limit,LimitMin,Penal,YPenal,UPenal]=subExtMatInfo(matInfo)

%matInfo={E.*[1,1e-8]; Yield.*[1,1e-8]; Limit.*[1,1e-8]; mu; Penal; YPenal; UPenal};

E=matInfo{1}(1); 
Emin=matInfo{1}(2);

Yield=matInfo{2}(1); 
YieldMin=matInfo{2}(2);

Limit=matInfo{3}(1);
LimitMin=matInfo{3}(2);

Penal=matInfo{5}; 
YPenal=matInfo{6};
UPenal=matInfo{7};

end