function Param=SetmmaParams(NumOPTEle,OriRelDensSet,MinRelDens)

Param.a=0; 
Param.d=0;
Param.a0=1;

Param.NumConstr=1;      % Constr of VolumeFrac
Param.NumVar=NumOPTEle; % Num of Vars（rho）
Param.MaxRelDens=1;     % Default
Param.MinRelDens=MinRelDens;
Param.Low=zeros(Param.NumVar,1); 
Param.Upp=zeros(Param.NumVar,1);
Param.Cmma=1000*ones(Param.NumConstr,1);
Param.rho2=ones(Param.NumVar,1);
Param.rho1=OriRelDensSet;

end

