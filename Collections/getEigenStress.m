function [stressVals,dir]=getEigenStress(stress)
%% get the stress tensor
numStrComp=numel(stress);
switch numStrComp
    case 6 % 3D
        sigma=[stress(1) stress(4) stress(5);
               stress(4) stress(2) stress(6);
               stress(5) stress(6) stress(3)];
    case 4 % 平面应变
        sigma=[stress(1) stress(4) 0;
               stress(4) stress(2) 0;
               0         0         stress(3)];
    case 3 % 平面应力
        sigma=[stress(1) stress(3) 0;
               stress(3) stress(2) 0;
               0         0         0];
    otherwise
        error('numStrComp must be 3, 4, or 6.');
end

%% get the eigen value
[dir,stressVals]=eig(sigma);
stressVals=diag(stressVals);

end