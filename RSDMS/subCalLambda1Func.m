function [deltaLambda1,plaStress_new]=subCalLambda1Func(info,deltaLambda1,totalStress,plaStress)
%% matInfo
E  = info.IntE;       % Young's modulus
nu = info.Poisson;    % Poisson's ratio
% Lame constant value
mu=E/(2*(1+nu));      % shear modulus
lambda=(E*nu)/((1+nu)*(1-2*nu));
elaMatrix_Prime=2*mu*eye(3)+lambda*ones(3); % 主应力空间下的等效弹性矩阵 Eij

%%
pos=cellfun(@(x)norm(x)>0,plaStress,'UniformOutput',1);
[stress,dir]=cellfun(@(x)getEigenStress(x),totalStress(pos),'UniformOutput',0);% dir{1}*diag(stress{1})*dir{1}'
gradFlowDir=cellfun(@(x)subFlowDirect(x),stress,'UniformOutput',0);
deltaStress=cellfun(@(x,d)subCalDeltaStress(x,info.Yield),stress,'UniformOutput',0);
eigenStress=cellfun(@(x)elaMatrix_Prime*x,gradFlowDir,'UniformOutput',0);

deltaLambda1_new=zeros(size(deltaLambda1));
deltaLambda1_new(pos)=cellfun(@(d,x)mean(d./x),deltaStress,eigenStress,'UniformOutput',1);
% deltaLambda1_new(pos)=cellfun(@(d,x)d./x,deltaStress,eigenStress,'UniformOutput',0);

deltaLambda1=deltaLambda1+deltaLambda1_new; 

plaStress_pos_tensor=cellfun(@(x,d)d*diag(x)*d',deltaStress,dir,'UniformOutput',0);
switch info.numStrComp
    case 6
        plaStress_pos=cellfun(@(x)[x(1,1);x(2,2);x(3,3);x(1,2);x(1,3);x(2,3)],plaStress_pos_tensor,'UniformOutput',0);
    case 4
        plaStress_pos=cellfun(@(x)[x(1,1);x(2,2);x(3,3);x(1,2)],plaStress_pos_tensor,'UniformOutput',0);
    case 3
        plaStress_pos=cellfun(@(x)[x(1,1);x(2,2);x(1,2)],plaStress_pos_tensor,'UniformOutput',0);
end
plaStress_new=plaStress;
plaStress_new(pos)=plaStress_pos;

end

