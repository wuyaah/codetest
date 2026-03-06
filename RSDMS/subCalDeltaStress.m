function deltaSig=subCalDeltaStress(stress,Yield)

hydro=mean(stress);
devi=stress-hydro;

% von Mises 等效应力(主应力空间)
eqStress=sqrt(0.5*((stress(1)-stress(2))^2+(stress(2)-stress(3))^2+(stress(1)-stress(3))^2));

% 偏应力缩至屈服面
sigYPrime=Yield/eqStress*devi;
sigY=sigYPrime+hydro;

% deltaSigma（主应力空间）
deltaSig=stress-sigY;

end


% 
% function res=subCalDeltaStress_Principal(totalStress,Yield,numStrComp)
% %% step1. Voigt -> 张量形式
% switch numStrComp
%     case 6 % 3D
%         sigma = [ totalStress(1) totalStress(4) totalStress(5);
%                   totalStress(4) totalStress(2) totalStress(6);
%                   totalStress(5) totalStress(6) totalStress(3) ];
%     case 4 % 平面应变
%         sigma = [ totalStress(1) totalStress(4) 0;
%                   totalStress(4) totalStress(2) 0;
%                   0              0             totalStress(3) ];
%     case 3 % 平面应力
%         sigma = [ totalStress(1) totalStress(3) 0;
%                   totalStress(3) totalStress(2) 0;
%                   0              0             0 ];
%     otherwise
%         error('numStrComp must be 3, 4, or 6.');
% end
% 
% %% step2. 求主应力与主方向
% [dirMat,stressVals]=eig(sigma);
% 
% %% step3. 在主应力空间计算屈服面投影
% stress=diag(stressVals);
% hydro=mean(stress);
% devi=stress-hydro;
% 
% % von Mises 等效应力(主应力空间)
% eqStress=sqrt(0.5*((stress(1)-stress(2))^2+(stress(2)-stress(3))^2+(stress(1)-stress(3))^2));
% 
% % 偏应力缩至屈服面
% sigYPrime=Yield/eqStress*devi;
% sigY=sigYPrime+hydro;
% 
% % deltaSigma（主应力空间）
% deltaSig_Principal=stress-sigY;
% 
% %% step4. 屈服函数梯度方向检查(可选)
% % gradFlow = (1/(2*eqStress)) * [2*stress(1)-stress(2)-stress(3);
% %                                2*stress(2)-stress(1)-stress(3);
% %                                2*stress(3)-stress(1)-stress(2)];
% % cosTheta = dot(deltaSig_Principal, gradFlow)/(norm(deltaSig_Principal)*norm(gradFlow));
% 
% %% step5. 回转到原坐标系
% % sigma_tensor   = dirMat * diag(stress) * dirMat';
% % sigmaY_tensor = dirMat * diag(sigY) * dirMat';
% % deltaTensor = sigma_tensor - sigmaY_tensor;
% deltaTensor=dirMat*diag(deltaSig_Principal)*dirMat'; 
% 
% %% step6. 张量形式 -> Voigt
% switch numStrComp
%     case 6
%         deltaStress=[deltaTensor(1,1);
%                      deltaTensor(2,2);
%                      deltaTensor(3,3);
%                      deltaTensor(1,2);
%                      deltaTensor(1,3);
%                      deltaTensor(2,3)];
%     case 4
%         deltaStress=[deltaTensor(1,1);
%                      deltaTensor(2,2);
%                      deltaTensor(3,3);
%                      deltaTensor(1,2)];
%     case 3
%         deltaStress=[deltaTensor(1,1);
%                      deltaTensor(2,2);
%                      deltaTensor(1,2)];
% end
% 
% %%
% res=deltaSig_Principal;
% % res=deltaStress;
% 
% 
% end

