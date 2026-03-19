function gradFlow=subFlowDirect(stress)
%% von Mises 等效应力(主应力空间)
eqStress=sqrt(0.5*((stress(1)-stress(2))^2+(stress(2)-stress(3))^2+(stress(1)-stress(3))^2));

%% flow direct
gradFlow=(1/(2*eqStress))*[2*stress(1)-stress(2)-stress(3);
                           2*stress(2)-stress(1)-stress(3);
                           2*stress(3)-stress(1)-stress(2)];

end