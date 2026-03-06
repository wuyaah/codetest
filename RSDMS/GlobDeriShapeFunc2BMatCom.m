function BMatCom=GlobDeriShapeFunc2BMatCom(GlobDeriShapeFunc,EleDim)
% Given is the derivated global shape function at a certain point
% [Nm,x Nm,y Nm,z] returns the component of B matrix as
% [Nm,x   0    0;
%   0    Nm,y  0;
%   0     0   Nm,z;
%  Nm,y  Nm,x  0
%   0    Nm,z Nm,y;
%  Nm,z   0   Nm,x;]

if EleDim==3
    Nm_x=GlobDeriShapeFunc(1);
    Nm_y=GlobDeriShapeFunc(2);
    Nm_z=GlobDeriShapeFunc(3);
    BMatCom=[Nm_x    0     0;
              0     Nm_y   0;
              0      0    Nm_z;
             Nm_y   Nm_x   0;
              0     Nm_z  Nm_y;
             Nm_z    0    Nm_x;];
elseif EleDim==2
    Nm_x=GlobDeriShapeFunc(1);
    Nm_y=GlobDeriShapeFunc(2);
    BMatCom=[Nm_x   0;
             0    Nm_y;
            Nm_y  Nm_x;];
end

end

