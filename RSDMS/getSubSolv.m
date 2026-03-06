function [xmma, ymma, zmma, lamma, xsimma, etamma, mumma, zetmma, smma] = getSubSolv(epsimin, alpha, beta, P0, Q0, P, Q, params)
% =========================================================================
% getSubSolv.m
% 用于求解 MMA（移动渐近线法）中的凸子问题
%
% 输入参数：
%   - Epsimin: 双重残差的最小阈值（收敛准则）
%   - alpha, belta: 每个变量的下/上渐近线
%   - P0, Q0, P, Q: 目标函数和约束的一阶近似矩阵
%   - params: 参数结构体，包含：
%       * NumVar: 变量个数
%       * NumConstr: 约束个数
%       * Low, Upp: 变量下限和上限
%       * a, b, a0, d, Cmma: MMA 内部参数
%
% 输出参数：
%   - xmma: 最优解（主变量）
%   - ymma, zmma: 双变量
%   - lamma: 拉格朗日乘子
%   - xsimma, etamma, mumma, zetmma, smma: KKT 条件相关变量
% =========================================================================
% 解包参数
NumVar    = params.NumVar;
NumConstr = params.NumConstr;
Upp       = params.Upp;
Low       = params.Low;
a         = params.a;
b         = params.b;
d         = params.d;
Cmma      = params.Cmma;
a0        = params.a0;

% 初始化参数
Epsi = 1;                       % 初始正则化项
NEye = ones(NumVar, 1);         % N维全1向量
MEye = ones(NumConstr, 1);      % M维全1向量

% 初始猜测解
x    = 0.5 * (alpha + beta);   % 初始化主变量为中点
y    = MEye;
z    = 1;
lam  = MEye;
xsi  = max(NEye, 1 ./ (x - alpha));         % 避免除以0
eta  = max(NEye, 1 ./ (beta - x));
mu   = max(MEye, 0.5 * Cmma);
zet  = 1;
s    = MEye;

Iter = 0;
%%
while Epsi > epsimin
    
    epsvecn = Epsi*NEye;
    epsvecm = Epsi*MEye;
    
    % 计算变量间距和倒数
    ux1 = Upp - x;    xl1 = x - Low;
    ux2 = ux1.^2;     xl2 = xl1.^2;
    uxinv1 = NEye ./ ux1;
    xlinv1 = NEye ./ xl1;

    % 一阶导数项
    plam = P0 + P' * lam;
    qlam = Q0 + Q' * lam;
    gvec = P * uxinv1 + Q * xlinv1;
    dpsidx = plam ./ ux2 - qlam ./ xl2;

    % 构造 KKT 残差向量
    rex = dpsidx - xsi + eta;
    rey = Cmma + d .* y - mu - lam;
    rez = a0 - zet - a' * lam;
    relam = gvec - a * z - y + s - b;
    rexsi = xsi .* (x - alpha) - epsvecn;
    reeta = eta .* (beta - x) - epsvecn;
    remu = mu .* y - epsvecm;
    rezet = zet * z - Epsi;
    res = lam .* s - epsvecm;
    
    % 合并残差
    residu1 = [rex; rey; rez];
    residu2 = [relam; rexsi; reeta; remu; rezet; res];
    residu = [residu1; residu2];
    residunorm = norm(residu);
    residumax = max(abs(residu));
    
    Iter2 = 0;

    % 判断是否满足残差阈值或达到最大迭代
    while residumax > 0.9 * Epsi && Iter2 < 100

        Iter  = Iter  + 1;
        Iter2 = Iter2 + 1;
        
        ux1 = Upp-x;
        xl1 = x-Low;
        ux2 = ux1.*ux1;
        xl2 = xl1.*xl1;
        ux3 = ux1.*ux2;
        xl3 = xl1.*xl2;
        uxinv1 = NEye./ux1;
        xlinv1 = NEye./xl1;
        uxinv2 = NEye./ux2;
        xlinv2 = NEye./xl2;
        plam = P0 + P'*lam ;
        qlam = Q0 + Q'*lam ;
        gvec = P * uxinv1 + Q * xlinv1;
        GG = P * spdiags(uxinv2,0,NumVar,NumVar) - Q * spdiags(xlinv2,0,NumVar,NumVar);
        dpsidx = plam./ux2 - qlam./xl2 ;
        delx = dpsidx - epsvecn./(x-alpha) + epsvecn./(beta-x);
        dely = Cmma + d.*y - lam - epsvecm./y;
        delz = a0 - a'*lam - Epsi/z;
        dellam = gvec - a*z - y - b + epsvecm./lam;
        diagx = plam./ux3 + qlam./xl3;
        diagx = 2*diagx + xsi./(x-alpha) + eta./(beta-x);
        diagxinv = NEye./diagx;
        diagy = d + mu./y;
        diagyinv = MEye./diagy;
        diaglam = s./lam;
        diaglamyi = diaglam+diagyinv;
        
        % --------- 求解线性系统 ---------
        if NumConstr < NumVar
             % 少约束情形（更常见）
            blam = dellam + dely./diagy - GG*(delx./diagx);
            bb = [blam' delz]';
            Alam = spdiags(diaglamyi,0,NumConstr,NumConstr) + GG*spdiags(diagxinv,0,NumVar,NumVar)*GG';
            AA = [Alam , a ; a' , -zet/z ];
            solut = AA\bb;
            dlam = solut(1:NumConstr);
            dz = solut(NumConstr+1);
            dx = -delx./diagx - (GG'*dlam)./diagx;
        else
            % 多约束情形
            diaglamyiinv = MEye./diaglamyi;
            dellamyi = dellam + dely./diagy;
            Axx = spdiags(diagx,0,NumVar,NumVar) + GG'*spdiags(diaglamyiinv,0,NumConstr,NumConstr)*GG;
            azz = zet/z + a'*(a./diaglamyi);
            axz = -GG'*(a./diaglamyi);
            bx = delx + GG'*(dellamyi./diaglamyi);
            bz = delz - a'*(dellamyi./diaglamyi);
            AA = [Axx,axz;axz',azz];
            bb = [-bx' -bz]';
            solut = AA\bb;
            dx = solut(1:NumVar);
            dz = solut(NumVar+1);
            dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
        end
        
        % 更新其他变量方向
        dy   = -dely./diagy + dlam./diagy;
        dxsi = -xsi + epsvecn./(x-alpha) - (xsi.*dx)./(x-alpha);
        deta = -eta + epsvecn./(beta-x) + (eta.*dx)./(beta-x);
        dmu  = -mu + epsvecm./y - (mu.*dy)./y;
        dzet = -zet + Epsi/z - zet*dz/z;
        ds   = -s + epsvecm./lam - (s.*dlam)./lam;

        % --------- 线搜索步长 ---------
        xx  = [ y' z lam' xsi' eta' mu' zet s']';
        dxx = [dy' dz dlam' dxsi' deta' dmu' dzet ds']';
        
        stepxx = -1.01*dxx./xx;
        stmxx  = max(stepxx);
        stepalfa = -1.01*dx./(x-alpha);
        stmalfa = max(stepalfa);
        stepbeta = 1.01*dx./(beta-x);
        stmbeta = max(stepbeta);
        stmalbe  = max(stmalfa,stmbeta);
        stmalbexx = max(stmalbe,stmxx);
        stminv = max(stmalbexx,1);
        steg = 1/stminv;
        
        % 保存旧变量
        xold   =  x;
        yold   =  y;
        zold   =  z;
        lamold =  lam;
        xsiold =  xsi;
        etaold =  eta;
        muold  =  mu;
        zetold =  zet;
        sold   =  s;
        
        % 反复缩小步长直到残差减小
        Iter3 = 0;
        resinew = 2*residunorm;
        while resinew > residunorm && Iter3 < 50

            Iter3 = Iter3+1;
            
            x   =   xold + steg*dx;
            y   =   yold + steg*dy;
            z   =   zold + steg*dz;
            lam = lamold + steg*dlam;
            xsi = xsiold + steg*dxsi;
            eta = etaold + steg*deta;
            mu  = muold  + steg*dmu;
            zet = zetold + steg*dzet;
            s   =   sold + steg*ds;

            % 重新计算残差
            ux1 = Upp-x;
            xl1 = x-Low;
            ux2 = ux1.*ux1;
            xl2 = xl1.*xl1;
            uxinv1 = NEye./ux1;
            xlinv1 = NEye./xl1;
            plam = P0 + P'*lam ;
            qlam = Q0 + Q'*lam ;
            gvec = P*uxinv1 + Q*xlinv1;
            dpsidx = plam./ux2 - qlam./xl2 ;
            
            rex = dpsidx - xsi + eta;
            rey = Cmma + d.*y - mu - lam;
            rez = a0 - zet - a'*lam;
            relam = gvec - a*z - y + s - b;
            rexsi = xsi.*(x-alpha) - epsvecn;
            reeta = eta.*(beta-x) - epsvecn;
            remu = mu.*y - epsvecm;
            rezet = zet*z - Epsi;
            res = lam.*s - epsvecm;
            
            residu1 = [rex' rey' rez]';
            residu2 = [relam' rexsi' reeta' remu' rezet res']';
            residu = [residu1' residu2']';
            resinew = sqrt(residu'*residu);
            steg = steg/2;
        end

    residunorm = resinew;
    residumax  = max(abs(residu));
    
    end

    % 更新 Epsi（减小正则化项）
    Epsi = 0.1*Epsi;

end

% 返回解
xmma   = x;
ymma   = y;
zmma   = z;
lamma  = lam;
xsimma = xsi;
etamma = eta;
mumma  = mu;
zetmma = zet;
smma   = s;

end
