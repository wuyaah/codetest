function [GPoint,weight]=GaussLegendreXW(Num,a,b)
% GaussLegendreXW: 计算Legendre高斯点(GPoint)和权重(w)
% Num: 节点数; a,b: 积分区间

Num=Num-1;
N1=Num+1; 
N2=Num+2;

x=transpose(linspace(-1,1,N1));
y=cos((2*(0:Num)'+1)*pi/(2*Num+2))+(0.27/N1)*sin(pi*x*Num/N2);

L=zeros(N1,N2);
y0=2;

% Newton法迭代求解Legendre多项式根
while max(abs(y-y0))>eps
    L(:,1)=1; %P0
    L(:,2)=y; %P1

    for k=2:N1
        L(:,k+1)=((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k;
    end
    
    Lp=N2*(L(:,N1)-y.*L(:,N2))./(1-y.^2); % 一阶导数
    y0=y;
    y=y0-L(:,N2)./Lp;
end

GPoint=y;
GPoint=(a*(1-GPoint)+b*(1+GPoint))/2; %%% ((b-a)*GPoint+(b+a))/2y为GPoint,区间变换为[a,b]
weight=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

[GPoint,index]=sort(GPoint);
weight=weight(index);

% disp(['Origial Cyclic Loading points: ', sprintf('%8.4f', sort(y)'), newline, ...
%   'Current Gauss-Legendre points: ', sprintf('%8.4f', GPoint'), newline, ...
%   'The real weight factors (w/2): ', sprintf('%8.4f', weight')]);

end

