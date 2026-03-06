function [Low, Upp] = getLowUppValue(xVal, params, NEye, asyincr, asydecr)
% getLowUppValue - 根据变量历史趋势调整MMA渐近上下界
%
% 输入：
%   xval     : 当前设计变量
%   Param    : 参数结构体，包含历史变量和上下界
%   NEye     : 单位向量（ones），大小为 NumVar×1
%   asyincr  : 渐近增加因子
%   asydecr  : 渐近减少因子
%
% 输出：
%   Low, Upp : 更新后的设计变量下界和上界

NumVar     = params.NumVar;
Upp        = params.Upp;
Low        = params.Low;
xold1      = params.rho1;
xold2      = params.rho2;
MaxRelDens = params.MaxRelDens * ones(NumVar, 1);
MinRelDens = params.MinRelDens * ones(NumVar, 1);

% 判断变量变化趋势：连续单调变化 → 增加渐近界； 震荡 → 减小渐近界
factor = NEye;
delta  = (xVal - xold1) .* (xold1 - xold2);
factor(delta > 0) = asyincr;   % 趋势明确 → 放宽
factor(delta < 0) = asydecr;   % 震荡趋势 → 收紧

% 更新上下界
Low = xVal - factor .* (xold1 - Low);
Upp = xVal + factor .* (Upp - xold1);

% 强制控制上下限范围，避免界限过宽或过窄
LowMin = xVal - 10 * (MaxRelDens - MinRelDens);
LowMax = xVal - 0.01 * (MaxRelDens - MinRelDens);
Low    = min(max(Low, LowMin), LowMax);

UppMin = xVal + 0.01 * (MaxRelDens - MinRelDens);
UppMax = xVal + 10   * (MaxRelDens - MinRelDens);
Upp    = max(min(Upp, UppMax), UppMin);

end

% function [Low,Upp]=getLowUppValue(xval,Param,NEye,asyincr,asydecr)
% 
% NumVar=Param.NumVar;
% Upp=Param.Upp; 
% Low=Param.Low;
% xold1=Param.rho1; 
% xold2=Param.rho2;
% MaxRelDens=Param.MaxRelDens*ones(NumVar,1);
% MinRelDens=Param.MinRelDens*ones(NumVar,1);
% 
% factor=NEye;
% delta=(xval-xold1).*(xold1-xold2);
% factor(delta>0)=asyincr;
% factor(delta<0)=asydecr;
% Low=xval-factor.*(xold1-Low);
% Upp=xval+factor.*(Upp-xold1);
% LowMin=xval-10*(MaxRelDens-MinRelDens);
% LowMax=xval-0.01*(MaxRelDens-MinRelDens);
% Low=min(max(Low,LowMin),LowMax);
% UppMin=xval+0.01*(MaxRelDens-MinRelDens);
% UppMax=xval+10*(MaxRelDens-MinRelDens);
% Upp=max(min(Upp,UppMax),UppMin);
% 
% end