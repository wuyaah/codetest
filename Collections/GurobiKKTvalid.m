function kktErr = GurobiKKTvalid(model,alphaPos,LinConstEnrich,QuadIndicesCell,result)
% ---------------------------------------------------------------
% GurobiKKTvalid  (with pi / qcpi existence check)
% ---------------------------------------------------------------

% ========== 检查 result.pi 和 result.qcpi 是否存在 ==========
if ~isfield(result,'pi')  || isempty(result.pi) || ...
   ~isfield(result,'qcpi') || isempty(result.qcpi)

    warning('result.pi 或 result.qcpi 不存在，跳过 KKT 验证计算。');

    % 返回 NaN 结构体
    kktErr.stationarity        = NaN;
    kktErr.maxStationarityErr  = NaN;
    kktErr.qcSlack             = NaN;
    kktErr.ineqComplement      = NaN;
    kktErr.eqComplement        = NaN;
    kktErr.PlaDiss             = NaN;
    return;
end
% ===============================================================


% ========== 基础维度 ==========
numVar     = size(LinConstEnrich,2);
numU1Comp  = numel(QuadIndicesCell{1});
numIneqs   = numel(QuadIndicesCell);

% ========== 拉格朗日乘子 ==========
lambda     = result.qcpi(:);      % 二次约束乘子 λ ≥ 0
eqLambda   = result.pi(:);        % 等式约束乘子 μ

% ========== 锥/二次约束相关数据 ==========
YieldDeg        = result.ProcessedRes.Yield(:);  % 多项式出力量
vertUlist       = result.vertUlist(:);           % 顶点 U 值
GradVertUlist   = 2 .* vertUlist;                % 对 U^2 的导数 = 2U

% ========== 目标函数梯度（线性部分） ==========
objFuncMat = model.obj(:);

% ========== 构造二次 Jacobian ==========
jacoColumn = reshape(repmat(1:numIneqs, numU1Comp, 1), [], 1);
jacoRow    = reshape(transpose(cell2mat(QuadIndicesCell)), [], 1);
jacoConst  = sparse(jacoRow, jacoColumn, GradVertUlist, numVar, numIneqs);

% ========== KKT：梯度驻留 ==========
stationarityErr = -objFuncMat + jacoConst * lambda + LinConstEnrich' * eqLambda;
maxStationarityErr = max(abs(stationarityErr));

% ========== 二次约束 slack consistency ==========
qcSlackErr = max(abs(1 - (result.qcslack + YieldDeg.^2)));

% ========== 不等式互补性 ==========
complementarityIneq = max(abs(lambda .* (1 - YieldDeg.^2)));

% ========== 等式互补性 ==========
eqConstraintValue = LinConstEnrich * result.x;
complementarityEq = max(abs(eqLambda .* eqConstraintValue));

PlaDiss = -1 + eqLambda' * LinConstEnrich(:,alphaPos);

% ========== 打包结构体输出 ==========
kktErr.stationarity        = stationarityErr;
kktErr.maxStationarityErr  = maxStationarityErr;
kktErr.qcSlack             = qcSlackErr;
kktErr.ineqComplement      = complementarityIneq;
kktErr.eqComplement        = complementarityEq;
kktErr.PlaDiss             = PlaDiss;

end
