function [NodeForceVec,dimension]=getNodeForceVec(LoadData,EleNodeLabelSet,NodeLabelList)

% 默认自由度 (如果找不到力的数据，将使用此维度生成0向量)
% 实体单元通常为3，壳/梁单元通常为6。为了安全起见，默认设为6，可根据需要修改
default_dim = 3; 

try 
    % 情况1：存在 NodeCFSet (集中力)
    if isfield(LoadData, 'NodeCFSet') && ~isempty(LoadData.NodeCFSet)
        NodeForce = LoadData.NodeCFSet;
        NodeForceVec = reshape(transpose(cell2mat(NodeForce)),[],1);
        dimension = size(NodeForce,2);

    % 情况2：不存在 NodeCFSet 但存在 NodeForceSet (分布力/节点力)
    elseif isfield(LoadData, 'NodeForceSet') && ~isempty(LoadData.NodeForceSet)
        NodeForce = LoadData.NodeForceSet;
        LoadNodes = LoadData.NodeForceLabel;
        dimension = size(NodeForce,1);
        
        % 初始化 从abaqus导出的节点力（内力）
        PosNodeForce = cellfun(@(x)ismember(EleNodeLabelSet,x),LoadNodes,'UniformOutput',0);   
        EleNodeForce = cellfun(@(x)cell2mat(squeeze(NodeForce(x,:,:))),num2cell(1:dimension),'UniformOutput',0);
        OriNodeForceVec = cellfun(@(x)cellfun(@(y)sum(x(y)),PosNodeForce),EleNodeForce,'UniformOutput',0);
        OriNodeForceVec = (-1).*cell2mat(OriNodeForceVec);
        
        % 初始化PosNodeForceVec
        PosNodeForceVec = repmat(ismember(NodeLabelList,cell2mat(LoadNodes)),1,dimension);
        PosNodeForceVec = reshape(transpose(PosNodeForceVec),[],1);
        
        % 初始化NodeForceVec（外力）
        NodeForceVec = zeros(numel(PosNodeForceVec),1); 
        NodeForceVec(PosNodeForceVec) = reshape(transpose(OriNodeForceVec),[],1);
        
    else
        % 情况3：两个字段都不存在，手动抛出错误以进入 catch 块处理为 0
        error('No force data found in LoadData');
    end

catch
    % --- 容错处理 ---
    % 如果上述过程出错或数据不存在，令 Force 为 0
    % 使用默认维度生成与节点列表匹配的零向量
    dimension = default_dim;
    NodeForceVec = zeros(numel(NodeLabelList) * dimension, 1);
end

% 清除极小值
NodeForceVec(abs(NodeForceVec)<1e-8)=0;

end