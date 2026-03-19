function [CMatFull,varPosFull]=GetConstrainedGlobalMatSparseList(dataPack,matName,dim,numStrComp,elemNodeLabelSet,nodeLabelList,nodeConstrainedDoF)

elemNGK=cellfun(@(x)x.NGK,dataPack);
nodeLabelCell=num2cell(nodeLabelList);
totalVar=sum(numStrComp.*elemNGK); %the number of components of residual stress
VirtualVar=totalVar+1;
% NGK=max(elemNGK); %临时保留

%% Return element in the given sequence that has the given node
[elemIDSharedbyNode,columInfo]=cellfun(@(x)find(elemNodeLabelSet==x),nodeLabelCell,'UniformOutput',0);
% numElemSharedbyNode=cellfun(@(x)numel(x),elemIDSharedbyNode,'UniformOutput',0);
% maxNumElemSharedByNode=max(cell2mat(numElemSharedbyNode));

%% Cellfun Reshape：GlobalCSparseForm；Iter：代表约束的自由度方向
elemMat2Node=cellfun(@(x)eval(['x.' matName]),dataPack,'UniformOutput',0);
GlobalCSparseForm=cellfun(@(x,o)getGlobalCSparseForm(elemMat2Node(x),o),elemIDSharedbyNode,columInfo,'UniformOutput',0);
if size(nodeConstrainedDoF,2)~=0 %存在位移约束
    for iter=1:size(nodeConstrainedDoF,1) %numel(NodeConstrainedDoF)
        selConstrNode=nodeConstrainedDoF{iter};
        selPos=cellfun(@(x)find(nodeLabelList==x),selConstrNode);
        GlobalCSparseForm(selPos)=cellfun(@(x) ...
            cellfun(@(c)AddConstr2CSparse(c,dim,iter),x,'UniformOutput',0), ...
            GlobalCSparseForm(selPos),'UniformOutput',0);
    end
end
%%

GlobalCSparseForm=cellfun(@(x)cell2mat(transpose(x)),GlobalCSparseForm,'UniformOutput',0);
numRelevantColumn=cellfun(@(x)size(x,2),GlobalCSparseForm,'UniformOutput',0);
maxNumRelevantColumn=max(cell2mat(numRelevantColumn));
ExtendZeroMat=cellfun(@(x)zeros(dim,(maxNumRelevantColumn-x)),numRelevantColumn,'UniformOutput',0);
CMatFull=cellfun(@(x,y)[x,y],GlobalCSparseForm,ExtendZeroMat,'UniformOutput',0);
CMatFull=cell2mat(CMatFull);

VarPos=cellfun(@(x)getSharedKronElemDof(x,elemNGK(x),numStrComp),elemIDSharedbyNode,'UniformOutput',0);
VarPos=cellfun(@(x)repmat(x,dim,1),VarPos,'UniformOutput',0);
VarPosComplementary=cellfun(@(x)VirtualVar.*ones(dim,(maxNumRelevantColumn-x)),numRelevantColumn,'UniformOutput',0);
varPosFull=cellfun(@(x,y)[x,y],VarPos,VarPosComplementary,'UniformOutput',0);
varPosFull=cell2mat(varPosFull);
              
end

%% %%%%%%%% 原始构建代码: GlobalCSparseForm %%%%%%%%
% GlobalCSparseForm=cell(size(NodeLabelList));
% for ii=1:numel(EleShareNode)
%     NodeCSparseForm=[];
%     AssociateEle=EleShareNode{ii};
%     NodeSeqinAssoEle=ColumInfo{ii};
%     for jj=1:1:size(AssociateEle,1)
%         NodeSeq=NodeSeqinAssoEle(jj);
%         NodeCMatComp=eval(['DataPack{AssociateEle(jj)}.' MatName '{NodeSeq}']);
%         NodeCSparseForm=[NodeCSparseForm;NodeCMatComp];
%     end
%     GlobalCSparseForm{ii}=NodeCSparseForm;
% end
% NodeLineConstrainedDoF1=cellfun(@(x)find(NodeLabelList==x),NodeConstrainedDoF{1},'UniformOutput',0);
% NodeLineConstrainedDoF2=cellfun(@(x)find(NodeLabelList==x),NodeConstrainedDoF{2},'UniformOutput',0);
% NodeLineConstrainedDoF3=cellfun(@(x)find(NodeLabelList==x),NodeConstrainedDoF{3},'UniformOutput',0);
% for ii=1:numel(NodeLineConstrainedDoF1)
%     GlobalCSparseForm{NodeLineConstrainedDoF1{ii}}(1:Dimension:end)=0;
% end
% for ii=1:numel(NodeLineConstrainedDoF2)
%     GlobalCSparseForm{NodeLineConstrainedDoF2{ii}}(2:Dimension:end)=0;
% end
% for ii=1:numel(NodeLineConstrainedDoF3)
%     GlobalCSparseForm{NodeLineConstrainedDoF3{ii}}(3:Dimension:end)=0;
% end

 