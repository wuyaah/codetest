function nodeCSparseForm=getGlobalCSparseForm(elemMat,ColumInfo)

nodeSeqinAssoElem=num2cell(ColumInfo);
nodeCSparseForm=cellfun(@(x,n)x{n},elemMat,nodeSeqinAssoElem,'UniformOutput',0);
% nodeCSparseForm=cell2mat(nodeCSparseForm);

% NodeCSparseForm=[];
% for iter=1:numElemShared
%     NodeSeqinAssoElem=ColumInfo(iter);
%     NodeCMatComp=SelMatName{iter}{NodeSeqinAssoElem};
%     NodeCSparseForm=[NodeCSparseForm;NodeCMatComp]; %#ok<AGROW>
% end

end