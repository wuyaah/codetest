function [Lambda,deltaLambda]=Paper_SortLambda2Ele(result)

numSet=result.numSet;
numElem=numSet.numElem;
NGK=numSet.NGK;

deltaLambda1Set=result.deltaLambda1Set;

deltaLambda=deltaLambda1Set{end};
deltaLambda=mat2cell(deltaLambda,NGK*ones(numElem,1),size(deltaLambda,2));

elemLambda=sum(cat(3,deltaLambda1Set{:}),3);
Lambda=mat2cell(elemLambda,NGK*ones(numElem,1),size(elemLambda,2));

end

% posZero=(lambda==0);
% numZero=numel(find(posZero==1));
% randZero=min(lambda(~posZero))*rand(numZero,1)/100;
% lambda(posZero)=randZero;
