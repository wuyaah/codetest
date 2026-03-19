function linConstraint=getEqConstrBesideAlphaCol(W1Data,numSet)
% Generate the linear constraint matrix with respect to the number of load
% vertices, Note that in this matrix the column corresponds to alpha is not 
% complete, it does not include theta values (difference between vertices)
% Therefore this column has to be included later!
% Required inputs:
% ABMatSparse:  AB matrix in a sparse form
% NumVert: Number of Vertices
% NumU1:   Number of variable U1
% W1Data:  C*sigma^E (W1data)
%%
%  A  B -W  0   0  0       U  
%[ I  0  0 -I   0  0 ][    X    ] = [0]
%  I  0  0  0   0  0     Alpha

global ABMatSparse

NumU1=numSet.NumU1;
NumVert=numSet.numVert;
linConstUX1=[ABMatSparse,-W1Data];
numTotalU1XAlpha=size(linConstUX1,2); % Number of [U1,X,alpha]

if NumVert==1
    linConstraint=linConstUX1; % Limit analysis
else
    
    [URow,UCol]=find(speye(NumU1));

    % Get list of row number
    URowIncrement=NumU1*(0:(NumVert-2));
    URowIncrementFull=kron(URowIncrement',ones(2*NumU1,1)); 
    URowBasic=repmat(URow,2*(NumVert-1),1);  %Required rows to specify equality condition for each vertex
    URowFull=URowIncrementFull+URowBasic;

    % Get list of column number
    U1Pos=repmat([1,0],1,(NumVert-1));
    UNColumn=reshape((numTotalU1XAlpha+1:numTotalU1XAlpha+(NumVert-1)*NumU1)',NumU1,(NumVert-1));
    UColunmFullMat=kron(U1Pos,UCol)+kron(UNColumn,[0 1]);
    UColunmFull=reshape(UColunmFullMat,size(UColunmFullMat,1)*size(UColunmFullMat,2),1);

    % Get list of value 
    ValuePos=repmat([1,-1],1,(NumVert-1));
    ValueFullMat=kron(ValuePos,ones(NumU1,1));
    ValueFull=reshape(ValueFullMat,size(ValueFullMat,1)*size(ValueFullMat,2),1);
    ULinConst=sparse(URowFull,UColunmFull,ValueFull);

    % LinConstUX1
    %NumTotalVar=size(ULinConst,2);
    %[LinConstUX1Row,LinConstUX1Co,LinConstUX1Val]=find(linConstUX1);
    %LinConstUX1Ext=sparse(LinConstUX1Row,LinConstUX1Co,LinConstUX1Val,size(linConstUX1,1),NumTotalVar);
    LinConstUX1Ext=[linConstUX1,sparse(size(linConstUX1,1),size(ULinConst,2)-size(linConstUX1,2))];
    linConstraint=[LinConstUX1Ext;ULinConst];

end

end