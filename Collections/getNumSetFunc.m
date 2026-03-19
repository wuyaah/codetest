function NumSet=getNumSetFunc(parameters,numVert)

NumSet.numVert=numVert;
NumSet.numElem=parameters{1};
NumSet.numNode=parameters{2};
% NumSet.NGKSet=parameters{3};
NumSet.NGK=parameters{3}(1);
NumSet.NGKSet=parameters{3};
NumSet.numStrComp=parameters{4};
NumSet.Dimension=parameters{7};
NumSet.NumTotalGauss=sum(NumSet.NGKSet);
NumSet.NumU1Comp=NumSet.numStrComp-1;                      % Number of [U] components at each gaussian point
NumSet.NumX1Comp=1;                       

if NumSet.numStrComp==3 % Plane stress case
    NumSet.NumU1Comp=NumSet.numStrComp;
    NumSet.NumX1Comp=0;
end

NumSet.NumU1=NumSet.NumTotalGauss*NumSet.NumU1Comp;        % Number of total [U] components
NumSet.NumX1=NumSet.NumTotalGauss*NumSet.NumX1Comp;        % Number of total [X] components

NumSet.totalVar=NumSet.numVert*NumSet.NumU1+NumSet.NumX1+2;     % Number of total variables
NumSet.numEqs=NumSet.Dimension*NumSet.numNode+(NumSet.numVert-1)*NumSet.NumTotalGauss*NumSet.NumU1Comp;
NumSet.numIneqs=NumSet.numVert*NumSet.NumTotalGauss;

end

