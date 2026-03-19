function InCompressCoe=getIncompressCoe(sigmaY,numStrComp)
% get coefficients of the incompressible condition for the back stress tensor writen in a form of [U,X]
% InCompressCoe is a 1-by-n vector, which 
% 3D case:  InCompressCoe*[U1;U2;X]=0; 
% CPE case: InCompressCoe*[U1;U2;X]=0;
% CPS case: InCompressCoe*{U1;U2]=0; using U2S matrix and

global OriFuncHandles
ConvertVariableFunc=OriFuncHandles.ConvertVariable;

if numStrComp==6       % 3D
    U2SigMat=ConvertVariableFunc.UX2Sig3D(sigmaY); 
    U2SigMatComp=[U2SigMat(1:3,1:2),U2SigMat(1:3,6)];
elseif numStrComp==4   % CPE
    U2SigMat=ConvertVariableFunc.UX2SigCPE(sigmaY); 
    U2SigMatComp=[U2SigMat(1:3,1:2),U2SigMat(1:3,4)];  
elseif numStrComp==3   % CPS
    U2SigMat=ConvertVariableFunc.UX2SigCPS(sigmaY);
    U2SigMatComp=U2SigMat(1:2,1:2);     
end
InCompressCoe=sum(U2SigMatComp,1);

end
