function [EffCalFunc,GradYieldFunc,UConvertSig,SigConvertU]=getConvertFunc(numStrComp)

global OriFuncHandles

if isempty(OriFuncHandles)
    OriFuncHandles=load('Ori.FuncHandles.mat');
end
GetDiffYieldFunc=OriFuncHandles.GetDiffYieldFunc;
GetEquivalent=OriFuncHandles.GetEquivalent;
ConvertVariable=OriFuncHandles.ConvertVariable;

if numStrComp==6       % 3D Model
    GradYieldFunc=GetDiffYieldFunc.GetDiffYield3D;
    EffCalFunc=GetEquivalent.GetEquivStress3D;
    UConvertSig=ConvertVariable.UX2Sig3D; 
    SigConvertU=ConvertVariable.Sig2UX3D;
elseif numStrComp==4   % Plane Strain
    GradYieldFunc=GetDiffYieldFunc.GetDiffYieldCPE;
    EffCalFunc=GetEquivalent.GetEquivStressCPE;
    UConvertSig=ConvertVariable.UX2SigCPE; 
    SigConvertU=ConvertVariable.Sig2UXCPE;
elseif numStrComp==3   % Plane Stress 
    GradYieldFunc=GetDiffYieldFunc.GetDiffYieldCPS;
    EffCalFunc=GetEquivalent.GetEquivStressCPS;
    UConvertSig=ConvertVariable.UX2SigCPS;
    SigConvertU=ConvertVariable.Sig2UXCPS;
end

end