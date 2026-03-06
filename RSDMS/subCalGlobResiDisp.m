function gResiDispRate = subCalGlobResiDisp(nodeForceRate, globKMat, constrainedDoFVec)

nodeForceRateMat = cell2mat(nodeForceRate);
gResiDispRate = getNodeDispVecNew(globKMat, nodeForceRateMat, constrainedDoFVec);

end
