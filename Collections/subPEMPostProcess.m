function result=subPEMPostProcess(stress,elemRsiStress,elemYield)

totalStressVec=cellfun(@cell2mat,stress,'UniformOutput',0);
resiStressVec=cellfun(@cell2mat,elemRsiStress,'UniformOutput',0);
effTotalStress=cellfun(@(x)cellfun(@subCalEffStress,x),stress,'UniformOutput',0);
effResiStress=cellfun(@(x)cellfun(@subCalEffStress,x),elemRsiStress,'UniformOutput',0);
Yield=cellfun(@(x,y)x./y,effTotalStress,elemYield,'UniformOutput',0);


% Store processed results in a structured format
ProcessedRes.TotalStress       = cell2mat(totalStressVec);
ProcessedRes.ResidualStress    = cell2mat(resiStressVec);
ProcessedRes.EffTotalStress    = cell2mat(effTotalStress);
ProcessedRes.EffResidualStress = cell2mat(effResiStress);
ProcessedRes.Yield             = cell2mat(Yield);
ProcessedRes.MaxYield          = max(max(ProcessedRes.Yield));

% Assign to output
result.ProcessedRes  = ProcessedRes;

end
