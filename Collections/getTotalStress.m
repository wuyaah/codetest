function totalStressCell=getTotalStress(DiffStressCell,BackStressCell,ListHardGauss)
% This function calculates total stress from cell of difference of stress, 
% cell of back stress, and indices list of gaussian points having hardening
% behavior

totalStressCell=DiffStressCell;  % Assign DiffStressCell as initial TotalStress

totalStressBlock=cellfun(@(d,b)d+b,DiffStressCell(ListHardGauss),BackStressCell,'UniformOutput',0);
totalStressCell(ListHardGauss)=totalStressBlock;
 
end



% function totalStressCell=getTotalStress(DiffStressCell,BackStressCell,ListHardGauss)
% % This function calculates total stress from cell of difference of stress, 
% % cell of back stress, and indices list of gaussian points having hardening
% % behavior
% 
% totalStressCell=DiffStressCell;  % Assign DiffStressCell as initial TotalStress
% 
%  for m=1:numel(ListHardGauss)
%      GaussIndex=ListHardGauss(m);
%      totcalStressBlock=DiffStressCell{GaussIndex}+BackStressCell{m};
%      totalStressCell{GaussIndex}=totcalStressBlock;
%  end
%  
% end

