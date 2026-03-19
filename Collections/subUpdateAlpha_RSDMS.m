function alpha=subUpdateAlpha_RSDMS(alphaSet,coePhiAlpha,outerLoop)
%% Descending of the procedure towards shakedown

Damp=0.5;
if outerLoop>2
    factor=(alphaSet(end)-alphaSet(end-1))/(coePhiAlpha(end)-coePhiAlpha(end-1));
    factor=Damp*factor;
else
    factor=1;
end
alpha=alphaSet(end)-factor*coePhiAlpha(end);

%% Descending arithmetic sequence, calculate "coePhi"
% alpha=alphaSet(end)-0.01;

end 
