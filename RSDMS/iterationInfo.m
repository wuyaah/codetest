function iterationInfo(VolumeFracAttainedSet,elaLimitSet,ObjSet,time,Change)

if nargin<5
    Change=0;
end
Iter=numel(ObjSet);

t=toc(time);
h=fix(t/3600);m=fix(mod(t,3600)/60);sec=fix(mod(mod(t,3600),60));
if ~isempty(elaLimitSet)
    fprintf('it.: %3i;  obj.: %.4f;  vol.: %.3f;  sd/ela.: %.3f;  ch.: %.5f;  time: %d:%d:%d\n', ...
         Iter, ObjSet(Iter), VolumeFracAttainedSet(Iter), ObjSet(Iter)/elaLimitSet(Iter), Change, h, m, sec);
else
    fprintf('it.: %3i;  obj.: %.4f;  vol.: %.3f;  ch.: %.5f;  time: %d:%d:%d\n', ...
         Iter, ObjSet(Iter), VolumeFracAttainedSet(Iter), Change, h, m, sec);
end
% disp

end