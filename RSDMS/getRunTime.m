function getRunTime(time)
% a=datetime('now');

t=toc(time);
h=fix(t/3600);
m=fix(mod(t,3600)/60);
sec=fix(mod(mod(t,3600),60));
fprintf('Current calculation time cost: %d:%d:%d\n',h,m,sec);

end
