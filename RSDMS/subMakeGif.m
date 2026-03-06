function subMakeGif(fileName,iter)

selFrame=getframe(gcf);  
rgbFrame=frame2im(selFrame);  
[indexFrame,coloredMap]=rgb2ind(rgbFrame,256);  
if iter==1
    imwrite(indexFrame,coloredMap,fileName,'gif','Loopcount',inf,'DelayTime',0.1);
else
    imwrite(indexFrame,coloredMap,fileName,'gif','WriteMode','append','DelayTime',0.1);
end

end