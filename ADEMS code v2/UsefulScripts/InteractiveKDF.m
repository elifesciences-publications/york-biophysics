%Plots the kernal density function of the variable AllI with a slider to
%manipulate the bandwidth

function InteractiveKDF(AllI)
[AllIdens,AllIx,bw] = ksdensity(AllI,'npoints',1000);
hplot=plot(AllIx, AllIdens,'r');
plotXrange=xlim;
titletext=strcat('KDF with bandwidth=',num2str(bw));
title(titletext)
h = uicontrol('style','slider','min', 0,'max',1,'units','pixel','position',[40 40 1000 23]);
addlistener(h,'ActionEvent',@(hObject, event) makeplot(hObject, event,AllI,hplot,bw,plotXrange));
function makeplot(hObject,event,AllI,hplot,bw,plotXrange)
n = get(hObject,'Value');
bw2=bw*n+0.001;
[AllIdens,AllIx] = ksdensity(AllI, 'bandwidth',bw2,'npoints',1000);
set(hplot,'xdata',AllIx);
set(hplot,'ydata',AllIdens);
xlim(plotXrange);
titletext=strcat('KDF with bandwidth=',num2str(bw2));
title(titletext)
drawnow;
