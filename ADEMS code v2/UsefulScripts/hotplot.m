%% HOTPLOT
% creates a heatmap of scatterplot data
% X,Y 1D vectors containing data normally plotted with a scatter plot
% binSize sets size bin in a 2D histogram of the data
% Resolution sets the number of points to interpolate over to make a nice smooth map
% SetbinCentres use this to set the location of the bins, usually so you
% can plot on same axes (can use binCentres output) OPTIONAL
% heatmapdata is the matrix of intensity values
% binCentres marks the location of the bins


function [heatmapdata,binCentres]=hotplot(X,Y, binSize,Resolution,SetbinCentres)

nucDS=cat(1,X,Y);
nucDS=nucDS';
if exist('SetbinCentres')==0
    [nucHot nucXY]=hist3(nucDS,[binSize binSize]);
else
    [nucHot nucXY]=hist3(nucDS,'Ctrs',SetbinCentres);
end
binCentres=nucXY;
nucXYq2=min(nucXY{2}):(max(nucXY{2})-min(nucXY{2}))/(Resolution-1):max(nucXY{2});
nucXYq1=min(nucXY{1}):(max(nucXY{1})-min(nucXY{1}))/(Resolution-1):max(nucXY{1});
[Xq,Yq]=meshgrid(nucXYq2,nucXYq1);
[X,Y]=meshgrid(nucXY{2},nucXY{1});
nucHotFine=interp2( X, Y, nucHot,Xq, Yq);
heatmapdata=nucHotFine;
surf(Xq,Yq,nucHotFine,'EdgeColor','none')
zlim([-10000,10000])
xlim([min(min(Xq)),max(max(Xq))])
ylim([min(min(Yq)),max(max(Yq))])
view(2)
end