% Reads in spots array and plots things, 
[TifFile TifDir]=uigetfile('*.tif','Select the frame average');
image_data=imread(strcat(TifDir,TifFile));
[DataFile DataDir]=uigetfile(strcat(TifDir,'*.mat'),'Select the correct tracking data');
 load(strcat(DataDir,DataFile),'SpotsCh1','SpotsCh2');
%image_data=imread(strcat(DataDir,'frame_average.tif'));
figure;
image_data=frame_average;
Spots=SpotsCh1;
%Spots=spotsBaseline;
clear spot_means stdx stdy
% subplot(2,3,2)
% imshow(image_data,[])
% title('Frame average with isolated spots superimposed')
% hold on
% plot(Spots(Spots(:,10)==0,1),Spots(Spots(:,10)==0,2),'.','color','r')
% subplot(2,3,1)
imshow(image_data,[])
title('Found spots superimposed on image data')
hold on
plot(Spots(Spots(:,10)==0,1),Spots(Spots(:,10)==0,2),'.','color','r')
% This seems to be a way of setting the size of reconstruction
ScaleFactor=4;
GaussFrame=zeros(ScaleFactor*size(image_data));
[Xpos,Ypos] = meshgrid(1:ScaleFactor*size(image_data,2),1:ScaleFactor*size(image_data,1));
PSFwidth=ScaleFactor*0.25;
for i=1:max(Spots(:,10))
    spot_means(i)=mean(Spots(Spots(:,10)==i,5));
    stdx(i)=std(Spots(Spots(:,10)==i,1));
    stdy(i)=std(Spots(Spots(:,10)==i,2));
    traj_color=rand(1,3);
   % subplot(2,3,1)
    plot(Spots(Spots(:,10)==i,1),Spots(Spots(:,10)==i,2),'.','color','r')
   % subplot(2,3,3)
  %  hold on
   % plot(Spots(Spots(:,10)==i,9),Spots(Spots(:,10)==i,5),'.','color',traj_color)
   % title('Trajectory intensities vs frame num')
   t=i;
   Intensity_rec=10000;
   % Intensity_rec=Spots(t,5);
   GaussFrame=GaussFrame+(Intensity_rec./(2.*pi.*PSFwidth.*PSFwidth))*exp(-(((Xpos-ScaleFactor*Spots(t,1)).^2)./(2.*PSFwidth^2)+((Ypos-ScaleFactor*Spots(t,2)).^2)./(2.*PSFwidth^2)));
end
[Idens,Ix] = ksdensity(Spots(Spots(:,10)>0,5));
[Imeansdens,Imeansx] = ksdensity(spot_means);
[Stdxdens,Stdxx] = ksdensity(stdx);
[Stdydens,Stdyx] = ksdensity(stdy);
 X=Spots(:,1:2);
 hold on
%scatter(meanClusterPos(:,1), -(meanClusterPos(:,2)-max(meanClusterPos(:,2))),60,'x','MarkerEdgeColor','k','LineWidth',2)
%scatter(meanClusterPos(:,1), (meanClusterPos(:,2)),60,'x','MarkerEdgeColor','k')
%[I,x]=ksdensity(distancesWithinClusters)
%figure; plot(x,I)
 %[nPositionsInCluster, distancesWithinClusters, distancesBetweenClusters meanClusterPos] = nearestNeighbourClustering(X, 0.2);
%[IDX,C]=kmeans(X,170);
%scatter(C(:,1),C(:,2),60,'x','MarkerEdgeColor','k','LineWidth',2)
%subplot(2,3,4)
% plot(Ix,Idens,'r')
% hold on
% plot(Imeansx, Imeansdens,'b')
% legend('Intensity KDF','Trajectory Mean Intensity KDF')
% subplot(2,3,5)
% plot(Stdxx, Stdxdens,'r')
% hold on
% plot(Stdyx, Stdydens,'b')
% legend('SDx KDF','SDy KDF')
% figure;
% plot(Stdxx, Stdxdens,'r')
% hold on
% plot(Stdyx, Stdydens,'b')
% legend('SDx KDF','SDy KDF')
mean((stdx.^2+stdy.^2).^0.5)

figure;
imshow(GaussFrame,[]);