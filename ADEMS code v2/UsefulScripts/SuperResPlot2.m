% Reads in spots array produced from ADEMS code and plots things, 
% This is a modified version of AJMW's SuperResPlot, this versions plots two lots of tracks
% from the same data e.g. blink on and blink off
%
% Requires user to select frame average( or image you want to plot on)as a .tiff and
% a .mat file of the localisations at a user prompt
% 
% Some functionality from SuperResPlot has just been commented out 


[TifFile TifDir]=uigetfile('*.tif','Select the frame average');
image_data=imread(strcat(TifDir,TifFile));
[DataFile DataDir]=uigetfile(strcat(TifDir,'*.mat'),'Select the blinkon tracking data');
 load(strcat(DataDir,DataFile),'SpotsOn')%,'SpotsCh2');
 Spots1=SpotsOn;
 [DataFile DataDir]=uigetfile(strcat(TifDir,'*.mat'),'Select the blinkoff tracking data');
 load(strcat(DataDir,DataFile),'SpotsOff')
 Spots2=SpotsOff;
 %[DataFile DataDir]=uigetfile(strcat(TifDir,'*.mat'),'Select the raw tracking data');
 %load(strcat(DataDir,DataFile),'ans')
 %Spots3=ans;
%image_data=imread(strcat(DataDir,'frame_average.tif'));
figure;
%image_data=frame_average;


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
plot(Spots1(Spots1(:,10)==0,1),Spots1(Spots1(:,10)==0,2),'.','color','r')
plot(Spots2(Spots2(:,10)==0,1),Spots2(Spots2(:,10)==0,2),'.','color','b')
%plot(Spots3(Spots3(:,10)==0,1),Spots3(Spots3(:,10)==0,2),'.','color','m')
% This seems to be a way of setting the size of reconstruction
ScaleFactor=4; 
GaussFrame=zeros(ScaleFactor*size(image_data));
[Xpos,Ypos] = meshgrid(1:ScaleFactor*size(image_data,2),1:ScaleFactor*size(image_data,1));
PSFwidth=ScaleFactor*0.25;
for i=1:max(Spots1(:,10))
    spot_means1(i)=mean(Spots1(Spots1(:,10)==i,5));
    stdx1(i)=std(Spots1(Spots1(:,10)==i,1));
    stdy1(i)=std(Spots1(Spots1(:,10)==i,2));
    %traj_color=rand(1,3);
   % subplot(2,3,1)
    plot(Spots1(Spots1(:,10)==i,1),Spots1(Spots1(:,10)==i,2),'.','color','c')
   % subplot(2,3,3)
  %  hold on
   % plot(Spots(Spots(:,10)==i,9),Spots(Spots(:,10)==i,5),'.','color',traj_color)
   % title('Trajectory intensities vs frame num')
   t=i;
   Intensity_rec=10000;
   % Intensity_rec=Spots(t,5);
   GaussFrame=GaussFrame+(Intensity_rec./(2.*pi.*PSFwidth.*PSFwidth))*exp(-(((Xpos-ScaleFactor*Spots1(t,1)).^2)./(2.*PSFwidth^2)+((Ypos-ScaleFactor*Spots1(t,2)).^2)./(2.*PSFwidth^2)));
end
for i=1:max(Spots2(:,10))
    spot_means2(i)=mean(Spots2(Spots2(:,10)==i,5));
    stdx2(i)=std(Spots2(Spots2(:,10)==i,1));
    stdy2(i)=std(Spots2(Spots2(:,10)==i,2));
    traj_color=rand(1,3);
   % subplot(2,3,1)
    plot(Spots2(Spots2(:,10)==i,1),Spots2(Spots2(:,10)==i,2),'.','color','y')
   % subplot(2,3,3)
  %  hold on
   % plot(Spots(Spots(:,10)==i,9),Spots(Spots(:,10)==i,5),'.','color',traj_color)
   % title('Trajectory intensities vs frame num')
   t=i;
   Intensity_rec=10000;
   % Intensity_rec=Spots(t,5);
   GaussFrame=GaussFrame+(Intensity_rec./(2.*pi.*PSFwidth.*PSFwidth))*exp(-(((Xpos-ScaleFactor*Spots2(t,1)).^2)./(2.*PSFwidth^2)+((Ypos-ScaleFactor*Spots2(t,2)).^2)./(2.*PSFwidth^2)));
end
%for i=1:max(Spots3(:,10))
 %   spot_means3(i)=mean(Spots3(Spots3(:,10)==i,5));
 %   stdx3(i)=std(Spots3(Spots3(:,10)==i,1));
 %   stdy3(i)=std(Spots3(Spots3(:,10)==i,2));
    %traj_color=rand(1,3);
   % subplot(2,3,1)
  %  plot(Spots3(Spots3(:,10)==i,1),Spots3(Spots3(:,10)==i,2),'.','color','w')
   % subplot(2,3,3)
  %  hold on
   % plot(Spots(Spots(:,10)==i,9),Spots(Spots(:,10)==i,5),'.','color',traj_color)
   % title('Trajectory intensities vs frame num')
 %  t=i;
 %  Intensity_rec=10000;
   % Intensity_rec=Spots(t,5);
 %  GaussFrame=GaussFrame+(Intensity_rec./(2.*pi.*PSFwidth.*PSFwidth))*exp(-(((Xpos-ScaleFactor*Spots3(t,1)).^2)./(2.*PSFwidth^2)+((Ypos-ScaleFactor*Spots3(t,2)).^2)./(2.*PSFwidth^2)));
%end
%%HM - - THIS GRAPH MEANS NOTHING ATM
%%[Idens,Ix] = ksdensity(Spots(Spots(:,10)>0,5));
%%[Imeansdens,Imeansx] = ksdensity(spot_means);
%%[Stdxdens,Stdxx] = ksdensity(stdx);
%%[Stdydens,Stdyx] = ksdensity(stdy);
%% X=Spots(:,1:2);
%% hold on
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
mean1=mean((stdx1.^2+stdy1.^2).^0.5);
mean2=mean((stdx2.^2+stdy2.^2).^0.5);
%mean3=mean((stdx3.^2+stdy3.^2).^0.5);
%((length(Spots1)*mean1)+(length(Spots2)*mean2)+(length(Spots3)*mean3))/(length(Spots1)+length(Spots2)+length(Spots3))
((length(Spots1)*mean1)+(length(Spots2)*mean2))/(length(Spots1)+length(Spots2))
figure;
imshow(GaussFrame,[]);