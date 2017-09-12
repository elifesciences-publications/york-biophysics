% Reads in spots array and plots things, uses frame average in the same
% % folder
[TifFile TifDir]=uigetfile('*.tif','Select the frame average');
image_data=imread(strcat(TifDir,TifFile));
[DataFile DataDir]=uigetfile(strcat(TifDir,'*.mat'),'Select the correct tracking data');
 load(strcat(DataDir,DataFile),'SpotsCh1','SpotsCh2');
%image_data=imread(strcat(DataDir,'frame_average.tif'));
%image_data=frame_average;
figure;
%Spots=SpotsCh1;
Spots=spotsBaseline;
clear spot_means stdx stdy
subplot(2,3,2)
imshow(image_data,[])
title('Frame average with isolated spots superimposed')
hold on
plot(Spots(Spots(:,10)==0,1),Spots(Spots(:,10)==0,2),'o','color','r')
subplot(2,3,1)
imshow(image_data,[])
title('Frame average with tracks superimposed')
hold on
for i=1:max(Spots(:,10))
    spot_means(i)=mean(Spots(Spots(:,10)==i,5));
    stdx(i)=std(Spots(Spots(:,10)==i,1));
    stdy(i)=std(Spots(Spots(:,10)==i,2));
    spot_size(i)=size(Spots(Spots(:,10)==i,5),1);
    traj_color=rand(1,3);
    subplot(2,3,1)
    plot(Spots(Spots(:,10)==i,1),Spots(Spots(:,10)==i,2),'-o','color',traj_color)
    subplot(2,3,3)
    hold on
    plot(Spots(Spots(:,10)==i,9),Spots(Spots(:,10)==i,5),'-o','color',traj_color)
    title('Trajectory intensities vs frame num')
    %pause
end
[Idens,Ix] = ksdensity(Spots(Spots(:,10)>0,5));
[Imeansdens,Imeansx] = ksdensity(spot_means);
[Stdxdens,Stdxx] = ksdensity(stdx);
[Stdydens,Stdyx] = ksdensity(stdy);

subplot(2,3,4)
plot(Ix,Idens,'r')
hold on
plot(Imeansx, Imeansdens,'b')
legend('Intensity KDF','Trajectory Mean Intensity KDF')
subplot(2,3,5)
plot(Stdxx, Stdxdens,'r')
hold on
plot(Stdyx, Stdydens,'b')
legend('SDx KDF','SDy KDF')

