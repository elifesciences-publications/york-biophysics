% Generates a movie of the calculated trajectories on the image data and
% uses matlabs built in player

[TifFile TifDir]=uigetfile('*.tif','Select the tiff stack');
cd(TifDir)
[DataFile DataDir]=uigetfile('*.mat','Select the correct tracking data');
pause on
load(strcat(DataDir,DataFile),'SpotsCh1','SpotsCh2');
SpotsCh1=sortrows(SpotsCh1,9);
cd(TifDir)
[numFrames, frame_Ysize, frame_Xsize, image_data, image_path] = extract_image_sequence_dataAWarray(TifFile(1:end-4), 0, min(SpotsCh1(:,9)), max(SpotsCh1(:,9)));
TrajColour=rand(max(SpotsCh1(:,10))+1,3);
TrajColour(:,1)=0;
TrajColour(1,:)=[1,0,0];
loopcounter=0;
for i=min(SpotsCh1(:,9)):max(SpotsCh1(:,9))
    
    loopcounter=loopcounter+1;
    TrajInFrame=[];
    plot1=subplot(2,2,1);
    imshow(image_data(:,:,i),[]);
    hold on
    FrameIndex=find(SpotsCh1(:,9)==i);
    % TrajInFrame=SpotsCh1(SpotsCh1(:,9)==i,10);
    if isempty(FrameIndex)==0
        for k=1:FrameIndex(end)
            if SpotsCh1(k,10)>0
                if ismember(SpotsCh1(k,10),SpotsCh1(SpotsCh1(:,9)==i,10))==1
                    
                    plot(SpotsCh1(k,1),SpotsCh1(k,2),'-o','color',TrajColour(SpotsCh1(k,10)+1,:))%,'-o','color',traj_color)
                    
                end
            end
        end
    end
    plot4=subplot(2,2,4);
    hold on
    title('Trajectory intensity vs frame number')
    axis([min(SpotsCh1(SpotsCh1(:,10)>0,9)) max(SpotsCh1(SpotsCh1(:,10)>0,9))...
        min(SpotsCh1(SpotsCh1(:,10)>0,5)) max(SpotsCh1(SpotsCh1(:,10)>0,5))])
    if isempty(FrameIndex)==0
        for k=1:FrameIndex(end)
            if SpotsCh1(k,10)>0
                if ismember(SpotsCh1(k,10),SpotsCh1(SpotsCh1(:,9)==i,10))==1
                    
                    
                    
                    plot(SpotsCh1(k,9),SpotsCh1(k,5),'-o','color',TrajColour(SpotsCh1(k,10)+1,:))
                end
            end
        end
    end
    
    
    
    
    
    title('Trajectories from start to current frame')
    hold off
    plot2=subplot(2,2,2);
    imshow(image_data(:,:,i),[]);
    hold on
    title('spots found in frame, unlinked in red')
    %plot(SpotsCh1(SpotsCh1(:,9)==i,1),SpotsCh1(SpotsCh1(:,9)==i,2),'o','color',TrajColour(SpotsCh1(k,10)+1,:))
    for q=min(FrameIndex):max(FrameIndex)
        subplot(2,2,2)
        plot(SpotsCh1(q,1),SpotsCh1(q,2),'o','color', TrajColour(SpotsCh1(q,10)+1,:))
        text(SpotsCh1(q,1)+3,SpotsCh1(q,2)+3,num2str(SpotsCh1(q,10)),'color',TrajColour(SpotsCh1(q,10)+1,:))
        if SpotsCh1(q,10)>0
            subplot(2,2,1)
            text(SpotsCh1(q,1)+3,SpotsCh1(q,2)+3,num2str(SpotsCh1(q,10)),'color',TrajColour(SpotsCh1(q,10)+1,:))
        end
    end
   plot3= subplot(2,2,3);
    imshow(image_data(:,:,i),[]);

    title('raw image data')
    
    %  pause
    F(loopcounter)=getframe(gcf);
    
    close
end
implay(F)
% figure;
% axes('pos',[0 0 1 1],'visible','off');
% movie(F,3,3)
