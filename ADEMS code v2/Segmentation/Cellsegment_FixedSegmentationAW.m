function [Mask, frame0,threshold_cell_auto] = Cellsegment_FixedSegmentationAW(image_label,AnalyseChannel,noFrames,threshold_cell,ALEX,AutoThresh,weight_factor)

%% segments images on basis of simple fixed threshold, either manual or autmated using peak and sd to pixle intensity distribution
%weight_factor=1 by default but can be used to shift thresholding a bit
%(multiple of fwhm)
all=0; %loads whole image if all=1
startFrame=1;
Mask=[];
image_label
if exist('weight_factor')==0
    weight_factor=1;
end

%% OPEN DATA

%Open tif data from image_label
if ALEX==0
    endFrame=startFrame+75;
    try
    [numFrames, frame_Ysize, frame_Xsize, image_data_wholeframe, ~] = ExtractImageSequence(image_label, all, startFrame, endFrame);
    catch
        endFrame=49;
        [numFrames, frame_Ysize, frame_Xsize, image_data_wholeframe, ~] = ExtractImageSequence(image_label, all, startFrame, endFrame);
    end
    firstLeft = LaserOn2(image_data_wholeframe, 1,0)
    if firstLeft+noFrames>endFrame
        firstLeft=endFrame-noFrames
    end
    sizewholeimage=size(image_data_wholeframe);
    image_data_left=image_data_wholeframe(:,1:sizewholeimage(2)/2,:);
    image_data_right=image_data_wholeframe(:,sizewholeimage(2)/2:sizewholeimage(2),:);
     if AnalyseChannel==0
         image_data=image_data_left;
     elseif AnalyseChannel==1
         image_data=image_data_right;
     end

else
    endFrame=startFrame+2*noFrames-1;

     [numFrames, frame_Ysize, frame_Xsize, image_data_wholeframe, image_path] = ExtractImageSequence(image_label, all, startFrame, endFrame);
    image_data_odd_frames=image_data_wholeframe(:,:,1:2:end);   % Extract all the odd elements
    image_data_even_frames=image_data_wholeframe(:,:,2:2:end);   % Extract all the even elements
    sizewholeimage=size(image_data_odd_frames);
    
    image_data_odd_frames_left=image_data_odd_frames(:,1:sizewholeimage(2)/2,:);
    image_data_odd_frames_right=image_data_odd_frames(:,sizewholeimage(2)/2:sizewholeimage(2),:);
    image_data_even_frames_left=image_data_even_frames(:,1:sizewholeimage(2)/2,:);
    image_data_even_frames_right=image_data_even_frames(:,sizewholeimage(2)/2:sizewholeimage(2),:);

    % Adam, I think ytou write something like this yourself, but I couldn't
    % find where so just wrote following below to work out start ALEX
    % frames for left or right channels etc...
    if AnalyseChannel==0 & (max(image_data_odd_frames_left(:))>max(image_data_even_frames_left(:)))%crop image to look just at the relevant channel
        image_data=image_data_odd_frames_left;
    elseif AnalyseChannel==0 & (max(image_data_odd_frames_left(:))<=max(image_data_even_frames_left(:)))
        image_data=image_data_even_frames_left;
    elseif AnalyseChannel==1 & (max(image_data_odd_frames_right(:))>max(image_data_even_frames_right(:)))
        image_data=image_data_odd_frames_right;
    elseif AnalyseChannel==1 & (max(image_data_odd_frames_right(:))<=max(image_data_even_frames_right(:)))
        image_data=image_data_even_frames_right;
    end
end
    
%% CALCULATE FRAME AVERAGE AND THEN DETERMINE CELL BOUNDARY FROM THIS

%Calculate Frame Average of the data
frame0 = FrameAverage(image_data, noFrames, firstLeft);
I_eq = mat2gray(frame0);
[histy,histx]=imhist(frame0);


    if AutoThresh==1
        [widthx, maxminvalue] = fwhm(histx,histy);
        threshold_cell_auto=maxminvalue+weight_factor*widthx %this seems to kinda work for Adam Mig1-GFP data for cell outlines
    elseif AutoThresh==2
        [widthx, maxminvalue] = fwhm(histx,log(histy)); %this seems to kinda work for Erik data for cell outlines, and for Adam data for nuclei
        threshold_cell_auto=maxminvalue+widthx;
        if threshold_cell_auto>1
            threshold_cell_auto=threshold_cell;
        end
        
            
    else
        threshold_cell_auto=threshold_cell; %threshold_cell of 0.10-0.12 seems reasonably robust for Adam data for edge of cell with Mig1-GFP images - 
%noFrames higher (eg 10) will ? find nucleus. Will find nucleus for
%Nrd1-mCherry data
    end

bw = im2bw(I_eq, threshold_cell_auto); %create BLOBs based on fixed threshold value, 0.1069 based on previous Gaussian fit to bg pixels in Adam data... not really happy that this is simply fixed though
SizeFrame=size(I_eq);
     
    %clean that up and then overlay the perimeter on the original image.
    bw2 = imfill(bw,'holes'); %fills in holes in image by flooding
    %se=ones(4,4); %square se
    se = strel('disk',4);%disk se
    bw3 = imopen(bw2, se );%morphological opening to remove objects smaller than se 
    bw4 = bwareaopen(bw3, 50); %Remove small objects from binary image here with <50pixels total (inner ROI has 5 pixel raidus, so take 4 whole pixels radius as limiting size)
   
cc = bwconncomp(bw4,4);
L_cells = labelmatrix(cc);
 
    
   % clear Mask
    for lbl = 1:max(L_cells(:)); %cycle through all objects in labe matrix... asumume only one nucleus per cell... could fall down here if problem with extended maxima of course!   

            Mask(:,:,lbl) = L_cells == lbl; %# find pixels belonging to current label. First object is the borders I think, so take the objects after these as indivudal cell masks  
    end
    Mask=Mask;
end
     

