function [CellObject, NucleusObject, frame0,threshold_cell_auto,threshold_nucleus_auto] = Cellsegment_CellAndNucleusAW(image_label,AnalyseChannel,noFrames,threshold_cell,threshold_nucleus,ALEX,AutoThresh)

%% segments images on basis peak and sd to background intensity pixels

%previous version had various segmentation options (Otzu, optimal, etc),
%which are now stripped out, and previous version had a workign ellipse
%fitting to the raw masks, but this better implemented outside of the
%segmentation function I think, if required.

% Open tif data, calculate frame average, determines
%cell perimeter using fixed threshold of peak intensity + n*sd noise (n set in function as 1.4), 
%extended maxima to get very bright objects (?centre of nuclei in cells),
%then watershed transform to link v bright objects and get final cell
%perimeter(s) and mask(s). 
%optional output for fitted ellipse or raw masks
%CellObject, NucleusObject are binary masks array stack of n+1 images each
%same dimension as 1/2 of image_label image selected by left
%(AnalyseImage=0) or right (=1), where n is number of cell objetcs found
%(assumes 1 bright object at centre of each which is nucleus); disregard
%first in the output stack, frame0 is just the frame average calcuated from raw
%image

%Assumes no ALEX currently, so best to implement ALEX separately?

%noFrames=5; % number of frames to average over, starting from start_frame. (default: 5, end_frame-start_frame)

all=0; %loads whole image if all=1
ncellwidth=1; %this is a fix!
nnucleuswidth=11;
startFrame=1;
%AnalyseChannel=1; %AnalyseChannel=0 is left channel only, AnalyseChannel=1 is right channel only


%show_output=0;
show_output=0;
CellObject=[];
NucleusObject=[];
CellObject_perim=[];
NucleusObject_perim=[];

%% OPEN DATA

%Open tif data from image_label
if ALEX==0
    endFrame=startFrame+75;
    [numFrames, frame_Ysize, frame_Xsize, image_data_wholeframe, image_path] = ExtractImageSequence(image_label, all, startFrame, endFrame);
    
    
  [firstLeft, firstRight, LeftAverage, RightAverage] = LaserOn2(image_data_wholeframe, 1,0);
    
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
frame0 = Frame_Average(image_data, noFrames, firstLeft);
I_eq = mat2gray(frame0);
[histy,histx]=imhist(frame0);


    switch AutoThresh
        case 1
        [widthx, maxminvalue] = fwhm(histx,histy);      
        threshold_cell_auto=maxminvalue+ncellwidth*widthx
        threshold_nucleus_auto= (1+threshold_cell_auto)/2 

        case 2
         [widthx, maxminvalue] = fwhm(histx,log(histy));        
        threshold_cell_auto=maxminvalue+ncellwidth*widthx
        threshold_nucleus_auto= (1+threshold_cell_auto)/2 
        otherwise
        threshold_cell_auto=threshold_cell %threshold_cell of 0.10-0.12 seems reasonably robust for Adam data for edge of cell with Mig1-GFP images - 
        threshold_nucleus_auto=threshold_nucleus %threshold_nucleus of ~0.5-0.7 pulls out nuclei from Adam GFP data for HG with noFrames=5, not for LG but if set 
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
    bw4_perim = bwperim(bw4); %get perimeter of all remaning BLOBS
    bw4_perim_filled=imfill(bw4_perim,'holes'); %fill in perimeter to generate mask
    %overlay1_perim = imoverlay(frame0, bw4_perim, [.3 1 .3]); %imoverlay is a new function, overlays here perimeter with frame average
    overlay1_perim = imoverlay(frame0, bw4_perim, [1 1 0]); %imoverlay is a new function, overlays here perimeter with frame average
 %from this line below is marker-based watershed segmentation method based on bright
 
   
 %objects in image
    mask_em = imextendedmax(I_eq, threshold_nucleus_auto); %extended maxima transform to identify v bright objects, not sure yet about setting threshold here
    %mask_em = imextendedmax(I_eq, 0.6);%mask_em = imextendedmax(I_eq, threshold_nucleus); %extended maxima transform to identify v bright objects, not sure yet about setting threshold here
    %0.6 works ok for fluorescennce to identify centre of nuclei
    mask_em_perim=bwperim(mask_em);
    cc = bwconncomp(mask_em,4);
number_nuclei  = cc.NumObjects;
L_nuclei = labelmatrix(cc);
   
     overlay1_perim_filled = imoverlay(frame0, bw4_perim_filled, [1 1 0]); %imoverlay is a new function, overlays here perimeter with frame average
    % overlay1_perim_filled is a plausiable mask provided number opf frames
    % averaged over is high enough to allow all particles to explore
    % relevant compartment of cell
    
    
    
    %clean that up and then overlay it
    mask_em = imclose(mask_em, ones(5,5)); %morphological closing - fills in gaps
    mask_em = imfill(mask_em, 'holes'); %fills in any holes
    mask_em = bwareaopen(mask_em, 40); % removes small objects, here <40 pixels in total
    %overlay2 = imoverlay(I_eq, bw4_perim | mask_em, [.3 1 .3]); %overlay perimeter and v bright objects 
    overlay2 = imoverlay(overlay1_perim, mask_em_perim, [0 1 1]); %overlay perimeter and v bright objects 
    
    
    
    %watershed bit... 
    I_eq_c = imcomplement(I_eq); %complement the image because we are about to apply the watershed transform, which identifies low points, not high points.
    I_mod = imimposemin(I_eq_c, ~bw4 | mask_em); %modify image so background and extended maxima pixels forced to be only local minima in image.
    L_cells = watershed(I_mod); %compute watershed transform... point of this I think to segment between cells touching eachother.. but does strange stuff when just one cell present soemtimes
    %Lperim=L; %perimeter labelmatrix
    
 %   clear CellObject NucleusObject CellObject_perim NucleusObject_perim
    cell_counter=1;
    nucleus_counter=1;
    for lbl = 1:max(L_cells(:)); %cycle through all cell objects in labe matrix...  

        if lbl>1; %generate perimeters of all objects found (mainly for visual display I think... assume same number of nuclei as cells./.. bound to go wrong here!
            CellObject_init = L_cells == lbl; %# find pixels belonging to current label. First object is the borders I think, so take the objects after these as indivudal cell masks
            NucleusObject_init = L_nuclei == lbl-1; %find pixels belong to equvalent label in nuclues labelmatrix L_nuclei... nb here the first labeled object is fine
         % For first "object" in label matrix set this as perimeters of all of the objects together 
        s_cell=regionprops(CellObject_init,'Area','EquivDiameter');
        s_nucleus=regionprops(NucleusObject_init,'Area','EquivDiameter');
        cell_area=s_cell.Area;
        cell_diameter=s_cell.EquivDiameter;
        nucleus_area=s_nucleus.Area;
        nucleus_diameter=s_nucleus.EquivDiameter;
        if cell_diameter>10
            CellObject(:,:,cell_counter)=CellObject_init;
            cell_counter=cell_counter+1;
        end
        if nucleus_diameter>5
            NucleusObject(:,:,nucleus_counter)=NucleusObject_init;
            nucleus_counter=nucleus_counter+1;
        end

            
            

        
        end
   

        
    end
    
 

    
    
if show_output==1
figure(1)
    subplot(2,1,1)
    imshow(I_eq,[])
    title('frame average')
    subplot(2,1,2)
    imshow(bw,[])
    title('pixel-thresholded')
    
    figure(2)
    subplot(2,1,1)
    imshow(overlay1_perim)
    title('overlaid perimeters')
    subplot(2,1,2)
    imshow(overlay1_perim_filled)
    title('overlaid filled perimeters')
    
    figure(3)
    subplot(2,1,1)
    imshow(mask_em);
    title('extended maxima transform')
    subplot(2,1,2)
    imshow(overlay2)
    title('permiter and v bright extended maxima objects')

    figure(4)
    subplot(2,1,1)
    imshow(label2rgb(L_cells))
    title('watershed segmentation between v bright objects and Otzu perimeters')
    subplot(2,1,2)
    imhist(frame0)
    title('pixel intensity distribution frame average image- this can be quite revealing...')

    figure(6)
    imshow(label2rgb(L_cells))
    title('warer-shed output for raw cell boundaries')
end




