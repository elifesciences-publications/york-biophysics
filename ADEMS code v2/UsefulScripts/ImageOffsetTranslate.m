function [xoffset, sigmax, yoffset, sigmay]=ImageOffsetTranslate(image_data_left,image_data_right,all,noFrames,startFrame)
%%Reads in brightfield image file, with assumption of same sample on left and right channels on eahc image
% and uses peak in normalized autocorerelation between frame average of two channels to find
% translation from left to right channel

%%earlier version had testing with mean fitlers etc, which is now stripped
%%out. Also stripped out frame averagiong in preference for looking at each
%%frame indivudally - if frame averaging is desired (not ideal anyway as
%%too much motion will blur images and reduce correlations) then this
%%should be better done separately before calling this function


%%INITIALISE PARAMETERS

%clear all
close all


OffsetPixelNum=13; %OffsetPixelNum is number of pixels in x and y to crop image channels on left border, central line between channels and right border, 

%to make sure we do not include either central black region between
%channels or the diffraction fringes in x, but also some smaller dark edges
%in y
subarray_halfwidth = 8; % (Default: 8 pixels). Halfwidth of image square subarray
% ROI, typically a square of size 17x17 pixels. is half (edge length of square -1, eg FitROIhalfwidth=8 gives edge length of 17) centred on initial 
%max position for normalized autocorrelation, then used for refined peak guess with 2D Guassian fit
inner_circle_radius = 5; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray.
gauss_mask_sigma = 2; % (Default: 2 pixels). Size in pixels of the applied Gaussian mask.
guess_sigma_Fit = 3; % starting guess for Gaussian fit of brightspot intensity (default = 3).
sigmaFit_min = 2;
sigmaFit_max = inner_circle_radius; % maximum acceptable sigma of gaussian fit to spot

Testfunction=0;
show_output=0; %if show_output=1, figures are displayed

if Testfunction==1 %if Testfunction=1 then go through testing subroutines...
    all=1; %Use this keyword to load entire image file, loads whole image if all=1
    CameraType=0; %CameraType=0 is Adam 128x128 camera, CameraType=1 is Erik 512x512 camera, but potentially with y-cropped images - this is for testing purposes only (I hope...)
    Test_with_same_image=1; %if Test_with_same_image=1 then simply sets the left channel image to be a copy of the right channel image
    if CameraType==0 %crop image to look just at the relevant channel
    image_label='2014-05-13T11h21m57smig1nrd134' %test brightfield for nrd1-mcherry yeast cells from Adam's setup
    startFrame=1;
    noFrames=10; % number of  frames to read in from file(default: 5, end_frame-start_frame)
    elseif CameraType==1%set this for one of Erik camera brightfield/DIC etc images
    image_label='BF_1_MMStack_Erik_07-05-14_LG_cell4.ome' %test mig1-gfp:nrd1-mcherry fluroescence data from Adam setup, ?561nm only laser (pre alex) startframe 24
    startFrame=1; 
    noFrames=50; % number of  frames to read in from file(default: 5, end_frame-start_frame)
    end
end
endFrame=startFrame+noFrames-1; 

%% OPEN DATA

%Open tif data from image_label which is whatever string is before the file
%extentions
% [numFrames, frame_Ysize, frame_Xsize, image_data_wholeframe, image_path] = ExtractImageSequence(image_label, all, startFrame, endFrame);
% %[numFrames, frame_Ysize, frame_Xsize, image_data, image_path] = extract_image_sequence_dataAWarray(image_label, all);
% disp('data loaded')
% 
% sizewholeimage=size(image_data_wholeframe) %assume colour channels are roughly left and right halves of camera image
% image_data_left=image_data_wholeframe(OffsetPixelNum+1:sizewholeimage(1)-OffsetPixelNum-1,OffsetPixelNum+1:sizewholeimage(2)/2-OffsetPixelNum-1,:);
% image_data_right=image_data_wholeframe(OffsetPixelNum+1:sizewholeimage(1)-OffsetPixelNum-1,sizewholeimage(2)/2+OffsetPixelNum+1:sizewholeimage(2)-OffsetPixelNum-1,:); 


%% cycle between consecutive frames
if noFrames>1 %this is for situations of typically single frame "sequences"
  %  countmax=sizewholeimage(3);
    countmax=startFrame+noFrames;
else
    countmax=1;
end

for i=1:countmax %number of frames in total in file
    clear frame_data_left frame_data_right Normcorr
    frame_data_left =image_data_left(:,:,i);
    frame_data_right =image_data_right(:,:,i);
    frame_L = mat2gray(frame_data_left);
    frame_R = mat2gray(frame_data_right);    

    I_eq_L=frame_L; %retain if do not use contrast enhancement
    I_eq_R=frame_R;
    
if Testfunction==1
    if Test_with_same_image==1
    I_eq_L=frame_R; %set left and right channels the same as a control test
    end
end

%% Pad out the right channel image on each edge by appropate size of the image channel itself (gues this is required when constructing correlation function with other image...)
SizeImageLeft=size(I_eq_L);
SizeImageRight=size(I_eq_R);

if max(SizeImageRight<SizeImageLeft==1)
    Autocorrelpad=max(abs(SizeImageLeft-SizeImageRight));
    ImagePadded_Right=zeros(SizeImageRight(1)+2*Autocorrelpad,SizeImageRight(2)+2*Autocorrelpad);
    ImagePadded_Right(Autocorrelpad+1:SizeImageRight(1)+Autocorrelpad,Autocorrelpad+1:SizeImageRight(2)+Autocorrelpad)=I_eq_R; %not sure this required as Autocorrelpad=0 seems to perform ok
else
    ImagePadded_Right=I_eq_R;
end
    
%% normilzed autoceorrelation between left channel and padded right channel

Normcorr = normxcorr2(I_eq_L,ImagePadded_Right);

[max_Normcorr, imax_Normcorr] = max(abs(Normcorr(:)));
[ypeak_init, xpeak_init] = ind2sub(size(Normcorr),imax_Normcorr(1));



% Iterative gaussian masking to determine spot centre
% [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
% findSpotCentre2_ML_v1(Normcorr,xpeak_init,ypeak_init,subarray_halfwidth,inner_circle_radius,gauss_mask_sigma, 0,0);
[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
findSpotCentre2(Normcorr,xpeak_init,ypeak_init,subarray_halfwidth,inner_circle_radius,gauss_mask_sigma,0.05, 0,0);
%Create variable for interpolation later               
% If a spot was found and it didn't clip
if noConvergenceFlag==0 && clipping_flag==0% Fit a 2D Gaussian to it
    [sdx, sdy, Icent] = fit2DgaussianFixedCenter(Normcorr,Ibg_avg, Isp, x_centre, y_centre,subarray_halfwidth, ...
                        guess_sigma_Fit, sigmaFit_max, sigmaFit_min, 0);
else
    disp('no convergence to Gaussian masking')
end

SizeNormcorr=size(Normcorr);
SizeI_eq_L=size(I_eq_L);
%xoffset_i(i)=SizeNormcorr(2)/2-x_centre+sizewholeimage(2)/2+0.5; %not sure I fully understand why adding 0.5 here is correct, but this works when doing same image in both channels as control
xoffset_i(i)=SizeNormcorr(2)/2-x_centre+0.5; %not sure I fully understand why adding 0.5 here is correct, but this works when doing same image in both channels as control
sigmax_i(i)=sdx;
yoffset_i(i)=SizeNormcorr(1)/2-y_centre+0.5;
sigmay_i(i)=sdy;

if show_output==1
    figure(1)
    subplot(2,1,1)
    imshow(I_eq_L,[])
    title('raw frame average left channel')
    subplot(2,1,2)
    imshow(I_eq_R,[])
    title('raw frame average right channel')
    figure(2)
    surf(Normcorr)
    shading flat
    
    title('normalized autocorrelation between template left channel and padded right channel')
% create ROI centred on initial autocorrelation peak estimate, and perform 2D Gaussian fit to refine estimate
%Ipeak=Normcorr(ypeak_init-subarray_halfwidth:ypeak_init+subarray_halfwidth,xpeak_init-subarray_halfwidth:xpeak_init+subarray_halfwidth);
%figure(3)
%surf(Ipeak)
%shading flat
%title('Peak of normalized autocorrelation between template left channel and padded right channel')
end
end
xoffset=mean(xoffset_i);
yoffset=mean(yoffset_i);
sigmax=std(xoffset_i);
sigmay=std(xoffset_i);




