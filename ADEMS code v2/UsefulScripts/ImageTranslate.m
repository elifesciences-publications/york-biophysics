function [yoffset_final, xoffset_final, yoffset_sd, xoffset_sd]=ImageTranslate(BFImage,BFImage2,show_output)
BFImage;
BFImage2;
exclude_region=3;
subarray_halfwidth = 8; % (Default: 8 pixels). Halfwidth of image square subarray
% ROI, typically a square of size 17x17 pixels. is half (edge length of square -1, eg FitROIhalfwidth=8 gives edge length of 17) centred on initial
%max position for normalized autocorrelation, then used for refined peak guess with 2D Guassian fit
inner_circle_radius = 5; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray.
gauss_mask_sigma = 2; % (Default: 2 pixels). Size in pixels of the applied Gaussian mask.
guess_sigma_Fit = 3; % starting guess for Gaussian fit of brightspot intensity (default = 3).
sigmaFit_min = 2;
sigmaFit_max = inner_circle_radius; % maximum acceptable sigma of gaussian fit to spot
%show_output=1;


for i=1:size(BFImage,3)
    
    % Compare only right channels as better contrast
    
    Frame1=mat2gray(BFImage(25:end-25,size(BFImage2,2)/2+3*exclude_region+10:end-exclude_region-10,i));
    %Frame2=mat2gray(BFImage2(exclude_region:end-exclude_region,size(BFImage2,2)/2+3*exclude_region:end-exclude_region,i));
    Frame2=mat2gray(BFImage2(1:end,size(BFImage2,2)/2:end,i));
    if show_output==1
        figure;
        imshow(Frame1)
        figure;
        imshow(Frame2)
    end
    corr_image=normxcorr2(Frame1,Frame2);
    [ypeak, xpeak] = find(corr_image==max(corr_image(:)));
    [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
        findSpotCentre2(corr_image,xpeak,ypeak,subarray_halfwidth,inner_circle_radius,gauss_mask_sigma,0.05, 0,0);
    yoffSet(i) = y_centre-size(Frame1,1)-0.5*(size(Frame2,1)-size(Frame1,1));
    %  yoffSet = ypeak+1.5*size(Frame1,1)-0.5*size(Frame2,1);
    xoffSet(i) = x_centre-size(Frame1,2)-0.5*(size(Frame2,2)-size(Frame1,2));
    % xoffSet = xpeak+1.5*size(Frame1,2)-0.5*size(Frame2,2);
    if show_output==1
        hFig = figure;
        hAx  = axes;
        imshow(Frame2,'Parent', hAx);
        imrect(hAx, [xoffSet(i), yoffSet(i), size(Frame2,2), size(Frame2,1)]);
    end
    
end
yoffset_final=mean(yoffSet);
xoffset_final=mean(xoffSet);
yoffset_sd=std(yoffSet);
xoffset_sd=std(xoffSet);
end