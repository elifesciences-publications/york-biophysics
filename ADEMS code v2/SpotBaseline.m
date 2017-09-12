subarray_halfwidth = 8; % (Default: 8 pixels). Halfwidth of image square subarray
% ROI, typically a square of size 17x17 pixels.
inner_circle_radius = 5; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray.
gauss_mask_sigma = 2; % (Default: 2 pixels). Size in pixels of the applied Gaussian mask.
clear spotsBaseline
[TifFile TifDir]=uigetfile('*.tif','Select the tiff stack');
[DataFile DataDir]=uigetfile('*.mat','Select the correct tracking data');
% DataStruct=load(strcat(DataDir,DataFile));
% SpotsCh1= DataStruct.(char(fieldnames(DataStruct)));
load(strcat(DataDir,DataFile),'SpotsCh1','SpotsCh2');
cd(TifDir)
[numFrames, frame_Ysize, frame_Xsize, image_data, image_path] = ...
    extract_image_sequence_dataAWarray(TifFile(1:end-4), 0, min(SpotsCh1(:,9)), max(SpotsCh1(:,9)));
SpotLastFrame=zeros(max(SpotsCh1(:,10)),1);
LFindex=zeros(max(SpotsCh1(:,10)),1);
y_estimate=zeros(max(SpotsCh1(:,10)),1);
x_estimate=zeros(max(SpotsCh1(:,10)),1);
spot_num=size(SpotsCh1,1)+1;
spotsBaseline=SpotsCh1;
%Loop over trajectory numbers
for k=1:max(SpotsCh1(:,10))
    [SpotLastFrame(k), LFindex(k)]=max(SpotsCh1(SpotsCh1(:,10)==k,9));
    SpotsCh1X=SpotsCh1(SpotsCh1(:,10)==k,1);
    SpotsCh1Y=SpotsCh1(SpotsCh1(:,10)==k,2);
    y_estimate(k,1)=SpotsCh1Y(end);
    x_estimate(k,1)=SpotsCh1X(end);
    %Loop over frames
    if SpotLastFrame(k)+10 < max(SpotsCh1(:,9))
    for p=SpotLastFrame(k)+1:SpotLastFrame(k)+11
        frame=image_data(:,:,p);
        [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
            findSpotCentre2noloop(frame,x_estimate(k),y_estimate(k),subarray_halfwidth,inner_circle_radius,gauss_mask_sigma,0.05, 0,0);
        snr1=Isp/(bg_noise_std*mask_pixels);
        spotsBaseline(spot_num,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, inner_circle_radius, inner_circle_radius, -1, p,k, snr1, SpotsCh1(1,12)];
        spot_num=spot_num+1;
    end
    end
end
datafilename=strcat(DataDir,'\','Baseline',DataFile);
save(datafilename,'spotsBaseline');