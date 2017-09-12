function [SpotsCh1, SpotsCh2, frame_average] = FinalCode12(image_label)

% Final code for opening tif data, calculating frame average, using this to determine
%cell boundary, then looping over user defined frames, finding spots,
%identifying their centres and accepting them if they meet criteria, then
%linking these spots together into trajectories.

% Number of frames to average over when calculating a frame average in
% used  as image to choose spots in cursor
% mode and output
noFrames=5; % number of frames to average over, starting from start_frame. (default: 5, end_frame-start_frame)
% PARAMETERS for finding spot centres 
% The total integrated spot intensity is a bgnd corrected one, inside a
% circular mask of radius inner_circle_radius.
subarray_halfwidth = 8; % (Default: 8 pixels). Halfwidth of image square subarray
%subarray_halfwidth = 7; % (Default: 8 pixels). Halfwidth of image square subarray
% ROI, typically a square of size 17x17 pixels.
inner_circle_radius = 5; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray.
%inner_circle_radius = 3; % (Default: 5 pixels). Radius of inner circular mask that moves inside the fixed square subarray.
gauss_mask_sigma = 1; % (Default: 2 pixels). Size in pixels of the applied Gaussian mask.
guess_sigma_Fit = 2; % starting guess for Gaussian fit of brightspot intensity (default = 3).

% PARAMETERS for deciding if we accept a spot centre 
sigmaFit_min = 0;  % minimum acceptable sigma of gaussian fit to spot, in pixels (2) (-3).
sigmaFit_max = inner_circle_radius; % maximum acceptable sigma of gaussian fit to spot, in pixels (4) (3).
SNR_min = 0.4; % minimum acceptable signal-to-noise ratio (at least 0.4)

% PARAMETERS for eliminating coincident spots:
d_min= 1; % distance (in pixels) for eliminating coincidences

% PARAMETERS for building trajectories:
% For linking spots in current and previous frames:
d_01_max = 5; % max distance in pixels between spot centres in current and previous frames, for linking them into a trajectory (5).
Iratio_01_min = 0.5; % min ratio of total spot intensities (after bgnd subtraction) (0.5).
Iratio_01_max = 3; % max ratio of total spot intensities (after bgnd subtraction) (frame k-1/frame k) (large enough value (3) to account for blinking).
SigmaRatio_01_min = 0.5; % min ratio of spot widths (sigma of Gaussian fit) (0.5).
SigmaRatio_01_max = 2; % max ratio of spot width (sigma of Gaussian fit) (2).
% Save parameters to results as structure "params":

% Parameter to exclude a region from the middle if splitting CCD
exclude_region = 0; %must be set to zero if using only one channel
%Channels to use
start_channel=1;
end_channel=1;

%Use this keyword to load entire image file
all=1;
% Or specify start and end frames
startFrame=1;
endFrame=999;

% Switch, if ALEX experiment=1
ALEX=0;

% Switch, =1 to determine laser on time with differential rather than max
% intensity
use_diff=1;

% for finding spots in image
disk_radius = 5;

% CSplit defines how the channels are split, =0 for whole frame (no split),
% 1 for left/right and 2 for up/down
CSplit=0;

%Initialise Spot variables
SpotsCh1=[];
SpotsCh2=[];

% There are 2 methods for finding candidate spots which work better for
% different datasets or if set =3 runs both and uses all spots
% (recommended)
CandidateFindMethod=2;


%gaussian=1 if running a gaussian filter over image data before finding
%spots
gaussian=0;
% Use this to specify interesting spots with the cursor, rather than
% autodetecting them
use_cursor=0;
%If this =1, then graphs will appear
show_output=0;
show_text_output=0;
if show_output==1
    pause on
end
%% CREATE FITTYPE AND OPTIONS

% Create fit type for constrained 2D Gaussian fit to spots
myfit = fittype('Ibg_avg+(Isp./(2.*pi.*sdx.*sdy))*exp(-(((x-x0).^2)./(2.*sdx^2)+((y-y0).^2)./(2.*sdy^2)))',...
    'problem', {'Ibg_avg','Isp','x0','y0'}, 'independent', {'x', 'y'}, 'dependent', 'z');
% Fit options:
options = fitoptions(myfit);

% Starting parameter values:
% Icent_start = 2000;
sdx_start = guess_sigma_Fit;
sdy_start = guess_sigma_Fit;
% x0_start = x_centre;
% y0_start = y_centre;

%options.StartPoint = [Icent_start, sdx_start, sdy_start x0_start, y0_start];
options.StartPoint = [sdx_start, sdy_start];

% Lower parameter values:
% Icent_lower = 0;
sdx_lower = sigmaFit_min;
sdy_lower = sigmaFit_min;
% x0_lower = x_centre - disk_radius;
% y0_lower = y_centre - disk_radius;

%options.Lower = [Icent_lower, sdx_lower, sdy_lower x0_lower, y0_lower];
options.Lower = [sdx_lower, sdy_lower];

% Upper parameter values:plot(
% Icent_upper = inf;
sdx_upper = sigmaFit_max;
sdy_upper = sigmaFit_max;
% x0_upper = x_centre + disk_radius;
% y0_upper = y_centre + disk_radius;

%options.upper = [Icent_upper,  sdx_upper, sdy_upper, x0_upper, y0_upper];
options.upper = [sdx_upper, sdy_upper];

%% OPEN DATA

%Open tif data from image_label which is whatever string is before the file
%extention
[numFrames, frame_Ysize, frame_Xsize, image_data, image_path] = ExtractImageSequence(image_label, all, startFrame, endFrame);
%[numFrames, frame_Ysize, frame_Xsize, image_data, image_path] = extract_image_sequence_dataAWarray(image_label, all);
disp('data loaded')

%% DETERMINE LASER ON FRAME

%Determine when laser turned on, detects both channels seperately if ALEX
[firstLeft, firstRight, LeftAverage, RightAverage] = LaserOn2(image_data, use_diff,ALEX,300);
disp('start determined')

%Number of frames to track
FramesToTrack=numFrames-firstLeft;
startFrame=firstLeft;

%% CALCULATE FRAME AVERAGE AND THEN DETERMINE CELL BOUNDARY FROM THIS

%Calculate Frame Average of the data
if startFrame<(endFrame-5)
    frame_average = FrameAverage(image_data, noFrames, startFrame);
    disp('average calculated')
else
    disp('WARNING FRAME AVERAGE NOT CALCULATED AS STARTFRAME TOO CLOSE TO ENDFRAME')
    frame_average=image_data(:,:,startFrame);
end





%%
for Ch=start_channel:end_channel
    % Initialise spot array and spot number
        spots=[];
    
    spot_num=1;
    
    if ALEX==0
        startFrame=firstLeft;
        endFrame=firstLeft+FramesToTrack;
         FrameInt=1;
    else
         FrameInt=2;
        if Ch==1
            startFrame=firstLeft;
            endFrame=firstLeft+FramesToTrack;
        else
            startFrame=firstRight;
            endFrame=firstRight+FramesToTrack;
        end
    end
    
    %Loop over frames
    for i=startFrame:FrameInt:endFrame
         if show_text_output==0
        h=waitbar(i/endFrame);
         end
        %% FIND CANDIDATE SPOTS
        % Divide image into two channels, left and right
        if Ch==1
            if show_text_output==1
            disp('Ch1')
            end
            switch CSplit
                case 0                   
                    frame=image_data(:,:,i);
                case 1
                    frame=image_data(:,1:round(size(image_data,2)/2-exclude_region),i);
                case 2
                    %Need to write this if needed
            end
            SpotsCh1=[];
        elseif Ch==2
             if show_text_output==1
            disp('Ch2')
             end
            switch CSplit
                case 0
                    frame=image_data(:,:,i);
                case 1
                    frame=image_data(:,round(size(image_data,2)/2+exclude_region):end,i);
                case 2
                    %Need to write this if needed
            end
            SpotsCh2=[];
        end
        
        % For cursor mode
        if use_cursor==1
            if i==startFrame
                pause on
                % Display frame average to choose spots
                if Ch==1
                   frame_averageCH=frame_average(:,1:round(size(image_data,2)/(end_channel-start_channel+1)-exclude_region));
                   %  frame_averageCH=frame_average;
                elseif Ch==2
                    frame_averageCH=frame_average(:,round(size(image_data,2)/2+exclude_region):end);
                end
                imshow(frame_averageCH,[],'InitialMagnification','fit')%HM modify magnification
                title('click a spot and hold alt key to select multiple spots, push any key when finished')
                datacursormode on
                pause
                dcm_obj = datacursormode(1);
                info_struct = getCursorInfo(dcm_obj);
                %Loop over spots chosen and pull out co-ordinates
                for q=1:size(info_struct,2)
                    Spot_coords=info_struct(q).Position;
                    y_estimate(q,1)=Spot_coords(2);
                    x_estimate(q,1)=Spot_coords(1);
                end
                close all
            end
        else
            % Create matrix of 1s where spots might be
            % Now 3 methods for doing this
            switch CandidateFindMethod
                case 1
            [result] = findSpots2(frame,2,disk_radius,gaussian,0);
                case 2
            [result] = findSpots3(frame,2,disk_radius,gaussian,0);
                case 3
                    [result1] = findSpots2(frame,2,disk_radius,gaussian,0);
                    [result2] = findSpots3(frame,2,disk_radius,gaussian,0);
                    result=result1+result2;
                    result(result>1)=1;
            end
            % Convert those to spot co-ordinates
            [y_estimate, x_estimate]=ind2sub(size(result), find(result));
        end
         if show_text_output==1
        disp('candidates found')
         end
        %Plot the candidate spots
        if show_output==1
            imshow(frame, [],'InitialMagnification','fit')%HM modify magnification
            hold on
            plot(x_estimate, y_estimate, 'o')
            hold off
            title('candidate spots')
            pause
        end
        %% FIT TO SPOTS AND REJECT
        %Loop over all found spots
        spots_temp=zeros(size(x_estimate,1),12);
        parfor j=1:size(x_estimate,1)
            if use_cursor==0
                % Iterative gaussian masking to determine spot centre
                [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
                    findSpotCentre2(frame,x_estimate(j),y_estimate(j),subarray_halfwidth,inner_circle_radius,gauss_mask_sigma,0.05, 0,0);
                %Create variable for interpolation later
                spot_background(j,:,:,:)=[x_centre, y_centre,  Ibg_avg];
                % Calculate spots signal to noise ratio
                snr1=Isp/(bg_noise_std*mask_pixels);
                
                % If a spot was found and it didn't clip
                if noConvergenceFlag==0 && clipping_flag==0
                    % Fit a 2D Gaussian to it

                    %    snr2=Icent/bg_noise_std;
                    % Only store spots with good enough snr
                    if snr1>SNR_min
                       [sdx, sdy, Icent] = fit2DgaussianFixedCenter2(frame,Ibg_avg, Isp, x_centre, y_centre,subarray_halfwidth, ...
                        guess_sigma_Fit, sigmaFit_max, sigmaFit_min, 0,myfit,options);
                        % The spot array, 10th field is trajectory number,
                        % initialised to 0
%                         spots(spot_num,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, sdx, sdy, Icent, i,0, snr1, firstLeft];
%                         spot_num=spot_num+1;
                        spots_temp(j,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, sdx, sdy, Icent, i,0, snr1, firstLeft];
                       % spot_num=spot_num+1;
                    end
                end
                %If in cursor mode
            else
                % Iterative gaussian masking to determine spot centre, will
                % return initial values if clipping_flag=1
                [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels,noConvergenceFlag]= ...
                    findSpotCentre2(frame,x_estimate(j),y_estimate(j),subarray_halfwidth,inner_circle_radius,gauss_mask_sigma, 0.05, 1,0);
                %Create variable for interpolation later
                spot_background(j,:,:,:)=[x_centre, y_centre,  Ibg_avg];
                % Calculate spots signal to noise ratio
                snr1=Isp/(bg_noise_std*mask_pixels);
                % Fit gaussians to all spots regardless and assign
                % trajectory numbers based on order in loop
                [sdx, sdy, Icent] = fit2DgaussianFixedCenter(frame,Ibg_avg, Isp, x_centre, y_centre,subarray_halfwidth, ...
                    guess_sigma_Fit, sigmaFit_max, sigmaFit_min, 0); %Is this meant to be fit2DgaussianFixed Center
                spots_temp(j,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, sdx, sdy, Icent, i,j, snr1, firstLeft];
               % spot_num=spot_num+1;
                
            end
        end

         spots_temp(spots_temp(:,1)==0,:)=[]; %HM, gets rid of empty rows
        
      
         if isempty('spots')==1
             
             spots=spots_temp;
         else
             if isempty('spots_temp')==0
               %  if spots_temp(1)>0
                     
                     spots=cat(1,spots,spots_temp);
                     spot_num=size(spots,1)+1;
               %  end
             end
         end
        if isempty(spots)==0
             if show_text_output==1
            disp('Centres found and fitted')
             end
            
            %% PLOT ALL FOUND SPOTS ON IMAGE
            %Plot all the found spots on image
            
            if show_output==1
                %To assess the gaussian masking, interpolate the background
                %based on the backgrounds found in the mask
                [xi,yi]=meshgrid(size(frame,1),size(frame,2));
                BGinterpolate=zeros(size(frame,1):size(frame,2));
                %If more than one spot
                if size(spot_background,1)>2
                    BGinterpolate=griddata(spot_background(:,1),spot_background(:,2),spot_background(:,3),xi,yi);
                    BGinterpolate(isnan(BGinterpolate))=mean(spot_background(:,3));
                    spot_BG_corrected=double(frame)-BGinterpolate;
                    %Subtract the one background if only one spot
                else
                    spot_BG_corrected=double(frame)-mean(spot_background(:,3));
                end
                subplot(1,2,1)
                imshow(frame,[],'InitialMagnification','fit')%HM modify magnification
                hold on
                title('Found elipses on original image')
                %    plot(spots(spots(:,9)==i,1),spots(spots(:,9)==i,2), 'o')
                %Plot elipses on original image
                for k=min(find(spots(:,9)==i)):max(find(spots(:,9)==i))
                    text(spots(k,1)+3,spots(k,2)+3,num2str(k),'color','b')
                   % rectangle('Position',[spots(k,1)-spots(k,6)/2,spots(k,2)-spots(k,7)/2,spots(k,6),spots(k,7)],'Curvature',[1,1],'EdgeColor','b')
                   rectangle('Position',[spots(k,1)-spots(k,6),spots(k,2)-spots(k,7),spots(k,6)*2,spots(k,7)*2],'Curvature',[1,1],'EdgeColor','b')
                end
                hold off
                %plot the BG corrected image from interpolated BG
                subplot(1,2,2)
                imshow(spot_BG_corrected,[],'InitialMagnification','fit')%HM modify magnification
                hold on
                title('interpolated BG corrected image')
                %plot(spot_background(:,1),spot_background(:,2), 'o')
                hold off
                pause
                close all
            end
            
            %% REMOVE COINCIDENT SPOTS
            %spotsold=spots;
            %Function calculated pairwise distances and merges any closer
            %than dmin
            [spots, new_spots]=MergeCoincidentSpots3(spots, i, d_min);
             spots_temp=zeros(size(new_spots,1),12);
            %If any coincidences were found
            if new_spots>0
                %Loop over all found coincident spots
                parfor p=1:size(new_spots,1)
                    %Re calculate centre
                    [x_centre, y_centre, clipping_flag, Ibg_avg, Isp, Idata, bg_noise_std, mask_pixels, noConvergenceFlag]= ...
                        findSpotCentre2(frame,spots(new_spots(p),1),spots(new_spots(p),2),subarray_halfwidth,inner_circle_radius,gauss_mask_sigma, 0.05,0,0);
                    snr1=Isp/(bg_noise_std*mask_pixels);
                    %If converged
                    if noConvergenceFlag==0
                         if snr1>SNR_min
                        [sdx, sdy, Icent] = fit2DgaussianFixedCenter2(frame,Ibg_avg, Isp, x_centre, y_centre,subarray_halfwidth, ...
                            guess_sigma_Fit, sigmaFit_max, sigmaFit_min, 0,myfit,options);
                        %    snr2=Icent/bg_noise_std;
                        % Only store spots above criteria
                       
                            % Store spot in the array
                            spots_temp(p,:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, sdx, sdy, Icent, i,0, snr1, firstLeft];
                         %   spots(new_spots(p),:)=[x_centre, y_centre, clipping_flag, Ibg_avg, Isp, sdx, sdy, Icent, i,0, snr1, firstLeft];
                        end
                    end
                end
                for p=1:size(new_spots,1)
                      if spots_temp(p,1)>0
                    spots(new_spots(p),:)=spots_temp(p,:);
                      end
                end
                %Remove spot number left over from merge
                spots(spots(:,1)==100000,:)=[];
                %Revalue spot number as spots have been removed
                spot_num=size(spots,1)+1;
                 if show_text_output==1
                disp('Coincident spots removed')
                 end
            end
            
            %% PLOT SPOTS WITHOUT COINCIDENT SPOTS
            if show_output==1
                imshow(frame,[],'InitialMagnification','fit')%HM modify magnification
                hold on
                title('elipses on original image-coincident spots')
                %    plot(spots(spots(:,9)==i,1),spots(spots(:,9)==i,2), 'o')
                for k=min(find(spots(:,9)==i)):max(find(spots(:,9)==i))
                    text(spots(k,1)+3,spots(k,2)+3,num2str(k),'color','b')
                  %  rectangle('Position',[spots(k,1)-spots(k,6)/2,spots(k,2)-spots(k,7)/2,spots(k,6),spots(k,7)],'Curvature',[1,1],'EdgeColor','b')
                     rectangle('Position',[spots(k,1)-spots(k,6),spots(k,2)-spots(k,7),spots(k,6)*2,spots(k,7)*2],'Curvature',[1,1],'EdgeColor','b')
                end
                hold off
                pause
                close all
            end
            
            
            %% LINK SPOTS INTO TRAJECTORIES
            if use_cursor==0
                %Start at 2nd frame
                if i>startFrame
                    [spots]=LinkSpots3(spots, i, i-FrameInt, d_01_max, Iratio_01_min,...
                        Iratio_01_max, SigmaRatio_01_min, SigmaRatio_01_max);
                end
            else
                % Don't need to link as cursor spots already linked
            end
             if show_text_output==1
            disp('trajectories determined')
             end
        else
             if show_text_output==1
            disp('No spots found')
             end
        end
    end
    %% FINAL PLOT
    if isempty(spots)==0
        if show_output==1
            if max(spots(:,10))>0
                subplot(2,3,2)
                imshow(image_data(:,:,startFrame),[],'InitialMagnification','fit')%HM modify magnification
                title('First frame with isolated spots superimposed')
                hold on
                plot(spots(spots(:,10)==0,1),spots(spots(:,10)==0,2),'o','color','r')
                subplot(2,3,1)
                imshow(image_data(:,:,startFrame),[],'InitialMagnification','fit')%HM modify magnification
                title('First frame with tracks superimposed')
                hold on
                for i=1:max(spots(:,10))
                    traj_color=rand(1,3);
                    subplot(2,3,1)
                    title('Trajectory intensities vs frame num')
                    if Ch==1
                        plot(spots(spots(:,10)==i,1),spots(spots(:,10)==i,2),'-o','color',traj_color)
                    elseif Ch==2
                        plot(spots(spots(:,10)==i,1)+round(size(image_data,2)/2+exclude_region),spots(spots(:,10)==i,2),'-o','color',traj_color)
                    end
                    subplot(2,3,3)
                    hold on
                    plot(spots(spots(:,10)==i,9),spots(spots(:,10)==i,5),'-o','color',traj_color)
                    
                    spot_means(i)=mean(spots(spots(:,10)==i,5));
                    spot_sd(i)=std(spots(spots(:,10)==i,5));
                end
                title('Trajectory intensities vs frame num')
                subplot(2,3,4)
                
                %    hist(spot_means)
                hist(spots(:,5))
                title('Histogram Trajectory Intensity Mean')
                subplot(2,3,5)
                
                hist(spot_sd)
                title('Histogram Trajectory Intensity SD')
            else
                 if show_text_output==1
                disp('No trajectories found')
                 end
            end
        end
        %Assign spots to final variables
        if Ch==1
            if isempty(spots)==0
                SpotsCh1=spots;
            else
                SpotsCh1=0;
            end
        elseif Ch==2
            if isempty(spots)==0
                SpotsCh2=spots;
                % Transform Channel2 so spots are in the right place
                SpotsCh2(:,1)=SpotsCh2(:,1)+round(size(image_data,2)/2+exclude_region);
            else
                SpotsCh2=0;
            end
        end
      %  clear spots
    end
end
close(h)
end