function [frame_average] = FrameAverage(image_data, noFrames, startFrame)

% Calculates a frame average of image_data array starting at startFrame
% over noFrames

frame_average=mean(image_data(:,:,startFrame:(startFrame+noFrames)),3);
frame_average=mat2gray(frame_average);
end