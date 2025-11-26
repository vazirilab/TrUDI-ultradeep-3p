%% Cellpose results refinement based on ROI properties
% To filter segmented cell by properties
% Input: Cellpose results (segmented stack data)
% Output: Filtered data
% Author: Kalyan Santra

%%
clc;
clc;
clearvars -except selpath;

%%
saveFlag = 1;

%%
% select mat files containing cellpose analysis
if exist("selpath", "var")
    [file, selpath] = uigetfile('*.mat', 'Select a file', selpath);
else
    [file, selpath] = uigetfile('*.mat', 'Select a file');
end

%% load data
disp('Loading data.')
fullFileName = [selpath, filesep, file];
load(fullFileName);

%%
[numRows, numCols, numPlanes] = size(roiLablesStack);
roiLablesStack_filter = zeros(size(ImData));
for NN = 1:numPlanes
    currentPlane = roiLablesStack(:, :, NN);
    props = regionprops(currentPlane, 'Centroid', 'EquivDiameter', 'Perimeter', 'Area', 'PixelIdxList');

    % dia = {props.EquivDiameter};
    % dia = cell2mat(dia);
    % sumSmallROI(NN)= sum(dia<5.5)
    minDiameter = 6; % pix
    maxDiameter = 10; % pix
    for MM = 1:numel(props)
        if props(MM).EquivDiameter >= minDiameter && props(MM).EquivDiameter <= maxDiameter
            % Add the region to the filtered mask
            roiLablesStack_filter(props(MM).PixelIdxList+(NN - 1)*numRows*numCols) = currentPlane(props(MM).PixelIdxList);
        end
    end
end

% figure;
% subplot(1, 2, 1);
% imshow(label2rgb(roiLablesStack(:, :, 100), 'jet', 'k', 'shuffle'));
% title('Original Plane 1');
%
% subplot(1, 2, 2);
% imshow(label2rgb(roiMask_filter(:, :, 100), 'jet', 'k', 'shuffle'));
% title('Filtered Plane 1');

%%
% get ROI mask
roiMask_filter = zeros(size(roiLablesStack_filter));
roiMask_filter(roiLablesStack_filter > 0) = 1;
% filter
ImData_filter = ImData .* roiMask_filter;

%% Save data
if saveFlag
    fname = file;
    selpath = pwd;
    % Save output image as .tif file
    data2save = ImData_filter;
    outfilename1 = [fname(1:end-4), '_cellposeStackAnalysis_filter', '.tif'];
    Folder_add = [selpath, filesep, 'Cellpose_stack_analysis_filter'];
    if ~exist(Folder_add, 'dir')
        mkdir(Folder_add)
    end
    Output_folder_path = [Folder_add, filesep];
    WriteTiff(data2save, [Output_folder_path, outfilename1])
end