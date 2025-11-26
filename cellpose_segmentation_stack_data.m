%% Anatomical segmentation using cellpose MATLAB version
% this version take stack data (double format) as input for segmentation of ROI
% plane by plane; Then the original data is filtered based on the ROI mask;
% Output: Process data as tiff and mat files.
% Author: Kalyan Santra

%%
clc;
clearvars -except selpath;
close all

%%
saveFlag = 1;

%%
if exist("selpath", "var")
    [file, selpath] = uigetfile('*.tif', 'Select a file', selpath);
else
    [file, selpath] = uigetfile('*.tif', 'Select a file');
end

%% Read image
disp('Loading data.')
fullFileName = [selpath, filesep, file];
ImData = ReadImageData(fullFileName);

%%
fig1 = figure();
SV = sliceViewer(rescale(ImData), Parent = fig1);

%% Run single plane cellpose. For diagnosis only
% tic
% disp('Running cellpose on prepared data.')
% planeNumber = 1;
% imagedata = ImData(:, :, planeNumber);
% roiLabels = cellpose_segment(imagedata);
% disp(['Number of ROI detected: ', num2str(max(roiLabels(:)))])
% %
% imageLableOverlay = labeloverlay(rescale(imagedata), roiLabels);
% imshow(imageLableOverlay)
% toc

%% Get segmented ROIs on raw data/Input data
% [data type converted to double; although internally the cellpose model convert to single or re-scale further]

[sz1, sz2, sz3] = size(ImData);
roiLablesStack = zeros(sz1, sz2, sz3);
overlayStack = zeros(sz1, sz2, sz3, 3); % three RGB color for overlay
numRoisStack = zeros(sz3, 1);
tic
for NN = 1:sz3
    fprintf('Running cellpose on plane: %d \n', NN);
    imagedata = ImData(:, :, NN);
    roiLabels = cellpose_segment(imagedata);
    roiLablesStack(:, :, NN) = roiLabels;
    numRoisStack(NN) = max(roiLabels(:));
    %
    imageLableOverlay = labeloverlay(rescale(imagedata), roiLabels);
    overlayStack(:, :, NN, :) = reshape(imageLableOverlay, sz1, sz2, 1, 3);
    toc
end

%%
fig2 = figure();
SV2 = sliceViewer(roiLablesStack, Parent = fig2);
fig3 = figure();
SV3 = sliceViewer(uint8(overlayStack), Parent = fig3);

%% Filter based on ROI mask
% get ROI mask
roiMask = zeros(size(roiLablesStack));
roiMask(roiLablesStack > 0) = 1;
% filter
ImData_out = ImData .* roiMask;

%% Save data
%saveFlag = 0;
if saveFlag
    fname = file;
    % Save output image as .tif file
    data2save = ImData_out;
    outfilename1 = [fname(1:end-4), '_cellposeStackAnalysis', '.tif'];
    Folder_add = [selpath, filesep, 'Cellpose_stack_analysis'];
    if ~exist(Folder_add, 'dir')
        mkdir(Folder_add)
    end
    Output_folder_path = [Folder_add, filesep];
    WriteTiff(data2save, [Output_folder_path, outfilename1])
    % Save output as .mat file
    outfilename2 = [fname(1:end-4), '_cellposeStackAnalysis', '.mat'];
    save([Output_folder_path, outfilename2], "ImData", "roiLablesStack", "numRoisStack", "overlayStack", "roiMask", "ImData_out")
end

%%
function [labels] = cellpose_segment(img)
    cp = cellpose(Model = "cyto2");
    labels = segmentCells2D(cp, img, ImageCellDiameter = 6, FlowErrorThreshold = 0.8, CellThreshold = -2);
end