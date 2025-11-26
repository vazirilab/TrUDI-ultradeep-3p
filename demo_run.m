%% Demo script to run selected pipeline
% The following section of the code demonstrate an example to process
% time-series data using CaImAn package
% Instruction to use the demo: "Run and Advance" each section separately
% Windows platform only
% Author: Kalyan Santra

%% Add data and code package paths
addpath(genpath(pwd));

%% Prepare raw data with default parameters
preprocess_scanimage_rawdata;
% Expected output: New sub-folder with pre-processed data as TIF file

%% Visualize motion analysis and corrected images
motion_correction_plots;
% Expected output: Figures with motion-corrected image and plot of motion
% correction metric

%% Apply motion correction with default parameters
motion_correction_batch;
% Expected output: New sub-folder with motion-corrected data as TIF and MAT file

%% Run CaImAn ROI detection on motion-corrected TIF data folder
caiman_roi_detection_batch;
% Expected output: New sub-folder with detected ROIs data for next step

%% Refine detected ROIs from CaImAn analysis based on specified ROI criteria on the detected ROI data folder
caiman_roi_refine_batch;
% Expected output: New sub-folder with refined ROIs data

%% Visualize CaImAn analysis traces and ROIs
caiman_results_plots;
% Expected output: Figures with ROI overlay on the image and time-trace
% plots