%% Temporal smoothing of motion and phase corrected data
% Input: motion and bidirectional scan phase corrected tiff file
% Output: Smooth data tiff file, setting mat file
% Author: Kalyan Santra

%%
clc;
clearvars -except selpath;

%% data smoothing options
saveflag = 1; %to save all outputs; default 1
smooth_window = 5; % rolling window size
smooth_dim = 3; % along T dimension
smooth_method = "movmean"; % options = { "movmean"; "gaussian"}

%% Get files
if exist("selpath", "var")
    selpath = uigetdir(selpath);
else
    selpath = uigetdir;
end
files = dir([selpath, filesep, '*.tif']);
nFiles = length(files);

%%
for KK = 1:size(files)
    fullFileName = [files(KK).folder, filesep, files(KK).name];
    imagedata = ReadImageData(fullFileName);

    %% smooth data save
    if saveflag
        ImDataSmooth = smoothdata(imagedata, smooth_dim, smooth_method, smooth_window);
        data2save = ImDataSmooth;

        %%
        switch smooth_method
            case "gaussian"
                fileNameAdd = ['_smoothGauss', num2str(smooth_window)];
                folderNameAdd = ['SmoothGauss', num2str(smooth_window)];
            case "movmean"
                fileNameAdd = ['_smoothMovMean', num2str(smooth_window)];
                folderNameAdd = ['SmoothMovMean', num2str(smooth_window)];
        end

        %%
        fname = files(KK).name;
        outfilename = [fname(1:end-4), fileNameAdd, '.tif'];
        Folder_add = [selpath, filesep, folderNameAdd];
        if ~exist(Folder_add, 'dir')
            mkdir(Folder_add)
        end
        Output_folder_path = [Folder_add, filesep];
        WriteTiff(data2save, [Output_folder_path, outfilename])
        disp(['File ', num2str(KK), ' of ', num2str(nFiles), ' done.'])
    end
end

%% Save the smooth settings
if saveflag
    SmoothSetting.smoothFlag = saveflag;
    SmoothSetting.smooth_window = smooth_window;
    SmoothSetting.smooth_dim = smooth_dim;
    SmoothSetting.smooth_method = smooth_method;
    Settings_outfilename = 'DataSmoothSetting.mat';
    Folder_add = [selpath, filesep, folderNameAdd, filesep, 'SmoothSetting'];
    if ~exist(Folder_add, 'dir')
        mkdir(Folder_add)
    end
    Output_folder_path = [Folder_add, filesep];
    save([Output_folder_path, Settings_outfilename], 'SmoothSetting')
    disp('Setting file saved.')
end

%%
