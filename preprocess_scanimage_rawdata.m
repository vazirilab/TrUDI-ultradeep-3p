%% Raw image data pre-processing pipeline:
% 1. read image files, 2. extract channel(s), 3. correct bidirectional scan phase
% 4. Mean intensity projection (optional) 5. Create stack (optional)
% Input: raw image tiff file from scanimage software
% Output: pre-processed image data tiff file
% Based on scanimage setting the output options should be selected
% Author: Kalyan Santra

%%
clc;
clearvars -except selpath;
close all

%% Options
saveflag = 1; % to save all outputs; default 1
bidiphaseCorrShift = 1; % pixel unit
numCh = 2; % Depends of scan image setting; default 2
SelectImageCh = 1; % for inspection only
SaveMean = 0; % To save mean data projection along time; default 0
CombineMean2Stack = 0; % To combine single plane data to stack (if applicable); default 0

%% Get files 
if exist("selpath", "var")
    selpath = uigetdir(selpath);
else
    selpath = uigetdir;
end
files = dir([selpath, filesep, '*.tif']);
nFiles = length(files);

%%
for KK = 1:nFiles
    fig = figure(1);
    clf(fig)
    tiledlayout(fig, 1, 3);

    %% Read image
    disp('Loading data.')
    fullFileName = [files(KK).folder, filesep, files(KK).name];

    ImgRaw = ReadImageData(fullFileName);
    [ImgCh1, ImgCh2] = Deinterleve(ImgRaw, numCh);

    %% Inspect correction
    if SelectImageCh == 1
        imagedata = ImgCh1;
    elseif SelectImageCh == 2
        imagedata = ImgCh2;
    else
        error('Check Image channel number and deinterleve output')
    end

    img = mean(imagedata, 3);
    img = kron(img, ones(2)); % two fold upscale
    ax1 = nexttile;
    surf(img)
    view(2)
    shading("flat")
    daspect([1, 1, 1])
    title('shift = 0')

    ImageDataOut1 = CorrectBidiPhaseManual(imagedata, bidiphaseCorrShift);
    img2 = mean(ImageDataOut1, 3);
    img2 = kron(img2, ones(2)); % two fold upscale
    ax2 = nexttile;
    surf(img2)
    view(2)
    shading("flat")
    daspect([1, 1, 1])
    title(['shift = ', num2str(bidiphaseCorrShift)])

    ImageDataOut2 = CorrectBidiPhaseManual(imagedata, -bidiphaseCorrShift);
    img3 = mean(ImageDataOut2, 3);
    img3 = kron(img3, ones(2)); % two fold upscale
    ax3 = nexttile;
    surf(img3)
    view(2)
    shading("flat")
    daspect([1, 1, 1])
    title(['shift = ', num2str(-bidiphaseCorrShift)])
    linkaxes([ax1, ax2, ax3]);

    %% Apply correction
    if saveflag
        ImageDataOutCh1 = CorrectBidiPhaseManual(ImgCh1, bidiphaseCorrShift);
        ImageDataOutMeanCh1 = mean(ImageDataOutCh1, 3);

        ImageDataOutCh2 = CorrectBidiPhaseManual(ImgCh2, bidiphaseCorrShift);
        ImageDataOutMeanCh2 = mean(ImageDataOutCh2, 3);

        if CombineMean2Stack
            ImageDataOutStackCh1(:, :, KK) = ImageDataOutMeanCh1;
            ImageDataOutStackCh2(:, :, KK) = ImageDataOutMeanCh2;
        end
        fname = files(KK).name;
    end

    %% Save processed image channel
    if saveflag && ~CombineMean2Stack
        % save ch1
        if SaveMean
            data2save = ImageDataOutMeanCh1;
        else
            data2save = ImageDataOutCh1;
        end
        outfilename = [fname(1:end-4), '_phaseCorrCh1', '.tif'];
        Folder_add = [selpath, filesep, 'ScanPhaseCorrCh1'];
        if ~exist(Folder_add, 'dir')
            mkdir(Folder_add)
        end
        Output_folder_path = [Folder_add, filesep];
        WriteTiff(data2save, [Output_folder_path, outfilename])

        % save ch2
        if SaveMean
            data2save = ImageDataOutMeanCh2;
        else
            data2save = ImageDataOutCh2;
        end
        outfilename = [fname(1:end-4), '_phaseCorrCh1', '.tif'];
        Folder_add = [selpath, filesep, 'ScanPhaseCorrCh2'];
        if ~exist(Folder_add, 'dir')
            mkdir(Folder_add)
        end
        Output_folder_path = [Folder_add, filesep];
        WriteTiff(data2save, [Output_folder_path, outfilename])
    end
    disp(['File ', num2str(KK), ' of ', num2str(nFiles), ' done.']);
end

%% save mean data to stack (if applicable)
if saveflag && CombineMean2Stack
    % save ch1
    data2save = ImageDataOutStackCh1;
    outfilename = [fname(1:25), '_phaseCorrCh1Stack', '.tif'];
    Folder_add = [selpath, filesep, 'ScanPhaseCorrCh1'];
    if ~exist(Folder_add, 'dir')
        mkdir(Folder_add)
    end
    Output_folder_path = [Folder_add, filesep];
    WriteTiff(data2save, [Output_folder_path, outfilename])

    % save ch2
    data2save = ImageDataOutStackCh2;
    outfilename = [fname(1:25), '_phaseCorrCh2Stack', '.tif'];
    Folder_add = [selpath, filesep, 'ScanPhaseCorrCh2'];
    if ~exist(Folder_add, 'dir')
        mkdir(Folder_add)
    end
    Output_folder_path = [Folder_add, filesep];
    WriteTiff(data2save, [Output_folder_path, outfilename])
end
