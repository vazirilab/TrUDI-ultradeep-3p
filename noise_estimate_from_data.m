%% Noise estimate from raw imaging data
% Input: Raw image tiff file
% Plots: scatter plot of the noise stat of pixels
%
% Author: Kalyan Santra

%%
clc;
clearvars -except selpath;
close all

%%
saveflag = 0;
bidiphaseCorrShift = 1; % pixel
numCh = 2;
SelectImageCh = 1;
improvePlotRendering = 1;
binSize = 150;
maxPointsPerBin = 40;
% Adjust binSize and maxPointsPerBin to improve scatter plot rendering while keeping extreme feature

%% Get files
if exist("selpath", "var")
    [file, selpath] = uigetfile('*.tif', 'Select a file', selpath);
else
    [file, selpath] = uigetfile('*.tif', 'Select a file');
end

%% Read image
disp('Loading data.')
fullFileName = [selpath, filesep, file];

%% Get data from ch and correct bidirectional scan phase
ImgRaw = ReadImageData(fullFileName);
if numCh == 2
    [ImgCh1, ImgCh2] = Deinterleve(ImgRaw, numCh);
    if SelectImageCh == 1
        imagedata = ImgCh1;
    elseif SelectImageCh == 2
        imagedata = ImgCh2;
    else
        error('Check Image channel number and deinterleve output')
    end
else
    imagedata = ImgRaw;
end
imagedata_corr = CorrectBidiPhaseManual(imagedata, bidiphaseCorrShift);

digi_offset_corr = 0; % Digitizer offset (if any): default: 0
if digi_offset_corr
    imagedata_corr = imagedata_corr - min(imagedata_corr(:)); %
end

%% calculate stat and plot results
stat = calculate_stats(imagedata_corr);
plot_results(stat, improvePlotRendering, binSize, maxPointsPerBin);

%%

%%
function stat = calculate_stats(imdata, dim, window, keep_length, trimEdgeWidth, negative_filter)
    % default options
    if nargin < 2
        dim = 3;
    end
    if nargin < 3
        window = 100; % number of from for moving stat
    end
    if nargin < 4
        keep_length = 100; % length of frame to keep
    end
    if nargin < 5
        trimEdgeWidth = 5; % width of the image edge to trim to remove phase correction artifacts
    end
    if nargin < 6
        negative_filter = 1;
    end
    sz = size(imdata);
    % calculate state
    Mean_img = movmean(imdata, window, dim);
    Std_img = movstd(imdata, window, 0, dim);
    Var_img = movvar(imdata, window, 0, dim);

    % Trim stat data calculated;
    t_idx = window:window + keep_length;
    x_idx = (trimEdgeWidth + 1):(sz(2) - trimEdgeWidth);
    y_idx = (trimEdgeWidth + 1):(sz(1) - trimEdgeWidth);

    Mean_img_trm = Mean_img(y_idx, x_idx, t_idx);
    Std_img_trm = Std_img(y_idx, x_idx, t_idx);
    Var_img_trm = Var_img(y_idx, x_idx, t_idx);

    % filter negative mean
    if negative_filter
        stat.mean = Mean_img_trm(Mean_img_trm > 0);
        stat.std = Std_img_trm(Mean_img_trm > 0);
        stat.var = Var_img_trm(Mean_img_trm > 0);
    else
        stat.mean = Mean_img_trm;
        stat.std = Std_img_trm;
        stat.var = Var_img_trm;
    end
end

%%
function plot_results(stat, improvePlotRendering, binSize, maxPointsPerBin)
    if nargin < 2
        improvePlotRendering = 1;
        binSize = 150;
        maxPointsPerBin = 40;
    end
    data_mean = stat.mean;
    data_std = stat.std;
    data_var = stat.var;
    %%%
    clr_fit = 'r';
    clr_model = 'g';
    %%% SNR vs Mean plot
    figure
    x = data_mean(:);
    y = (data_mean(:) ./ data_std(:)).^1;
    if improvePlotRendering
        [x_ds, y_ds] = downSampleScatteredPoints2D(x, y, 'BinSize', binSize, 'MaxPerBin', maxPointsPerBin);
        x_few = linspace(min(x), max(x), 1001);
        scatter(x_ds, y_ds)
    else
        scatter(x, y)
    end
    %P = polyfit(x, y, 1); % check linear fit
    %P_para = polyfit(y, x, 2); % check parabolic fit
    try
        [fitresult, ~] = createFit_SNR(x, y);
        hold on;
        if improvePlotRendering
            yfit_few = feval(fitresult, x_few);
            hAx = plot(x_few, yfit_few);
            pixel_count = x_few;
        else
            hAx = plot(fitresult);
            pixel_count = x;
        end
        hAx.LineWidth = 2;
        hAx.Color = clr_fit;
        model_gamma = 16.6;
        SNR_model = (1 / model_gamma) * sqrt(pixel_count);
        plot(pixel_count, SNR_model, LineStyle = "--", Color = clr_model)
        hold off
        ylim([0, 4])
        xlim([0, 1500])
        xlabel('Mean pixel count')
        ylabel('Pixel SNR')
        lgd = legend({'Data'; 'Data Fit'; 'Model'}, "Box", "off", Location = "northwest");
        %grid on
        title('Pixel SNR: experiment vs. model')

        legend_entries = findobj(lgd, 'Type', 'line'); % includes lines and markers

        for k = 1:length(legend_entries)
            if strcmp(get(legend_entries(k), 'Marker'), '.')
                set(legend_entries(k), 'MarkerSize', 32); 
            end
        end
    catch
    end

    %%% Variance vs Mean
    figure
    x = data_mean(:);
    y1 = data_var(:);
    if improvePlotRendering
        [x1_ds, y1_ds] = downSampleScatteredPoints2D(x, y1, 'BinSize', binSize, 'MaxPerBin', maxPointsPerBin);
        x_few = linspace(min(x), max(x), 1001);
        scatter(x1_ds, y1_ds)
    else
        scatter(x, y1)
    end
    P1 = polyfit(x, y1, 1); %check linear fit
    hold on;
    if improvePlotRendering
        yfit1_few = polyval(P1, x_few);
        plot(x_few, yfit1_few, LineStyle = '-', color = clr_fit, LineWidth = 2);
        pixel_count = x_few;
    else
        yfit1 = polyval(P1, x);
        plot(x, yfit1, LineStyle = '-', color = clr_fit, LineWidth = 2);
        pixel_count = x;
    end

    model_gamma = 16.6;
    var_model = (model_gamma.^2) * pixel_count;
    plot(pixel_count, var_model, LineStyle = "--", Color = clr_model)

    hold off
    %grid on
    xlim([0, 1500])
    ylim([0, 5e5])
    xlabel('Mean pixel count')
    ylabel('Variance')
    lgd = legend({'Data'; 'Data Fit'; 'Model'}, "Box", "off", Location = "northwest");
    ax = gca;
    ax.YAxis.Exponent = 0;
    ytickformat('%1.0e')
    title('Variance of pixel count: experiment vs. model')
    ylim([0, 5.2e5])
    %yticklabels(arrayfun(@(v) sprintf('%.0e', v), (0:5)*1e5, 'UniformOutput', false))

    legend_entries = findobj(lgd, 'Type', 'line'); 
    for k = 1:length(legend_entries)
        if strcmp(get(legend_entries(k), 'Marker'), '.')
            set(legend_entries(k), 'MarkerSize', 32); 
        end
    end

    %%% distribution of SNR
    figure
    clr = lines(10);
    histogram(y, Normalization = 'probability', EdgeColor = "none", FaceColor = clr(2, :))
    xlabel('Pixel SNR')
    ylabel('Probability')
    % title(sprintf('SNR distribution: Mean SNR = %0.2f', mean(y)))
    title({'Histogram of experimental pixel SNR'; ' for 2 nJ pulse energy'})
    ytickformat('%0.3f')
    legend({['Mean pixel SNR = ', num2str(mean(y), 2)]}, "Box", "off", "IconColumnWidth", 0.0001)
    %grid on
end

%%
function [fitresult, gof] = createFit_SNR(x, y)
    [xData, yData] = prepareCurveData(x, y);

    % Set up fittype and options.
    ft = fittype('power1');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Display = 'Off';
    opts.StartPoint = [0.01, 0.8];

    % Fit model to data.
    [fitresult, gof] = fit(xData, yData, ft, opts);
end
