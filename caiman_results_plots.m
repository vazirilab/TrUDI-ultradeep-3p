%% Script for plotting CaImAn analysis results
% Input:  tiff image files and mat files for CaImAn analysis results from
% same folder
% Output: Plots and settings files
% Author: Kalyan Santra

%%
clc
clearvars -except selpath;
close all;

%% Options
Save_fig = 0;
save_settings = 0;
save_stat = 0;
x_axis_factor = (1 / 8.595094664371773) * (1 / 60); % frame rate (used high precision values)

maxNumTrace2show = 50;
maxFrame2show = 5000;
frameStartOffset = 0;

y_axis_magnify = 1;

% Plotting options
detrend_data = 0; % def: true (also if rolling_avg_trace and calc_DFF are false, detrend_data is true by default )
std_calc_full_window = 0; % def: false
rolling_avg_trace = 0; % def: false
rolling_zscore = 1;
calc_DFF = 0; % if true, trace plot will be DF/F and all zscore parameter will be ignored; if false the z-score will be calculated base on options; (No need for detrending/ optional)
d_prime_calc = 0; % If calc_DF/F true

%%
if exist("selpath", "var")
    selpath = uigetdir(selpath);
else
    selpath = uigetdir;
end

Files = dir([selpath, filesep, '*.mat']);
[parentFolder, ~, ~] = fileparts(selpath); % get parent folder
FilesTiff = dir([parentFolder, filesep, '*.tif']);

for k = 1 %1:size(Files, 1)
    FileName = Files(k).name;
    fullFileName = [Files(k).folder, filesep, Files(k).name];
    disp(FileName);
    load(fullFileName)
    if length(FileName) < 35
        short_filename = strrep(FileName(1:end-4), '~', '');

    else
        short_filename = strrep(FileName(1:30), '.', '');
    end

    %% Calculate DF/F, z-Score etc
    if exist('F_raw_refined', 'var')
        YY = F_raw_refined;
        YY2 = F_inferred_refined;
    else
        error('Filtered ROI data not found.')
        %warning('Filtered ROI data not found. Using all ROIs data.')
        YY = F_raw;
        YY2 = F_inferred;
    end
    % temporal smooth
    if rolling_avg_trace
        smooth_window = 10;
        smooth_dim = 2; % along T dimension
        smooth_method = "movmean";
        YY = smoothdata(YY, smooth_dim, smooth_method, smooth_window);
        YY2 = smoothdata(YY2, smooth_dim, smooth_method, smooth_window);
    end

    % detrend data
    if ~rolling_zscore && ~calc_DFF
        detrend_data = 1;
    end
    if detrend_data
        detrend_ptcl_p = 20;
        detrend_ptcl_window = 5000;
        detrend_ptcl_shift = detrend_ptcl_window;
        detrend_ptcl_mode = 1;
        YY = prctfilt_detrend(YY, detrend_ptcl_p, detrend_ptcl_window, detrend_ptcl_shift, detrend_ptcl_mode);
        YY2 = prctfilt_detrend(YY2, detrend_ptcl_p, detrend_ptcl_window, detrend_ptcl_shift, detrend_ptcl_mode);

    end
    % Calculate DF/F
    if calc_DFF
        Y = YY;
        Y2 = YY2;
        if d_prime_calc
            Y = YY .* sqrt(Addtnl_ref.F0_refined) * sqrt(0.2/2) * sqrt(1/2); % sqrt(tau/2), tau =  0.2 for 8s,   sqrt(1/2) for SBR
            Y2 = YY2 .* sqrt(Addtnl_ref.F0_refined) * sqrt(0.2/2) * sqrt(1/2);
        end
    else
        % Calculate standard deviation and z-score
        dim_zscore = 2;
        if ~rolling_zscore
            if std_calc_full_window
                YY_std = std(YY, 0, dim_zscore);
            else % filter by percentile 20% default
                std_ptcl_p = 20;
                YY_std = std_below_percentile(YY, std_ptcl_p);
            end
            zscore = YY ./ YY_std; % baseline subtraction is not needed for detrended data
            zscore_infer = YY2 ./ YY_std; % note: same std as raw
            Y = zscore;
            Y2 = zscore_infer;
        else
            window_zscore = 5000;
            [Rzscore, Rmean, Rstd] = zscore_calc_rolling(YY, dim_zscore, window_zscore);
            Y = Rzscore;
            Y2 = (YY2 - Rmean) ./ Rstd;
        end
    end

    %% Calc. trace stat
    pks_all = [];
    pks_thresold_ptcl = 90;
    pks_thresold = prctile(Y, pks_thresold_ptcl, "all");
    for KK = 1:size(Y, 1)
        % pks = findpeaks(Y(KK,:), "MinPeakHeight",0.5,"MinPeakProminence",0.5 );
        pks = findpeaks(Y(KK, :), "MinPeakHeight", pks_thresold);
        pks_all = [pks_all, pks];
    end
    fprintf('Max peak: %0.2f ;  mean: %0.2f ; median: % 0.2f \n', max(pks_all), mean(pks_all), median(pks_all))
    figure('Name', ['PeakStat_', short_filename]);
    histogram(pks_all)
    ylabel('Counts')
    if calc_DFF
        if d_prime_calc
            xlabel('d-prime peaks')
        else
            xlabel('DF/F peaks')
        end

    else
        xlabel('Z-scores peaks')
    end
    if save_stat
        Folder_add = [Files(1).folder, filesep, 'PostProcessStat'];
        if ~exist(Folder_add, 'dir')
            mkdir(Folder_add)
        end
        if calc_DFF
            if d_prime_calc
                setting_ops1 = 'dPrime';
            else
                setting_ops1 = 'DFF';
            end
            setting_ops2 = '';
        else
            setting_ops1 = 'zScore_';
            if rolling_zscore
                setting_ops2 = 'rolling';
            else
                setting_ops2 = 'detrended';
            end
        end

        Stat_filename = [Folder_add, filesep, 'PostProcessStat_', short_filename, '_options_', setting_ops1, setting_ops2];
        save(Stat_filename, "pks_thresold", "pks_thresold_ptcl", "pks_all")
    end

    %% Plot trace
    figure('Name', ['Int_', short_filename]);

    [sz1, sz2] = size(Y)
    if sz1 > 0
        frm_start = frameStartOffset + 1;
        frm_end = min(sz2, frameStartOffset+maxFrame2show);
        Frames = (1:min(sz2, maxFrame2show)) * x_axis_factor;
        Offset_factor = y_axis_magnify;
        MaxY = max(Y(:));
        MinY = min(Y(:));

        %% plot maximum of 5 traces
        numROI2plot = min(sz1, maxNumTrace2show);
        SkipList = [];  
        NoSkipList = [];
        ROICount = 0;
        accumarrayY = [];
        for NN = 1:sz1
            if ~ismember(NN, SkipList)
                y_offset = ROICount * MaxY * Offset_factor;
                plot(Frames, Y(NN, frm_start:frm_end)+y_offset, 'LineWidth', 3)
                hold on
                %plot(Frames, Y2(NN, frm_start:frm_end)+y_offset, 'color', [0.3, 0.6, 0.7], 'LineWidth', 1.0)

                ROICount = ROICount + 1;
                NoSkipList = [NoSkipList, NN];
                cutoff_thr(ROICount, :) = ones(1, sz2) + y_offset;
                zero_thr(ROICount, :) = zeros(1, sz2) + y_offset;
                traceNumLabel_y(ROICount) = y_offset;
                zscore_trace = Y(NN, :);
                zscore_filter = zscore_trace(zscore_trace > 1); % select activity as z-score greater than 1;
                accumarrayY = [accumarrayY; zscore_filter'];
            end

            if ROICount == numROI2plot
                break
            end
        end
        disp(['Max Z-score:', num2str(max(accumarrayY(:)))])
        disp(['Mean Z-score:', num2str(mean(accumarrayY(:)))])

        %%
        MaxYScale = MaxY * numROI2plot * Offset_factor + 1;

        plot(Frames, cutoff_thr(:, frm_start:frm_end), 'color', 'r', 'LineWidth', 2, 'LineStyle', ':')
        plot(Frames, zero_thr(:, frm_start:frm_end), 'color', 'k', 'LineWidth', 2, 'LineStyle', ':')

        % numbering trace
        for ii = 1:length(traceNumLabel_y)
            traceNum_x = round(-0.04*length(Frames)) * x_axis_factor;
            traceNum_y = traceNumLabel_y(ii) + (MaxY * Offset_factor * 0.1);
            text(traceNum_x, traceNum_y, [num2str(ii)], 'FontSize', 20, 'Color', [1, 0.5, 0], Rotation = 0, HorizontalAlignment = 'center', VerticalAlignment = 'bottom');

        end

        %%
        ax = gca;
        scale_bar_flag = 1;
        if scale_bar_flag
            scale_bar_height = 5;
            scale_bar_xpos = round(length(Frames)/100);
            scale_bar_xpos2 = round(scale_bar_xpos) * 2;
            x = ([0, 0] - scale_bar_xpos2) * x_axis_factor;
            y = [0, scale_bar_height] + MinY * 3;
            line(ax, x, y, 'LineWidth', 4, 'Color', 'k', 'LineStyle', '-');
            text(x(1)*4, (y(2) - y(1))/2+y(1), 'z-score', 'FontSize', 20, 'Color', 'k', Rotation = 90, HorizontalAlignment = 'center', VerticalAlignment = 'bottom');
            text(x(1)*1, (y(2) - y(1))/2+y(1), [num2str(y(2)-y(1))], 'FontSize', 20, 'Color', 'k', Rotation = 90, HorizontalAlignment = 'center', VerticalAlignment = 'bottom');

            scale_bar_width = 1 / x_axis_factor;
            x = ([0, scale_bar_width] - scale_bar_xpos2) * x_axis_factor;
            y = [0, 0] + MinY * 3;
            line(ax, x, y, 'LineWidth', 4, 'Color', 'k', 'LineStyle', '-');
            %text((x(2) - x(1))/2+x(1), (y(2) - y(1))/2+y(1), [num2str(x(2)-x(1))], 'FontSize', 20, 'Color', 'k', Rotation = 0, HorizontalAlignment = 'center', VerticalAlignment = 'bottom');
            text((x(2) - x(1))/2+x(1), 0+5*MinY, [num2str(x(2)-x(1)), ' min.'], 'FontSize', 20, 'Color', 'k', Rotation = 0, HorizontalAlignment = 'center', VerticalAlignment = 'bottom');

        end
        hold off
        set(ax, 'Clipping', 'off');
        axis_label_flag = 1;
        if axis_label_flag
            ylabel('Z-score')
            xlabel('Time (Min.)')
        else
            ylabel(' ')
            xlabel(' ')
        end
        xlim([0, max(Frames)])

        ylim([MinY, MaxYScale])

        %
        %ax.XTickLabel = {[]};
        % ax.YTick = [0, 200, 400];
        % ax.YTickLabel = {'0'; '200'; '400'};
        ax.YTickLabel = {[]};
        ax.FontSize = 20;

        axis off

        hFig = gcf;
        NewPos = hFig.Position;
        NewPos(3) = NewPos(3) .* 0.75;
        NewPos(4) = NewPos(4) .* 2;
        set(hFig, 'position', NewPos)
        if Save_fig

            Folder_add = [Files(k).folder, filesep, 'Time_trace_All_SI'];
            if ~exist(Folder_add, 'dir')
                mkdir(Folder_add)
            end
            OutFileName = [Folder_add, filesep, 'Trace_', short_filename];
            savefig(hFig, OutFileName)
            %close(hFig);
        end

        %%
        if 1
            pause(1)
            figure('Name', ['ROIs_Contour_', short_filename]);

            %%
            TiffFileName = [FilesTiff(k).folder, filesep, FilesTiff(k).name];

            ImRaw = ReadImageData(TiffFileName);

            summaryImage0 = mean(ImRaw, 3);

            Im_pixels_X = size(summaryImage0, 2);
            Im_pixels_Y = size(summaryImage0, 1);

            %% new summary image
            % ROIs = full(reshape(spatial_comp_refined, [], Im_pixels_X)); % get spatial foot print.
            ROIs = full(spatial_comp_refined);
            ROIs(ROIs > 0) = 1; % binarize
            Y_zs = Y;
            Y_zs(Y_zs > 1) = 1;
            Y_zs(Y_zs < 1) = 0;

            SpTempFilter = ROIs * Y_zs; % matrix multiplication

            zScoreFilter = ones(size(SpTempFilter));
            zScoreFilter(~(SpTempFilter > 0)) = NaN;
            zScoreFilter_re = reshape(zScoreFilter, Im_pixels_Y, Im_pixels_X, sz2);

            ImRaw_mod = zScoreFilter_re .* ImRaw;
            summaryImage1 = mean(ImRaw_mod, 3, 'omitnan');
            summaryImage_diff = summaryImage0 - summaryImage1;
            % mix
            summaryImage2 = summaryImage1;
            summaryImage2(isnan(summaryImage1)) = 0;
            summaryImage3 = summaryImage0;
            summaryImage3(~isnan(summaryImage1)) = 0;
            weight_mix = 1.1;
            summaryImage = summaryImage3 + (weight_mix .* summaryImage2); % weighted summary image
            %summaryImage = summaryImage0;
            %options.plot_bck_image = 1;
            [CC, jsf, im] = plot_contours_custom(spatial_comp_refined, summaryImage, options, 1, length(NoSkipList), [], [], NoSkipList, []);

            % adjust color map
            ax2 = gca;
            axCLmt = ax2.CLim;
            axCLmt_rng = axCLmt(2) - axCLmt(1);
            axCLmt_new = round([axCLmt_rng * 0.02, axCLmt_rng * 0.5]);
            ax2.CLim = axCLmt_new;

            %%
            hFig2 = gcf;
            if Save_fig
                Folder_add = [Files(k).folder, filesep, 'ROIsContour_All_SI'];
                if ~exist(Folder_add, 'dir')
                    mkdir(Folder_add)
                end
                OutFileName = [Folder_add, filesep, 'RoisContour_', short_filename];
                savefig(hFig2, OutFileName)
                %close(hFig);
            end
        end

    end
end

%% save settings files
if save_settings
    Folder_add = [Files(1).folder, filesep, 'PostProcessSettings'];
    if ~exist(Folder_add, 'dir')
        mkdir(Folder_add)
    end
    if calc_DFF
        setting_ops1 = 'DFF';
        setting_ops2 = '';
    else
        setting_ops1 = 'zScore_';
        if rolling_zscore
            setting_ops2 = 'rolling';
        else
            setting_ops2 = 'detrended';
        end
    end

    Settings_filename = [Folder_add, filesep, 'PostProcess_', short_filename, '_options_', setting_ops1, setting_ops2];
    save(Settings_filename)
end
