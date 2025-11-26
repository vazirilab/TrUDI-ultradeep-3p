%% ROI refine pipeline using CaImAn-MATLAB
% Based on demo_script.m (CaImAn-MATLAB)
% and EZcalcium (https://github.com/porteralab/EZcalcium)
% Input:  mat files with ROI data (from detection steps: see caiman_roi_detection_batch.m)
% Output: Updated mat files with refined ROI data appended
% Author: Kalyan Santra

%%
clearvars -except selpath;

%% Add CaImAn-MATLAB library
%addpath(genpath(["External", filesep, "CaImAn-MATLAB"]))
%% Get files
if exist("selpath", "var")
    selpath = uigetdir(selpath);
else
    selpath = uigetdir;
end
files = dir([selpath, filesep, '*.mat']);
nFiles = length(files);

%% ROI refinement key settings

roiref_settings.roundness_check = true;
roiref_settings.roundness = 0.8;

roiref_settings.oblongness_check = true;
roiref_settings.oblongness = 2;

roiref_settings.width_check = true;
roiref_settings.width_min = 5;
roiref_settings.width_max = 10;

roiref_settings.borderline_allow = 1;

roiref_settings.export_mat = true;

%% run analysis pipeline
run_roi_refine_pipeline(files, roiref_settings)

%%
function run_roi_refine_pipeline(files, roiref_settings)
    nFiles = length(files);
    for NN = 1:nFiles
        filename = [files(NN).folder, filesep, files(NN).name];

        %% run pipeline
        roiref = roiref_settings_init(roiref_settings);
        [roidata, border] = load_data(filename, roiref);
        disp(['Starting ROI refinement for ', files(NN).name]);
        roidata_keep = run_exclusion_check(roidata, border, roiref);
        disp(['Exporting data for ', files(NN).name]);
        export_data(roidata_keep, roiref)
    end
end

%%
function [roidata, border] = load_data(fullFileName, roiref)
    border = calc_border(roiref);
    roidata = load(fullFileName);
    roidata.path = fullFileName;

    if ~isfield(roidata, 'F_raw')
        warning('No ROI data found, check ROI data.');
        return
    end

    roidata.map_size = size(roidata.Cn);
    roidata.num_of_roi = size(roidata.F_raw, 1);
    roidata.num_of_frames = size(roidata.F_raw, 2);

    roidata = calc_roidata(roidata, roiref);

    nA = full(sqrt(sum(roidata.spatial_comp.^2))');
    [K, ~] = size(roidata.F_raw);
    roidata.spatial_comp = roidata.spatial_comp / spdiags(nA, 0, K, K);

    roidata.sx = min([roidata.options.sx, floor(roidata.map_size(1)/2), floor(roidata.map_size(2)/2)]);
    roidata.int_x = zeros(roidata.num_of_roi, 2*roidata.sx);
    roidata.int_y = zeros(roidata.num_of_roi, 2*roidata.sx);
    roidata.cm = com(roidata.spatial_comp, roidata.map_size(1), roidata.map_size(2));
end

%%
function roidata = run_exclusion_check(roidata, border, roiref)
    if roidata.num_of_roi < 1
        warning('No ROI data can be found. Check ROI analysis.')
        return
    end
    N = roidata.num_of_roi;
    keep = true(N, 1);
    for cur_roi = 1:N
        border_count = 0;

        if roiref.zscore_check && roidata.active_ROI(cur_roi) == 0
            keep(cur_roi) = false; continue
        end

        if roiref.baseline_stability_check
            if roidata.baseline_stability(cur_roi) > border.baseline_stability
                keep(cur_roi) = false; continue
            elseif roidata.baseline_stability(cur_roi) > roiref.baseline_stability
                border_count = border_count + 1;
            end
        end

        if roiref.roundness_check
            if roidata.roundness(cur_roi) < border.roundness
                keep(cur_roi) = false; continue
            elseif roidata.roundness(cur_roi) < roiref.roundness
                border_count = border_count + 1;
            end
        end

        if roiref.oblongness_check
            if roidata.oblongness(cur_roi) > border.oblongness
                keep(cur_roi) = false; continue
            elseif roidata.oblongness(cur_roi) > roiref.oblongness
                border_count = border_count + 1;
            end
        end

        if roiref.width_check
            if roidata.width(cur_roi) < border.width_min || roidata.width(cur_roi) > border.width_max
                keep(cur_roi) = false; continue
            elseif roidata.width(cur_roi) < roiref.width_min || roidata.width(cur_roi) > roiref.width_max
                border_count = border_count + 1;
            end
        end

        if roiref.area_check
            if roidata.area(cur_roi) < border.area_min || roidata.area(cur_roi) > border.area_max
                keep(cur_roi) = false; continue
            elseif roidata.area(cur_roi) < roiref.area_min || roidata.area(cur_roi) > roiref.area_max
                border_count = border_count + 1;
            end
        end

        if strcmp(roiref.skew_kurtosis_style, 'Raw')
            if roiref.skewness_check
                if roidata.skewness_raw(cur_roi) < border.skewness_min || roidata.skewness_raw(cur_roi) > border.skewness_max
                    keep(cur_roi) = false; continue
                elseif roidata.skewness_raw(cur_roi) < roiref.skewness_min || roidata.skewness_raw(cur_roi) > roiref.skewness_max
                    border_count = border_count + 1;
                end
            end

            if roiref.kurtosis_check
                if roidata.kurtosis_raw(cur_roi) < border.kurtosis_min || roidata.kurtosis_raw(cur_roi) > border.kurtosis_max
                    keep(cur_roi) = false; continue
                elseif roidata.kurtosis_raw(cur_roi) < roiref.kurtosis_min || roidata.kurtosis_raw(cur_roi) > roiref.kurtosis_max
                    border_count = border_count + 1;
                end
            end
        elseif strcmp(roiref.skew_kurtosis_style, 'Deconvolved')
            if roiref.skewness_check
                if roidata.skewness_deconv(cur_roi) < border.skewness_min || roidata.skewness_deconv(cur_roi) > border.skewness_max
                    keep(cur_roi) = false; continue
                elseif roidata.skewness_deconv(cur_roi) < roiref.skewness_min || roidata.skewness_deconv(cur_roi) > roiref.skewness_max
                    border_count = border_count + 1;
                end
            end

            if roiref.kurtosis_check
                if roidata.kurtosis_deconv(cur_roi) < border.kurtosis_min || roidata.kurtosis_deconv(cur_roi) > border.kurtosis_max
                    keep(cur_roi) = false; continue
                elseif roidata.kurtosis_deconv(cur_roi) < roiref.kurtosis_min || roidata.kurtosis_deconv(cur_roi) > roiref.kurtosis_max
                    border_count = border_count + 1;
                end
            end
        end

        if roiref.saturated_check
            if roidata.saturated_frames(cur_roi) > border.saturated_frames
                keep(cur_roi) = false; continue
            elseif roidata.saturated_frames(cur_roi) > roiref.saturated_frames
                border_count = border_count + 1;
            end
        end

        if border_count > roiref.borderline_allow
            keep(cur_roi) = false;
        elseif border_count > 0
            % do nothing
        else
            % do nothing
        end
    end
    % roidata.active_ROI   = keep;
    roidata.included_ROIs = find(keep(:)');
end

%%
function export_data(roidata_keep, roiref)
    roidata = roidata_keep;
    included_ROIs = roidata.included_ROIs;
    F_raw_refined = roidata.F_raw(roidata.included_ROIs, :);
    F_inferred_refined = roidata.F_inferred(roidata.included_ROIs, :);
    S_deconv_refined = roidata.S_deconv(roidata.included_ROIs, :);
    Z_mod = roidata.Z_mod;
    Z_mod_refined = roidata.Z_mod(roidata.included_ROIs, :);
    spatial_comp_refined = roidata.spatial_comp(:, roidata.included_ROIs);
    ROI_center_refined = roidata.ROI_center(roidata.included_ROIs, :);
    Addtnl_ref.F_raw_binMask_refined = roidata.F_raw_binMask(roidata.included_ROIs, :);
    Addtnl_ref.F_dff_refined = roidata.F_dff(roidata.included_ROIs, :);
    Addtnl_ref.Df_refined = roidata.AddtnlOut.Df(roidata.included_ROIs, :);
    Addtnl_ref.F0_refined = roidata.AddtnlOut.F0(roidata.included_ROIs, :);


    roiref_settings = roiref;
    % Save settings
    filename_mat_source = roidata.path;
    [fpath, fname] = fileparts(filename_mat_source);
    fpath_p = fileparts(fpath);
    add_folder_mat = [fpath_p, filesep, 'RoiRefine'];
    if ~exist(add_folder_mat, 'dir')
        mkdir(add_folder_mat)
    end
    filename_mat = fullfile(add_folder_mat, [fname, '.mat']);

    copyfile(filename_mat_source, filename_mat)

    if roiref.export_mat
        try
            if ismember('F_raw_refined', who('-file', filename_mat))
                warning(['File ', filename_mat, ' already have saved ROI refinement data. File name will be modified to avoid overwriting.'])
                filename_mat = [filename_mat(1:end-4), '_', datestr(now, 30), '.mat'];
                save(filename_mat, 'included_ROIs', 'F_raw_refined', 'F_inferred_refined', 'S_deconv_refined', 'Z_mod', 'Z_mod_refined', 'spatial_comp_refined', 'ROI_center_refined', 'roiref_settings', 'Addtnl_ref');
            else
                save(filename_mat, 'included_ROIs', 'F_raw_refined', 'F_inferred_refined', 'S_deconv_refined', 'Z_mod', 'Z_mod_refined', 'spatial_comp_refined', 'ROI_center_refined', 'roiref_settings', 'Addtnl_ref', '-append');
            end
        catch ME
            disp('If Error: copy the non-refined mat files into RoiRefine folder and then run export buttons')
            rethrow(ME)
        end
    else
        save(filename_mat, 'roiref_settings', '-append')
    end

    if roiref.export_csv
        writematrix(F_raw_refined, [filename_mat(1:end-4), '_F_raw_refined.csv'])
        writematrix(F_inferred_refined, [filename_mat(1:end-4), '_F_inferred_refined.csv'])
        writematrix(S_deconv_refined, [filename_mat(1:end-4), '_S_deconv_refined.csv'])
        writematrix(ROI_center_refined, [filename_mat(1:end-4), '_ROI_center_refined.csv'])
    end

    if roiref.export_xlsx
        filename_xlsx = [filename_mat(1:end-4), '.xlsx'];

        writematrix(F_raw_refined, filename_xlsx, 'Sheet', 'F_raw_refined')
        writematrix(F_inferred_refined, filename_xlsx, 'Sheet', 'F_inferred_refined')
        writematrix(S_deconv_refined, filename_xlsx, 'Sheet', 'S_deconv_refined')
        writematrix(ROI_center_refined, filename_xlsx, 'Sheet', 'ROI_center_refined')
    end
    disp('Data exported successfully.')
end

%%

%%
function roidata = calc_roidata(roidata, roiref)
    % Check if baseline stability frame window is valid
    baseline_loc = roiref.baseline_location / 100;
    if roidata.num_of_frames * baseline_loc < roiref.baseline_window
        roiref.baseline_window = floor(roidata.num_of_frames*baseline_loc);
    end

    % Calculate baseline stability - first vs last
    roidata.Z_mod = calc_zscore(roidata.F_raw, roiref.baseline_window); % calculate Z_mod

    first_frame_start = 1;
    first_frame_end = floor(baseline_loc*roidata.num_of_frames);
    last_frame_start = roidata.num_of_frames - first_frame_end + 1;
    last_frame_end = roidata.num_of_frames;

    [~, first_baseline, ~, ~] = calc_zscore(roidata.Z_mod(:, first_frame_start:first_frame_end), roiref.baseline_window);
    [~, last_baseline, ~, ~] = calc_zscore(roidata.Z_mod(:, last_frame_start:last_frame_end), roiref.baseline_window);
    roidata.baseline_stability = abs(first_baseline'-last_baseline');

    % Calculate significant activity
    significant_frames = roidata.Z_mod > roiref.zscore;
    roidata.active_ROI = max(movsum(significant_frames, roiref.zscore_window, 2, 'Endpoints', 'discard'), [], 2) == roiref.zscore_window;

    % Calculate saturated frames
    saturated_value = max(roidata.F_raw, [], 1);
    roidata.saturated_frames = zeros(1, roidata.num_of_roi);
    for i = 1:roidata.num_of_roi
        roidata.saturated_frames(i) = numel(find(roidata.F_raw(i, :) == saturated_value(i)));
        if roidata.saturated_frames(i) == 1 % if only one frame is at max value, assume it is not saturated
            roidata.saturated_frames(i) = 0;
        end
    end

    % Binarize image, extract stats
    roidata.area = zeros(1, roidata.num_of_roi);
    roidata.roundness = zeros(1, roidata.num_of_roi);
    roidata.width = zeros(1, roidata.num_of_roi);
    roidata.oblongness = zeros(1, roidata.num_of_roi);
    for i = 1:roidata.num_of_roi
        single_ROI = full(reshape(roidata.spatial_comp(:, i), roidata.map_size(1), roidata.map_size(2)));

        single_ROI = imbinarize(single_ROI, 0); % binary image with logical values
        image_stats = regionprops(single_ROI, 'Area', 'Perimeter');
        if length(image_stats) > 1 % check if ROI has discontiguous regions
            [~, region_ind] = max([image_stats.Area]);
            image_stats = image_stats(region_ind); % take only the biggest region
        end
        roidata.roundness(i) = (4 * pi * image_stats.Area) / (image_stats.Perimeter^2);

        single_ROI = double(single_ROI); % convert to numerical values
        image_stats = regionprops(single_ROI, 'Area', 'MajorAxisLength', 'MinorAxisLength');
        roidata.area(i) = image_stats.Area;
        roidata.width(i) = mean([image_stats.MajorAxisLength, image_stats.MinorAxisLength]);
        roidata.oblongness(i) = image_stats.MajorAxisLength / image_stats.MinorAxisLength;
    end

    roidata.kurtosis_raw = kurtosis(roidata.F_raw');
    roidata.kurtosis_deconv = kurtosis(roidata.S_deconv');

    roidata.skewness_raw = skewness(roidata.F_raw');
    roidata.skewness_deconv = skewness(roidata.S_deconv');
end

%%
function border = calc_border(roiref)
    border_percent = roiref.borderline / 100;

    border.zscore = roiref.zscore * (1 - border_percent);
    border.baseline_stability = roiref.baseline_stability * (1 + border_percent);
    border.roundness = roiref.roundness * (1 - border_percent);
    border.oblongness = roiref.oblongness * (1 + border_percent);
    border.saturated_frames = roiref.saturated_frames * (1 + border_percent);
    border.width_min = roiref.width_min * (1 - border_percent);
    border.width_max = roiref.width_max * (1 + border_percent);
    border.area_min = roiref.area_min * (1 - border_percent);
    border.area_max = roiref.area_max * (1 + border_percent);
    border.skewness_min = roiref.skewness_min * (1 - border_percent);
    border.skewness_max = roiref.skewness_max * (1 + border_percent);
    border.kurtosis_min = roiref.kurtosis_min * (1 - border_percent);
    border.kurtosis_max = roiref.kurtosis_max * (1 + border_percent);
end

%%
function [Z_mod, F_baseline_mean, F_baseline_std, baseline_start_frame] = calc_zscore(F, frames)
    [F_baseline_std, baseline_start_frame] = min(movstd(F, frames, 0, 2, 'Endpoints', 'discard'), [], 2);
    F_baseline_mean = mean(F(:, baseline_start_frame:baseline_start_frame+frames-1), 2);
    Z_mod = (F - F_baseline_mean) ./ F_baseline_std;
end

%%
function roiref = roiref_settings_init(roiref_settings)
    P = roiref_settings;
    defaults = struct( ...
        'zscore_check', false, 'zscore', 3, 'zscore_window', 5, ...
        'baseline_location', 25, 'baseline_window', 25, ...
        'baseline_stability_check', false, 'baseline_stability', 10, ...
        'roundness_check', false, 'roundness', 0.6, ...
        'oblongness_check', false, 'oblongness', 2.5, ...
        'saturated_check', false, 'saturated', 0, 'saturated_style', 'Frames', 'saturated_frames', 0, ...
        'width_check', false, 'width_min', 4, 'width_max', 40, ...
        'area_check', false, 'area_min', 100, 'area_max', 1000, ...
        'skewness_check', false, 'skewness_min', 0, 'skewness_max', 100, ...
        'kurtosis_check', false, 'kurtosis_min', 0, 'kurtosis_max', 100, ...
        'skew_kurtosis_style', 'Raw', ...
        'borderline_allow', 0, 'borderline', 10, ...
        'export_mat', true, 'export_csv', false, 'export_xlsx', false);
    f = fieldnames(defaults);
    for i = 1:numel(f)
        if ~isfield(P, f{i})
            P.(f{i}) = defaults.(f{i});
        end
    end
    roiref = P;
end