%% ROI detection pipeline using CaImAn-MATLAB
% Modified from demo_script.m (CaImAn-MATLAB)
% Input: motion corrected tiff files
% Output: mat files with ROI data (all detected ROIs)
% Author: Kalyan Santra

%%
clearvars -except selpath;

%% Add CaImAn-MATLAB library
addpath(genpath(['External', filesep, 'CaImAn-MATLAB']))

%% Get files
if exist("selpath", "var")
    selpath = uigetdir(selpath);
else
    selpath = uigetdir;
end
files = dir([selpath, filesep, '*.tif']);
nFiles = length(files);

%% Custom settings
roidet_settings.input_components = 50; % number of components to find
roidet_settings.input_width = 8;
roidet_settings.frame_rate = 8.6;

roidet_settings.autoregression = 'Decay';
roidet_settings.init_method = 'Greedy';
roidet_settings.search_method = 'Ellipse';
roidet_settings.deconv_method = 'Constrained FOOPSI - SPGL1';

roidet_settings.input_merge_thresh = 0.9;
roidet_settings.input_fudge = 0.99;
roidet_settings.input_space_down = 1;
roidet_settings.input_time_down = 1;
roidet_settings.input_time_iter = 5;

roidet_settings.manual_refine = 0;
roidet_settings.use_classifier = 0;
roidet_settings.check_contours = 0;
roidet_settings.check_map = 0;

%% run analysis pipeline
run_roi_detection_pipeline(files, roidet_settings)

%%
function run_roi_detection_pipeline(files, roidet_settings)

    %% Set parameters

    K = roidet_settings.input_components; % number of components to be found
    tau = round(roidet_settings.input_width/2); % std of gaussian kernel (half size of neuron)
    %p = 2;
    switch roidet_settings.autoregression
        case 'No dynamics'
            p = 0;
        case 'Decay'
            p = 1;
        case 'Rise and Decay'
            p = 2;
    end

    switch roidet_settings.init_method
        case 'Greedy'
            init_method = 'greedy';
        case 'Sparse NMF'
            init_method = 'sparse_NMF';
    end

    switch roidet_settings.search_method
        case 'Ellipse'
            search_method = 'ellipse';
        case 'Dilate'
            search_method = 'dilate';
    end

    switch roidet_settings.deconv_method
        case 'Constrained FOOPSI - CVX'
            deconv_method = 'constrained_foopsi';
            method = 'cvx';
        case 'Constrained FOOPSI - SPGL1'
            deconv_method = 'constrained_foopsi';
            method = 'spgl1';
        case 'MCEM FOOPSI - CVX'
            deconv_method = 'MCEM_foopsi';
            method = 'cvx';
        case 'MCEM FOOPSI - SPGL1'
            deconv_method = 'MCEM_foopsi';
            method = 'spgl1';
        case 'MCMC'
            deconv_method = 'MCMC';
            method = 'cvx';
    end

    %%
    nFiles = length(files);
    for NN = 1:nFiles
        disp(['Starting ROI detection for ', files(NN).name]);
        filename = [files(NN).folder, filesep, files(NN).name];

        %% Read data files
        Y = read_file(filename);
        if ~isa(Y, 'single')
            Y = single(Y);
        end % convert to single
        Y = Y - min(Y(:));
        [d1, d2, T] = size(Y); % dimensions of dataset
        d = d1 * d2; % total number of pixels

        %%
        options = CNMFSetParms( ... % Check CNMFSetParms.m for more details
            'd1', d1, 'd2', d2, ... % Dimensions of datasets
            'ssub', roidet_settings.input_space_down, ... % Spatial downsampling
            'tsub', roidet_settings.input_time_down, ... % Temporal downsampling
            'init_method', init_method, ... % Initialization method
            'search_method', search_method, ... % Search locations when updating spatial components
            'dist', 3, ... % Expansion factor of ellipse
            'deconv_method', deconv_method, ... % Activity deconvolution method
            'method', method, ... % Method for solving convex problems
            'temporal_iter', roidet_settings.input_time_iter, ... % Number of block-coordinate descent steps
            'fudge_factor', roidet_settings.input_fudge, ... % Bias correction for autoregression coefficients
            'merge_thr', roidet_settings.input_merge_thresh, ... % Merging threshold
            'p', p, ... % Order of AR dynamics
            'nb', 1, ... % Number of background components
            'min_SNR', 1.5, ... % Minimum SNR threshold
            'cnn_thr', 0.2, ... % Threshold for CNN classifier
            'bas_nonneg', 0, ... % Allow a negative baseline
            'fr', roidet_settings.frame_rate, ... % frame rate
            'decay_time', 0.5, ...
            'rolling_length', ceil(roidet_settings.frame_rate*5));

        %% Data pre-processing
        [P, Y] = preprocess_data(Y, p, options);

        %% fast initialization of spatial components using greedyROI and HALS
        [Ain, Cin, bin, fin, center] = initialize_components(Y, K, tau, options, P); % initialize
        Cn = correlation_image(Y);

        %% manually refine components (optional)
        if roidet_settings.manual_refine
            [Ain, Cin, ~] = manually_refine_components(Y, Ain, Cin, center, Cn, tau, options);
        end

        %% update spatial components
        Yr = reshape(Y, d, T);
        [A, b, Cin] = update_spatial_components(Yr, Cin, fin, [Ain, bin], P, options);

        %% update temporal components
        P.p = 0; % set autoregression temporarily to zero for speed
        [C, f, P, S, YrA] = update_temporal_components(Yr, A, b, Cin, fin, P, options);

        %% classify components

        if roidet_settings.use_classifier

            % classify components
            rval_space = classify_comp_corr(Y, A, C, b, f, options);
            ind_corr = rval_space > options.space_thresh; % components that pass the correlation test
            % this test will keep processes

            % further classification with cnn_classifier
            try % matlab 2017b or later is needed
                [ind_cnn, ~] = cnn_classifier(A, [d1, d2], 'cnn_model', options.cnn_thr);
            catch
                ind_cnn = true(size(A, 2), 1); % components that pass the CNN classifier
            end

            % event exceptionality

            fitness = compute_event_exceptionality(C+YrA, options.N_samples_exc, options.robust_std);
            ind_exc = (fitness < options.min_fitness);

            % select components

            keep = (ind_corr | ind_cnn) & ind_exc;

            % display kept and discarded components
            figure;
            subplot(121);
            montage(extract_patch(A(:, keep), [d1, d2], [30, 30]), 'DisplayRange', [0, 0.15]);
            title('Kept Components');
            subplot(122); montage(extract_patch(A(:, ~keep), [d1, d2], [30, 30]), 'DisplayRange', [0, 0.15])
            title('Discarded Components');

            A = A(:, keep);
            C = C(keep, :);
        end

        %% merge found components
        [Am, Cm, ~, ~, Pm, ~] = merge_components(Yr, A, b, C, f, P, S, options);

        %% refine estimates excluding rejected components

        Pm.p = p; % restore autoregression value
        [A2, b2, C2] = update_spatial_components(Yr, Cm, f, [Am, b], Pm, options);
        [C2, f2, P2, S2, YrA2] = update_temporal_components(Yr, A2, b2, C2, f, Pm, options);

        %% Order ROI and plot

        [A_or, C_or, S_or, ~, srt] = order_ROIs(A2, C2, S2, P2); % order ROIs

        % plot ROI contours
        if roidet_settings.check_contours
            figure;
            plot_contours(A_or, Cn, options, 1);
        end

        % display ROI map
        if roidet_settings.check_map
            plot_components_GUI(Yr, A_or, C_or, b2, f2, Cn, options);
        end

        %% Save data
        % Calculate traces to be saved
        [F_raw, F_inferred, F_raw_binMask, F_dff, Bg_estimated, AddtnlOut] = construct_traces(Yr, A_or, C_or, b2, f2, YrA2(srt, :), options); % modified KS
        S_deconv = S_or;
        spatial_comp = A_or;
        ROI_center = zeros(size(spatial_comp, 2), 2);
        for cur_ROI = 1:size(spatial_comp, 2)
            single_ROI = full(reshape(spatial_comp(:, cur_ROI), size(Cn, 1), size(Cn, 2)));
            s = regionprops(true(size(single_ROI)), single_ROI, 'WeightedCentroid');
            ROI_center(cur_ROI, :) = s.WeightedCentroid;
        end

        % Save MAT file
        [fpath, fname] = fileparts(filename);
        add_folder_mat = [fpath, filesep, 'RoiDet'];
        if ~exist(add_folder_mat, 'dir')
            mkdir(add_folder_mat)
        end
        filename_mat = fullfile(add_folder_mat, [fname, '.mat']);

        if isfile(filename_mat)
            if ismember('F_raw', who('-file', filename_mat))
                warning(['File ', filename_mat, ' already have saved ROI detection data. Will save as a new filename to avoid overwriting.'])
                filename_mat = [filename(1:end-4), '_', datestr(now, 30), '.mat'];
                save(filename_mat, 'Cn', 'spatial_comp', 'F_raw', 'F_inferred', 'S_deconv', 'ROI_center', 'options', 'roidet_settings', 'F_raw_binMask', 'F_dff', 'Bg_estimated,', 'AddtnlOut');
            else
                save(filename_mat, 'Cn', 'spatial_comp', 'F_raw', 'F_inferred', 'S_deconv', 'ROI_center', 'options', 'roidet_settings', 'F_raw_binMask', 'F_dff', 'Bg_estimated', 'AddtnlOut', '-append');
            end
        else
            save(filename_mat, 'Cn', 'spatial_comp', 'F_raw', 'F_inferred', 'S_deconv', 'ROI_center', 'options', 'roidet_settings', 'F_raw_binMask', 'F_dff', 'Bg_estimated', 'AddtnlOut');
        end

    end
end

%%
function [F_raw, F_inferred, F_raw_binMask, F_dff, Bg_estimated, AddtnlOut] = construct_traces(Y, A, C, b, f, YrA, options)
    Y = double(Y);
    b = double(b);
    C = double(C);
    f = double(f);
    nA = full(sqrt(sum(A.^2))');
    [K, ~] = size(C);
    A = A / spdiags(nA, 0, K, K);
    C = bsxfun(@times, C, nA(:));
    AY = mm_fun(A, Y);

    Y_r = (AY - (A' * A) * C - full(A'*b) * f) + C;

    [~, Df] = extract_DF_F(Y, A, C, [], options, AY);

    F_raw = Y_r ./ Df;
    F_inferred = C ./ Df;

    %%
    Bg_estimated.sp_comp = b;
    Bg_estimated.temp_comp = f;

    AddtnlOut.Df = Df; % overall background level by component

    A_full = full(A);
    A_bin = zeros(size(A_full));
    A_bin(A > 0) = 1;
    F_raw_binMask = A_bin' * Y;

    [F_dff, F0] = detrend_df_f(A, b, C, f, YrA, options);
    AddtnlOut.F0 = F0;
end