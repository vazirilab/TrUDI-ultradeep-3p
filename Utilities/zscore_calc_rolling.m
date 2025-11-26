function [rolling_zscore, rolling_mean, rolling_std] = zscore_calc_rolling(input_data, dim, window)
    % Need input data in k x T format, where k in the components, T is the
    % time steps
    if ndims(input_data) > 2
        error('Dimension greater than 2 is not supported in the current version.');
    end
    if dim ~= 2
        error('Check if the correct dimension is selected for rolling statistics.')
    end
    sz = size(input_data, dim);

    if window > sz + 1
        window = sz + 1;
    end

    if mod(window, 2) == 0
        window = window + 1;
    end
    disp(window)
    % Calculate the rolling mean and rolling standard deviation
    rolling_mean = movmean(input_data, window, dim);
    rolling_std = movstd(input_data, window, 0, dim);

    %% De-noise z-score at the end points using larger window than the default 'Shrink' window
    rolling_mean(:, 1:(window - 1)/2) = repmat(rolling_mean(:, (window - 1)/2+1), [1, (window - 1) / 2]);
    rolling_std(:, 1:(window - 1)/2) = repmat(rolling_std(:, (window - 1)/2+1), [1, (window - 1) / 2]);

    %%
    % Calculate the rolling Z-score
    rolling_zscore = (input_data - rolling_mean) ./ rolling_std;
end

%%
