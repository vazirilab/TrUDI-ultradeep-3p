function std_below_percentile = std_below_percentile(X, p)
    % Calculate std of filter data
    % INPUTS
    % X:        input data
    % p:        percentile
    % Author: Kalyan Santra

    %%
    threshold = prctile(X, p, 2);
    threshold_matrix = repmat(threshold, 1, size(X, 2));
    X_below_threshold = X;
    X_below_threshold(X >= threshold_matrix) = NaN;
    std_below_percentile = std(X_below_threshold, 0, 2, 'omitnan');
end
