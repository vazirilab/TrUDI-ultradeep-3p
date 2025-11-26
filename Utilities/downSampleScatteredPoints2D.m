function [x_ds, y_ds] = downSampleScatteredPoints2D(x, y, varargin)
% Downsamples (x,y) based on density in 2D space.
%
% Usage:
%   [x_ds, y_ds] = downSampleScatteredPoints2D(x, y)
%   [x_ds, y_ds] = downSampleScatteredPoints2D(x, y, 'BinSize', 100, 'MaxPerBin', 50)
%
% Options:
%   'BinSize'    - Number of bins in x and y (default: 100)
%   'MaxPerBin'  - Max points to keep in each bin (default: 50)
% Author: Kalyan Santra

    p = inputParser;
    addParameter(p, 'BinSize', 100);
    addParameter(p, 'MaxPerBin', 50);
    parse(p, varargin{:});
    opts = p.Results;

    % Normalize into bin indices
    xEdges = linspace(min(x), max(x), opts.BinSize+1);
    yEdges = linspace(min(y), max(y), opts.BinSize+1);
    [~, xBin] = histc(x, xEdges);
    [~, yBin] = histc(y, yEdges);

    % Combine into 2D bin ID
    binID = xBin + opts.BinSize * yBin;

    % Group by bin
    [binGroups, ~, binIdx] = unique(binID);
    keepIdx = [];

    for i = 1:length(binGroups)
        thisBin = find(binIdx == i);
        if numel(thisBin) > opts.MaxPerBin
            % Subsample dense bins
            keepIdx = [keepIdx; randsample(thisBin, opts.MaxPerBin)];
        else
            % Keep all sparse bins
            keepIdx = [keepIdx; thisBin];
        end
    end

    x_ds = x(keepIdx);
    y_ds = y(keepIdx);
end