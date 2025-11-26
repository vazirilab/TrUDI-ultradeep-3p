function ImageDataOut = CorrectBidiPhaseManual(ImageData, bidiphase)
    %%% This function correct bi-directional phase offset in x-direction
    % Note: x-direction is 2nd dimension (col) if ReadImageData is used as
    % data importer
    % Author: Kalyan Santra 
    %%
    % check dimension
    if ndims(ImageData) > 3
        error('Check data format. Expect data in Y x X x T format or in Y x X format')
    end

    if nargin < 2
        bidiphase = 1;
    end

    ImageDataOut = shift(ImageData, bidiphase);

    function ImageData = shift(ImageData, bidiphase)
        % Shift last (which is 2nd axis) axis of "frames" by bidirectional phase offset
        % in-place; This fixes the X axis scan phase offset 
        if bidiphase > 0
            ImageData(1:2:end, bidiphase+1:end, :) = ImageData(1:2:end, 1:end-bidiphase, :);
        else
            ImageData(1:2:end, 1:end+bidiphase, :) = ImageData(1:2:end, -bidiphase+1:end, :);
        end
    end
end
