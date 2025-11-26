function varargout = Deinterleve(Img, numCh)
    %%% This function de-interleave the data in number of channel provide
    % as input;
    % Example use 2-ch: [ImgCh1, ImgCh2] = Deinterleve(Img, 2);
    % Author: Kalyan Santra
    %%
    if nargin < 1
        numCh = 1;
    end
    % check if number of output arguments matches number of channels
    if nargout ~= numCh
        error('Number of channels does not matches the number of output.')
    end
    if rem(size(Img, 3), numCh) ~= 0
        error('Number of total frames should be interger multiple of number of channels.')
    end
    varargout = cell(numCh, 1);
    for k = 1:numCh
        varargout{k} = Img(:, :, k:numCh:end);
    end
end
