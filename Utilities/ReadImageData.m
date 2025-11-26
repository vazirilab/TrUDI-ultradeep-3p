function imagedata = ReadImageData(filename)
    %%% This function reads the tiff file using ScanImageTiffReader library [platform dependent]
    % Req: ScanImageTiffReader
    % Author: Kalyan Santra
    %%    
    import ScanImageTiffReader.ScanImageTiffReader;
    hTif = ScanImageTiffReader(filename);
    imagedata = hTif.data();
    % Note: the above import switch the 1st and 2nd dimension of the data; this
    % is fixed here; 

    % check dimension and apply page transpose 
    if ndims(imagedata) > 3
        error('Check data format. Expect data in Y x X x T format or in Y x X format')
    end
    imagedata = pagetranspose(imagedata); % this correct the 1st and 2nd dimension flip

    % check input data type
    disp(['Input data type is ', class(imagedata), '. The final output type will be double.'])
    if isa(imagedata, 'uint16')
        imagedata = im2int16(imagedata);
    end
    imagedata = double(imagedata);
end