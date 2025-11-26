function WriteTiff(imagedata, filename)
    %%% This function write data as tiff file in int16 format
    % Author: Kalyan Santra
    %%
    % Warning message if data is out of range 
    mn = min(imagedata(:));
    mx = max(imagedata(:));
    if mn < -32768
        warning('Write image data array to Tiff. Minimum value %0.3f is outside the range [-32768 32767]. All values smaller than -32768 will be converted to -32768.', mn)
    end
    if mx > 32767
        warning('Write image data array to Tiff. Maximum value %0.3f is outside the range [-32768 32767]. All values greater than 32767 will be converted to 32767.', mx)
    end

    np2 = nextpow2(mx-mn);
    if np2 < 8
        warning('Data spans less than %d-bit out of 16-bit. Cosider rescaling or normalization of data.', np2)
    end

    % Implement tiff write
    imagedata = int16(imagedata);
    %faster and datatype preserved
    t = Tiff(filename, 'w');
    for fm = 1:size(imagedata, 3)
        t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
        % t.setTag('Compression',Tiff.Compression.None);
        t.setTag('BitsPerSample', 16);
        t.setTag('SamplesPerPixel', 1);
        t.setTag('SampleFormat', Tiff.SampleFormat.Int);
        % t.setTag('ExtraSamples',Tiff.ExtraSamples.Unspecified);
        t.setTag('ImageLength', size(imagedata, 1));
        t.setTag('ImageWidth', size(imagedata, 2));
        % t.setTag('TileLength',32);
        % t.setTag('TileWidth',32);
        t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);

        t.write(imagedata(:, :, fm));
        t.writeDirectory
    end
    t.close();
end