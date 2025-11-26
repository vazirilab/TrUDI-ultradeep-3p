function AddScaleBar(Ax, barSize, pixelSize, fontSize, lineWidth, sbColor,text_ON)
    % This functions add scale bar to the image 
    % Author: Kalyan Santra 
    %%
    % Defaults
    arguments
        Ax
        barSize = 50;
        pixelSize = 1.8
        fontSize = 30
        lineWidth = 8;
        sbColor = [1, 1, 1];
        text_ON = false;
    end
    % pixelSize = 1.8; % um
    % realLength = 50; % um
    Magnification = 1;
    interpolationRatio = 1;
    % sb_clr = [1, 1, 1];
    % Ax = gca;
    scaleBar(Ax, barSize, Magnification, pixelSize, sbColor, fontSize, lineWidth, interpolationRatio,text_ON)

    function scaleBar(Ax, realLength, Magnification, pixelSize, color, FontSz, LnWidth, interpolationRatio,text_ON)
    

        %%
        % for axis handle
        Im = findobj(Ax,'Type', 'Image');
        % [r, c, ~] = size(Im.CData);
        sz = size(Im.CData); r = sz(1); c = sz(2);
        %interpolationRatio = 1;%if image is not interpolated
        col1 = color;
        pixelCount = realLength * Magnification / (pixelSize * interpolationRatio);
        c_offset = round(c*(0.05));
        x = [c - pixelCount, c]- c_offset; 
        y = round([r * .95, r * .95]);
        line(Ax, x, y, 'LineWidth', LnWidth, 'Color', col1); hold on
        if text_ON
            text(x(1), round(r*.90), [num2str(realLength), '\mum'], 'FontWeight', 'bold', 'FontSize', FontSz, 'Color', col1);
        end
        hold off
    end
end