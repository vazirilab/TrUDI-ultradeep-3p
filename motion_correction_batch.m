%% Motion correction pipeline using NoRMCorre
% Modified from demo.m
% Output: Motion corrected tiff files and motion correction setting files
% Author: Kalyan Santra

%%
clearvars -except selpath;

%% Add NoRMCorre library
addpath(genpath(['External', filesep, 'NoRMCorre']))

%% Get files
if exist("selpath", "var")
    selpath = uigetdir(selpath);
else
    selpath = uigetdir;
end
files = dir([selpath, filesep, '*.tif']);
nFiles = length(files);

%% Custom settings
motcor_settings.bin_width = 200;
motcor_settings.max_shift = 15;
motcor_settings.us_fac = 50;
motcor_settings.initial_batch_size = 200;
motcor_settings.correct_bidir = false; % Manual correction applied

%% run motion corr pipeline
run_motion_corr_pipeline(files, motcor_settings)

%%
function run_motion_corr_pipeline(files, motcor_settings)
    nFiles = length(files);
    for NN = 1:nFiles
        disp(['Starting motion correction for ', files(NN).name]);
        filename = [files(NN).folder, filesep, files(NN).name];
        [fpath, fname] = fileparts(filename);
        filename_mcor = fullfile(fpath, [fname, '_mcor.tif']);

        %%
        tic;
        Y = read_file(filename);
        toc; % read the file (optional, you can also pass the path in the function instead of Y)
        % Y = single(Y); % convert to single precision
        % T = size(Y, ndims(Y));
        % Y = Y - min(Y(:));

        %% set parameters (rigid motion correction)
        options_rigid = NoRMCorreSetParms( ...
            'd1', size(Y, 1), ...
            'd2', size(Y, 2), ...
            'bin_width', motcor_settings.bin_width, ...
            'max_shift', motcor_settings.max_shift, ...
            'us_fac', motcor_settings.us_fac, ...
            'init_batch', motcor_settings.initial_batch_size, ...
            'output_type', 'tiff', ... % save tiff
            'tiff_filename', filename_mcor, ...
            'correct_bidir', motcor_settings.correct_bidir);

        %% run motion correction
        tic;
        [~, ~, template, options_mc] = normcorre_batch_even(filename, options_rigid);
        toc

        %% Save settings and organize files
        filename_mat = [filename_mcor(1:end-4), '.mat'];
        save(filename_mat, 'motcor_settings', 'template', 'options_mc')
        pause(1)
        add_folder_tiff = [fpath, filesep, 'MCorr'];
        if ~exist(add_folder_tiff, 'dir')
            mkdir(add_folder_tiff)
        end
        movefile(filename_mcor, add_folder_tiff);

        add_folder_mat = [fpath, filesep, 'MCorr_mat'];
        if ~exist(add_folder_mat, 'dir')
            mkdir(add_folder_mat)
        end
        movefile(filename_mat, add_folder_mat);
    end
end

%%
%{
%% now try non-rigid motion correction (also in parallel)
options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);
tic; [M2,shifts2,template2,options_nonrigid] = normcorre_batch(Y,options_nonrigid); toc


%% compute metrics

nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);

[cY,mY,vY] = motion_metrics(Y,10);
[cM1,mM1,vM1] = motion_metrics(M1,10);
[cM2,mM2,vM2] = motion_metrics(M2,10);
T = length(cY);
%% plot metrics
figure;
    ax1 = subplot(2,3,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax2 = subplot(2,3,2); imagesc(mM1,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    ax3 = subplot(2,3,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
    linkaxes([ax1,ax2,ax3],'xy')
%% plot shifts        

shifts_r = squeeze(cat(3,shifts1(:).shifts));
shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
shifts_x = squeeze(shifts_nr(:,1,:))';
shifts_y = squeeze(shifts_nr(:,2,:))';

patch_id = 1:size(shifts_x,2);
str = strtrim(cellstr(int2str(patch_id.')));
str = cellfun(@(x) ['patch # ',x],str,'un',0);

figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')

%% plot a movie with the results

figure;
for t = 1:1:T
    subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(122);imagesc(M2(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.02);
end
%}