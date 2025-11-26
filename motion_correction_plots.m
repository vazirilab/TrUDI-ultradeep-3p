%% Motion correction pipeline using NoRMCorre
% Modified from demo.m
% To visualize motion correction of a selected data
% Author: Kalyan Santra

%%
clearvars -except selpath;

%% Add NoRMCorre libray
addpath(genpath(['External', filesep, 'NoRMCorre']))

%%
if exist("selpath", "var")
    selpath = uigetdir(selpath);
else
    selpath = uigetdir;
end

files = dir([selpath, filesep, '*.tif']);

nFiles = length(files);
Select_file = 1;
fullFileName = [files(Select_file).folder, filesep, files(Select_file).name];
disp(files(Select_file).name)

tic;
Y = read_file(fullFileName);
toc; % read the file (optional, you can also pass the path in the function instead of Y)
Y = single(Y); % convert to single precision
T = size(Y, ndims(Y));
Y = Y - min(Y(:));

%% set parameters (first try out rigid motion correction)

options_rigid = NoRMCorreSetParms('d1', size(Y, 1), 'd2', size(Y, 2), 'bin_width', 200, 'max_shift', 15, 'us_fac', 50, 'init_batch', 200, 'correct_bidir', false);

%% perform motion correction
tic;
[M1, shifts1, template1, options_rigid] = normcorre(Y, options_rigid);
toc

%% now try non-rigid motion correction (also in parallel)
% options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);
% tic; [M2,shifts2,template2,options_nonrigid] = normcorre_batch(Y,options_nonrigid); toc

%% compute metrics

nnY = quantile(Y(:), 0.005);
mmY = quantile(Y(:), 0.995);

[cY, mY, vY] = motion_metrics(Y, 10);
[cM1, mM1, vM1] = motion_metrics(M1, 10);
% [cM2,mM2,vM2] = motion_metrics(M2,10);
T = length(cY);

%% plot metrics
figure;
ax1 = subplot(2, 1, 1);
%imagesc(mY, [nnY, mmY]);
imagesc(mY, [min(mY(:)), max(mY(:))]);
axis equal;
axis tight;
axis off;
title('Before correction')
colormap gray
AddScaleBar(ax1, 50, 1.8, 20, 8, [1, 1, 1], 0)
ax2 = subplot(2, 1, 2);
%imagesc(mM1, [nnY, mmY]);
imagesc(mM1, [min(mM1(:)), max(mM1(:))]);
axis equal;
axis tight;
axis off;
title('After correction') %
colormap gray
AddScaleBar(ax2, 50, 1.8, 20, 8, [1, 1, 1], 0)
% ax3 = subplot(2,3,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
% subplot(2, 3, 4);
% % plot(1:T, cY, 1:T, cM1, 1:T, cM2);
% plot(1:T, cY, 1:T, cM1);
% legend('raw data', 'corrected');
% title('correlation coefficients', 'fontsize', 14, 'fontweight', 'bold')
% subplot(2, 3, 5);
% scatter(cY, cM1);
% hold on;
% plot([0.9 * min(cY), 1.05 * max(cM1)], [0.9 * min(cY), 1.05 * max(cM1)], '--r');
% axis square;
% xlabel('raw data', 'fontsize', 14, 'fontweight', 'bold');
% ylabel('rigid corrected', 'fontsize', 14, 'fontweight', 'bold');
% subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
%xlabel('rigid corrected', 'fontsize', 14, 'fontweight', 'bold');
%ylabel('non-rigid corrected', 'fontsize', 14, 'fontweight', 'bold');
%linkaxes([ax1,ax2,ax3],'xy')
linkaxes([ax1, ax2], 'xy')

hFig = gcf;
NewPos = hFig.Position;
NewPos(3) = NewPos(3) .* 0.75;
NewPos(4) = NewPos(4) .* 2;
set(hFig, 'position', NewPos)

%% plot shifts
pixel_res = 1.8; % um
shifts_r = squeeze(cat(3, shifts1(:).shifts));
maxShift = max(abs(shifts_r(:)));
% shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
% shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
% shifts_x = squeeze(shifts_nr(:,1,:))';
% shifts_y = squeeze(shifts_nr(:,2,:))';

% patch_id = 1:size(shifts_x, 2);
% str = strtrim(cellstr(int2str(patch_id.')));
% str = cellfun(@(x) ['patch # ', x], str, 'un', 0);

figure;
ax1 = subplot(311);
% plot(1:T, cY, 1:T, cM1, 1:T, cM2);
plot(1:T, cY, 1:T, cM1);
legend({'raw data', 'corrected'}, 'Box', 'off');
title('Correlation coefficients')
ylabel('Correlation coefficients')
ylim([0, 1])
%set(gca, 'Xtick', [])

%plot(shifts_x);
ax2 = subplot(312);
plot(shifts_r(:, 1)*pixel_res, 'linewidth', 2);
ylabel('Displacement (\mum)')
title('Displacements along x')
%set(gca, 'Xtick', [])
ylim([-maxShift, maxShift]*pixel_res)

%plot(shifts_x);
ax3 = subplot(313);
plot(shifts_r(:, 2)*pixel_res, 'linewidth', 2);
ylabel('Displacement (\mum)')
title('Displacements along y')
%set(gca, 'Xtick', [])
ylim([0, 1])

xlabel('Frames')
linkaxes([ax1, ax2, ax3], 'x')
ylim([-maxShift, maxShift]*pixel_res)

hFig = gcf;
NewPos = hFig.Position;
NewPos(3) = NewPos(3) .* 0.75;
NewPos(4) = NewPos(4) .* 2;
set(hFig, 'position', NewPos)

figure;
% plot(1:T, cY, 1:T, cM1, 1:T, cM2);
% plot(1:T, cY, 1:T, cM1);
%plot(1:T, -cY+cM1);
plot(cY, cM1, 'o')
hold on
plot([min(cY), max(cM1)], [min(cY), max(cM1)], '-r')

% legend({'raw data', 'corrected'}, 'Box', 'off');
title('Correlation coefficients')
xlabel('Raw correlation')
ylabel('Corrected correlation')
xlim([min(cY), max(cM1)])
ylim([min(cY), max(cM1)])
% set(gca, 'Xtick', [])

%% plot a movie with the results

% figure;
% for t = 1:1:T
%     subplot(121);
%     imagesc(Y(:, :, t), [nnY, mmY]);
%     xlabel('raw data', 'fontsize', 14, 'fontweight', 'bold');
%     axis equal;
%     axis tight;
%     title(sprintf('Frame %i out of %i', t, T), 'fontweight', 'bold', 'fontsize', 14); colormap('bone')
%     % subplot(122);imagesc(M2(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
%     title(sprintf('Frame %i out of %i', t, T), 'fontweight', 'bold', 'fontsize', 14); colormap('bone')
%     set(gca, 'XTick', [], 'YTick', []);
%     drawnow;
%     pause(0.02);
% end-