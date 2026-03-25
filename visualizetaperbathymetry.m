function visualizetaperbathymetry(ddir_big, ddir_small, ddir_small_old)
% VISUALIZETAPERBATHYMETRY(ddir_big, ddir_small, ddir_small_old)
%
% Make summary plots from simulation outputs to test whether the new 
% bathymetry tapering method is better or worse than the old, simple cosine
% taper. See TAPERBATHYMETRY and TESTTAPERBATHYMETRY for the taper and
% setup implementation.
%
% INPUT:
% ddir_big              FK-SPECFEM3D run with big, original bathymetry
%                       acting as the ground truth
% ddir_small            FK-SPECFEM3D run with small, sliced bathmetry and
%                       tapered using the new method
% ddir_small_old        FK-SPECFEM3D run with small, sliced bathmetry and
%                       tapered using the old method
%
% SEE ALSO:
% TAPERBATHYMETRY, TESTTAPERBATHYMETRY
%
% Last modified by sirawich-at-prineton.edu, 03/25/2026

% read the interfaces (bathymetry)
itfs_big        = loadinterfacefiles3d(fullfile(ddir_big, 'DATA', ...
    'meshfem3D_files', 'interfaces.dat'));
itfs_small      = loadinterfacefiles3d(fullfile(ddir_small, 'DATA', ...
    'meshfem3D_files', 'interfaces.dat'));
itfs_small_old  = loadinterfacefiles3d(fullfile(ddir_small_old, 'DATA', ...
    'meshfem3D_files', 'interfaces.dat'));
bath_big        = itfs_big{1}.Z;
bath_small      = itfs_small{1}.Z;
bath_small_old  = itfs_small_old{1}.Z;

% compute gradient
[gx, gy] = gradient(bath_big, itfs_big{1}.SPACING_XI);
G_big = sqrt(gx.^2 + gy.^2);
[gx, gy] = gradient(bath_small, itfs_small{1}.SPACING_XI);
G_small = sqrt(gx.^2 + gy.^2);
[gx, gy] = gradient(bath_small_old, itfs_small_old{1}.SPACING_XI);
G_small_old = sqrt(gx.^2 + gy.^2);

% masking
G_big_sliced = G_big(81:177, 81:177);
mask = ones(size(G_big_sliced));
mask(25:73, 25:73) = 0;
mask = logical(mask);
G_big_sliced(25:73, 25:73) = NaN;
G_small(25:73, 25:73) = NaN;
G_small_old(25:73, 25:73) = NaN;

% read the seismograms
[tmh3d_big, pmh3d_big] = read_seismogram(fullfile(ddir_big, ...
    'OUTPUT_FILES','MH.P0009.HXP.semp'));
[tmh3d_small, pmh3d_small] = read_seismogram(fullfile(ddir_small, ...
    'OUTPUT_FILES', 'MH.P0009.HXP.semp'));
[tmh3d_small_old, pmh3d_small_old] = read_seismogram(fullfile(...
    ddir_small_old, 'OUTPUT_FILES', 'MH.P0009.HXP.semp'));

% compute the seismogram's sampling rate
fs_big = 1 / (tmh3d_big(2) - tmh3d_big(1));
fs_small = 1 / (tmh3d_small(2) - tmh3d_small(1));
fs_small_old = 1 / (tmh3d_small_old(2) - tmh3d_small_old(1));

% 2-pole 1-pass lowpass filter to 10 Hz
pmh3d_bigf = lowpass(pmh3d_big, fs_big, 10, 2, 1, 'butter', 'linear');
pmh3d_smallf = lowpass(pmh3d_small, fs_small, 10, 2, 1, 'butter', ...
    'linear');
pmh3d_small_oldf = lowpass(pmh3d_small_old, fs_small_old, 10, 2, 1, ...
    'butter', 'linear');

% running CC
t_running = (0:0.1:tmh3d_small(end))';
cc_running = nan(size(t_running).*[1 2]);
lag_running = nan(size(t_running).*[1 2]);
for ii = 50:length(t_running)
    xbig = pmh3d_bigf(and(tmh3d_big >= 0, tmh3d_big <= t_running(ii)));
    xsmall = pmh3d_smallf(and(tmh3d_small >= 0, ...
        tmh3d_small <= t_running(ii)));
    xsmall_old = pmh3d_small_oldf(and(tmh3d_small_old >= 0, ...
        tmh3d_small_old <= t_running(ii)));
    [r, lags] = xcorr(xbig, xsmall, round(fs_big), 'coeff');
    cc_running(ii,1) = max(r);
    lag_running(ii,1) = lags(r == max(r)) / fs_big;
    [r, lags] = xcorr(xbig, xsmall_old, round(fs_big), 'coeff');
    cc_running(ii,2) = max(r);
    lag_running(ii,2) = lags(r == max(r)) / fs_big;
end
cc_running(1:49,1) = cc_running(50,1);
cc_running(1:49,2) = cc_running(50,2);
lag_running(1:49,1) = lag_running(50,1);
lag_running(1:49,2) = lag_running(50,2);

% plot the bathymetry and the absolute gradient
figure(6)
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 8.5 7])
subplot('Position', [1/17 9/14 4/17 4/14])
imagesc([-1 1] * 32, [-1 1] * 32, bath_big/1000)
hold on
plot([-1 1 1 -1 -1]*6, [-1 -1 1 1 -1]*6, 'Color', 'k', 'LineWidth', 1.5)
axis xy
axis tight
axis equal
grid on
colormap('kelicol')
xlim([-1 1]*12)
ylim([-1 1]*12)
hold on
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'on', 'CLim', [-5 -3])
set(gca, 'XTick', get(gca, 'YTick'))
nolabels(gca, 1)
ylabel('Northing (km)')
title('Big bathymetry (zoomed in)')

subplot('Position', [6/17 9/14 4/17 4/14])
imagesc([-1 1] * 12, [-1 1] * 12, bath_small/1000)
hold on
plot([-1 1 1 -1 -1]*6, [-1 -1 1 1 -1]*6, 'Color', 'k', 'LineWidth', 1.5)
axis xy
axis tight
axis equal
grid on
colormap('kelicol')
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'on', 'CLim', [-5 -3])
set(gca, 'XTick', get(gca, 'YTick'))
nolabels(gca, 3)
title('Small bathymetry, new method')

subplot('Position', [11/17 9/14 5.5/17 4/14])
imagesc([-1 1] * 12, [-1 1] * 12, bath_small_old/1000)
hold on
plot([-1 1 1 -1 -1]*6, [-1 -1 1 1 -1]*6, 'Color', 'k', 'LineWidth', 1.5)
axis xy
axis tight
axis equal
grid on
colormap('kelicol')
cb = colorbar;
set(get(cb, 'Label'), 'String', 'elevation (km)')
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'on', 'CLim', [-5 -3])
set(gca, 'XTick', get(gca, 'YTick'))
nolabels(gca, 3)
title('Small bathymetry, old method')

subplot('Position', [1/17 4.25/14 4/17 4/14])
im = imagesc([-1 1] * 12, [-1 1] * 12, G_big_sliced);
set(im, 'AlphaData', mask)
axis xy
axis tight
axis equal
grid on
colormap('kelicol')
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'on', 'CLim', [0 max(G_big_sliced,[],'all')])
set(gca, 'XTick', get(gca, 'YTick'))
xlabel('Easting (km)')
ylabel('Northing (km)')
title('Gradient magnitude')

subplot('Position', [6/17 4.25/14 4/17 4/14])
im = imagesc([-1 1] * 12, [-1 1] * 12, G_small);
set(im, 'AlphaData', mask)
axis xy
axis tight
axis equal
grid on
colormap('kelicol')
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'on', 'CLim', [0 max(G_big_sliced,[],'all')])
set(gca, 'XTick', get(gca, 'YTick'))
nolabels(gca, 2)
xlabel('Easting (km)')
title('Gradient magnitude')

subplot('Position', [11/17 4.25/14 5.5/17 4/14])
im = imagesc([-1 1] * 12, [-1 1] * 12, G_small_old);
set(im, 'AlphaData', mask)
axis xy
axis tight
axis equal
grid on
colormap('kelicol')
cb = colorbar;
set(get(cb, 'Label'), 'String', 'Gradient magnitude')
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'on', 'CLim', [0 max(G_big_sliced,[],'all')])
set(gca, 'XTick', get(gca, 'YTick'))
nolabels(gca, 2)
xlabel('Easting (km)')
title('Gradient magnitude')

ax7 = subplot('Position', [1/17 1/14 4/17 2/14]);

histogram(G_big_sliced(mask))
grid on
xlim tight
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'on')
xlabel('Gradient magnitude')
ylabel('counts')

subplot('Position', [6/17 1/14 4/17 2/14])
histogram(G_small(mask))
grid on
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'on', ...
    'XLim', get(ax7, 'XLim'), 'YLim', get(ax7, 'YLim'))
nolabels(gca, 2)
xlabel('Gradient magnitude')

subplot('Position', [11/17 1/14 4/17 2/14])
histogram(G_small_old(mask))
grid on
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'on', ...
    'XLim', get(ax7, 'XLim'), 'YLim', get(ax7, 'YLim'))
nolabels(gca, 2)
xlabel('Gradient magnitude')

set(gcf, 'Renderer', 'painters')

%% plot the seismograms
figure(2)
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 8 9.5])

% title plot
subplot('Position', [1/16 18.4/19 14/16 0.01/19])
title('Title Goes Here')
set(gca, 'FontSize', 12, 'Color', 'none')
set(get(gca, 'XAxis'), 'Visible', 'off')
set(get(gca, 'YAxis'), 'Visible', 'off')

subplot('Position', [1/16 16/19 14/16 2/19])
plot(tmh3d_big, pmh3d_bigf, 'LineWidth', 1, 'Color', 'k')
grid on
xlim([0 50])
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)
nolabels(gca, 1)
title('Acoustic pressure')

subplot('Position', [1/16 13.5/19 14/16 2/19])
plot(tmh3d_small, pmh3d_smallf, 'LineWidth', 1, 'Color', [0.2 0.6 1])
grid on
xlim([0 50])
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)
nolabels(gca, 1)

subplot('Position', [1/16 11/19 14/16 2/19])
plot(tmh3d_small_old, pmh3d_small_oldf, 'LineWidth', 1, 'Color', [1 0.4 0.2])
grid on
xlim([0 50])
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)
nolabels(gca, 1)

subplot('Position', [1/16 8.25/19 14/16 2/19])
plot(tmh3d_small, pmh3d_smallf - pmh3d_bigf, 'LineWidth', 1, 'Color', [0.2 0.6 1])
hold on
plot(tmh3d_small, pmh3d_small_oldf - pmh3d_bigf, 'LineWidth', 1, 'Color', [1 0.4 0.2])
grid on
xlim([0 50])
ylim([-1 1]*max(ylim))
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)
xlabel('time (s)')
title('Difference')

subplot('Position', [1/16 4.25/19 14/16 2/19])
plot(t_running, cc_running(:,1), 'LineWidth', 1, 'Color', [0.2 0.6 1])
hold on
plot(t_running, cc_running(:,2), 'LineWidth', 1, 'Color', [1 0.4 0.2])
grid on
xlim([0 50])
a=1.15; ylim(ylim * [1+a 1-a; 1-a 1+a]/2)
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)
nolabels(gca, 1)
title('Running CC')

subplot('Position', [1/16 1.5/19 14/16 2/19])
plot(t_running, lag_running(:,1), 'LineWidth', 1, 'Color', [0.2 0.6 1])
hold on
plot(t_running, lag_running(:,2), 'LineWidth', 1, 'Color', [1 0.4 0.2])
grid on
xlim([0 50])
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)
xlabel('running time window length (s)')
title('Running time shift')

set(gcf, 'Renderer', 'painters')
end