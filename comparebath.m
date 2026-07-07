function comparebath
% COMPAREBATH
%
% Compares bathymetry of CRUST1.0 and GEBCO_2020 as well as the sediment
% thickness of CRUST1.0 and GlobSed in a rectanglar box: 160 < LONS < 260,
% -40 < LATS < 10
%
% Last modified by spipatprathanporn@ucsd.edu, 07/06/2026

LONS = 160:1/240:260;
LATS = -40:1/240:10;

% crust1.0 bathymetry (ETOPO1) vs GEBCO bathymetry
elev_CRUST = readcrust1([], 'bnds', 2, LONS, LATS) * 1000;

[lons, lats, elev] = bathymetry([], [159.8 260.2], [-40.2 10.2], false);
[xx, yy] = meshgrid(LONS, LATS);
elev_GEBCO = interp2(lons, lats, double(elev'), xx, yy, 'linear');

% plot bathymetry comparison
figure(1)
set(gcf, 'Units', 'inches', 'Position', [0 1 8 12])
clf
subplot('Position', [0.08 0.69 0.84 0.26])
imagesc(LONS, LATS, elev_CRUST)
axis xy
axis tight
axis equal
grid on
cb = cax2dem([-7000 4000]);
set(get(cb, 'Label'), 'String', 'elevation (m)')
set(cb, 'TickDirection', 'out')
title('CRUST1.0: interpolated to 15x15 arc-second grid')
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)

subplot('Position', [0.08 0.37 0.84 0.26])
imagesc(LONS, LATS, elev_GEBCO)
axis xy
axis tight
axis equal
grid on
cb = cax2dem([-7000 4000]);
set(get(cb, 'Label'), 'String', 'elevation (m)')
set(cb, 'TickDirection', 'out')
title('GEBCO: interpolated to 15x15 arc-second grid')
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)

subplot('Position', [0.08 0.05 0.84 0.26])
imagesc(LONS, LATS, elev_CRUST - elev_GEBCO)
% imagesc(LONS, LATS, max(0, log10(abs(elev_CRUST - elev_GEBCO))) .* ...
%     sign(elev_CRUST - elev_GEBCO))
axis xy
axis tight
axis equal
grid on
colormap(gca, kelicol);
cb = colorbar;
set(get(cb, 'Label'), 'String', 'elevation difference (m)')
set(cb, 'TickDirection', 'out')
clim([-1 1] * 4000)
% clim([-4 4])
% set(cb, 'TickLabels', {-10000, -1000, -100, -10, '|x|<=1', 10, 100, 1000, 10000})
title('CRUST1.0 - GEBCO: interpolated to 15x15 arc-second grid')
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)

set(gcf, 'Renderer', 'painters')


% sediment thickness comparison
LONS = 160:1/12:260;
LATS = -40:1/12:10;

sed_CRUST = readcrust1([], 'thick', 10, LONS, LATS) * 1000;
fname_GlobSed = fullfile(getenv('IFILES'), 'EARTHMODELS', ...
    'PHYSICAL', 'GlobSed_v2.nc');
lons_GlobSed = ncread(fname_GlobSed, 'lon');
lats_GlobSed = ncread(fname_GlobSed, 'lat');
[~, i_lon_east] = min(abs(lons_GlobSed - LONS(1)));
[~, i_lon_west] = min(abs(lons_GlobSed - (LONS(end)-360)));
[~, i_lat1] = min(abs(lats_GlobSed - LATS(1)));
sed_GlobSed_east = ncread(fname_GlobSed, 'z', [i_lon_east i_lat1], ...
    [Inf length(LATS)])';
sed_GlobSed_west = ncread(fname_GlobSed, 'z', [2 i_lat1], ...
    [i_lon_west-1 length(LATS)])';
sed_GlobSed = [sed_GlobSed_east sed_GlobSed_west];
sed_Diff = sed_CRUST - sed_GlobSed;

% plot sediment thickness comparison
figure(2)
set(gcf, 'Units', 'inches', 'Position', [8 1 8 12])
clf
subplot('Position', [0.08 0.69 0.84 0.26])
im = imagesc(LONS, LATS, log10(max(1,sed_CRUST)));
set(im, 'AlphaData', ~isnan(sed_GlobSed))
axis xy
axis tight
axis equal
grid on
clim([0 4])
colormap(gca, kelicol);
cb = colorbar;
set(cb, 'TickLabels', {'<=1', 3.16, 10, 31.6, 100, 316, 1000, 3162, 10000})
set(get(cb, 'Label'), 'String', 'sediment thickness (m)')
set(cb, 'TickDirection', 'out')
title('CRUST1.0: interpolated to 5x5 arc-minute grid')
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)

subplot('Position', [0.08 0.37 0.84 0.26])
im = imagesc(LONS, LATS, log10(max(1,sed_GlobSed)));
set(im, 'AlphaData', ~isnan(sed_GlobSed))
axis xy
axis tight
axis equal
grid on
colormap(gca, kelicol);
clim([0 4])
cb = colorbar;
set(cb, 'TickLabels', {'<=1', 3.16, 10, 31.6, 100, 316, 1000, 3162, 10000})
set(get(cb, 'Label'), 'String', 'sediment thickness (m)')
set(cb, 'TickDirection', 'out')
title('GlobSed: 5x5 arc-minute grid')
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)

subplot('Position', [0.08 0.05 0.84 0.26])
im = imagesc(LONS, LATS, max(log10(abs(sed_Diff)), 0) .* sign(sed_Diff));
set(im, 'AlphaData', ~isnan(get(im, 'CData')))
axis xy
axis tight
axis equal
grid on
colormap(gca, kelicol);
cb = colorbar;
set(get(cb, 'Label'), 'String', 'sediment thickness difference (m)')
set(cb, 'TickDirection', 'out')
clim([-4 4])
set(cb, 'TickLabels', {-10000, -1000, -100, -10, '|x|<=1', 10, 100, 1000, 10000})
title('CRUST1.0 - GlobSed: interpolated to 5x5 arc-minute grid')
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)

set(gcf, 'Renderer', 'painters')
end