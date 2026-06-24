function plotcrust1(crust1dir)
% plotcrust1(crust1dir)
%
% Plots maps of CRUST1.0 model parameters.
%
% CRUST1.0 reference:
% Laske, G., Masters, G., Ma Z. & Pasyanos M. (2013). Update on
%     CRUST1.0--A 1-degree global of Earth's crust. Paper presented at the
%     Geophys. Res. Abstracts
%
% The model is available at https://igppweb.ucsd.edu/~gabi/crust1.html
%
% INPUT:
% crust1dir         file directory of CRUST1.0 model
%
% Last modified by sirawich@princeton.edu, 06/18/2026

defval('crust1dir', fullfile(getenv('IFILES'), 'EARTHMODELS', ...
    'PHYSICAL', 'crust1.0'))

% constants
NUMLAYERS = 9;
NUMLONS = 360;
NUMLATS = 180;
LONLIM = [-180 180];
LATLIM = [-90 90];

% layer names (see crust1dir/readme)
LAYERNAMES = {'water', ...
               'ice', ...
               'upper sediments', ...
               'middle sediments', ...
               'lower sediments', ...
               'upper crystalline crust', ...
               'middle crystalline crust', ...
               'lower crystalline crust', ...
               'mantle'};
% interface names
BNDNAMES = {'top of water', ...
             'bottom of water', ...
             'bottom of ice', ...
             'bottom of sediments 1', ...
             'bottom of sediments 2', ...
             'bottom of sediments 3', ...
             'bottom of crystalline crust 1', ...
             'bottom of crystalline crust 2', ...
             'bottom of crystalline crust 3 (Moho)'};

% figure positions
FIGPOS = repmat([0.025 0.06 0.29 0.23], [NUMLAYERS 1]);
FIGPOS = FIGPOS + 1/3 * [0 2 0 0; ...
                         1 2 0 0; ...
                         2 2 0 0; ...
                         0 1 0 0; ...
                         1 1 0 0; ...
                         2 1 0 0; ...
                         0 0 0 0; ...
                         1 0 0 0; ...
                         2 0 0 0];

%% RHO: Density
% read the data files
fid = fopen(fullfile(crust1dir, 'crust1.rho'));
rho = fscanf(fid, '%f', [NUMLAYERS Inf]);
rho = rho';
fclose(fid);

% make a plot
figure(1)
set(gcf, 'Units', 'inches', 'Position', [0 1 12 6])
clf
for ii = 1:NUMLAYERS
    subplot('Position', FIGPOS(ii,:))
    rho1deg = flipud(reshape(rho(:,ii), [NUMLONS NUMLATS])');
    % zero density is not physical
    rho1deg(rho1deg <= 0) = NaN;
    im = imagesc(LONLIM, LATLIM, rho1deg);
    set(im, 'AlphaData', ~isnan(rho1deg))
    cb = colorbar;
    colormap('kelicol')
    set(cb, 'TickDirection', 'out')
    set(get(cb, 'Label'), 'String', 'density (g/cm^3)')
    grid on
    axis xy
    axis tight
    axis equal
    title(sprintf('layer %d: %s', ii, LAYERNAMES{ii}))
    set(gca, 'TickDir', 'out', 'XTick', LONLIM(1):60:LONLIM(2), ...
        'YTick', LATLIM(1):30:LATLIM(2), 'FontSize', 10)
end
set(gcf, 'Renderer', 'painters')
figdisp([mfilename '_rho'], [], [], 2, [], 'epstopdf')

%% VP
% read the data files
fid = fopen(fullfile(crust1dir, 'crust1.vp'));
vp = fscanf(fid, '%f', [NUMLAYERS Inf]);
vp = vp';
fclose(fid);

% make a plot
figure(2)
set(gcf, 'Units', 'inches', 'Position', [0 1 12 6])
clf
for ii = 1:NUMLAYERS
    subplot('Position', FIGPOS(ii,:))
    vp1deg = flipud(reshape(vp(:,ii), [NUMLONS NUMLATS])');
    vp1deg(vp1deg <= 0) = NaN;
    im = imagesc(LONLIM, LATLIM, vp1deg);
    set(im, 'AlphaData', ~isnan(vp1deg))
    cb = colorbar;
    colormap('kelicol')
    set(cb, 'TickDirection', 'out')
    set(get(cb, 'Label'), 'String', 'v_p (km/s)')
    grid on
    axis xy
    axis tight
    axis equal
    title(sprintf('layer %d: %s', ii, LAYERNAMES{ii}))
    set(gca, 'TickDir', 'out', 'XTick', LONLIM(1):60:LONLIM(2), ...
        'YTick', LATLIM(1):30:LATLIM(2), 'FontSize', 10)
end
set(gcf, 'Renderer', 'painters')
figdisp([mfilename '_vp'], [], [], 2, [], 'epstopdf')

%% VS
% read the data files
fid = fopen(fullfile(crust1dir, 'crust1.vs'));
vs = fscanf(fid, '%f', [NUMLAYERS Inf]);
vs = vs';
fclose(fid);

% make a plot
figure(3)
set(gcf, 'Units', 'inches', 'Position', [0 1 12 6])
clf
for ii = 1:NUMLAYERS
    subplot('Position', FIGPOS(ii,:))
    vs1deg = flipud(reshape(vs(:,ii), [NUMLONS NUMLATS])');
    if ii > 1
        vs1deg(vs1deg <= 0) = NaN;
    end
    im = imagesc(LONLIM, LATLIM, vs1deg);
    set(im, 'AlphaData', ~isnan(vs1deg))
    cb = colorbar;
    colormap('kelicol')
    set(cb, 'TickDirection', 'out')
    set(get(cb, 'Label'), 'String', 'v_s (km/s)')
    grid on
    axis xy
    axis tight
    axis equal
    title(sprintf('layer %d: %s', ii, LAYERNAMES{ii}))
    set(gca, 'TickDir', 'out', 'XTick', LONLIM(1):60:LONLIM(2), ...
        'YTick', LATLIM(1):30:LATLIM(2), 'FontSize', 10)
end
set(gcf, 'Renderer', 'painters')
figdisp([mfilename '_vs'], [], [], 2, [], 'epstopdf')

%% BNDS: interface elevation
% read the data files
fid = fopen(fullfile(crust1dir, 'crust1.bnds'));
bnds = fscanf(fid, '%f', [NUMLAYERS Inf]);
bnds = bnds';
fclose(fid);

% make a plot
figure(4)
set(gcf, 'Units', 'inches', 'Position', [0 1 12 6])
clf
for ii = 1:NUMLAYERS
    subplot('Position', FIGPOS(ii,:))
    imagesc(LONLIM, LATLIM, flipud(reshape(bnds(:,ii), [NUMLONS NUMLATS])'));
    cb = colorbar;
    colormap('kelicol')
    set(cb, 'TickDirection', 'out')
    set(get(cb, 'Label'), 'String', 'elevation (km)')
    grid on
    axis xy
    axis tight
    axis equal
    title(sprintf('interface %d: %s', ii, BNDNAMES{ii}))
    set(gca, 'TickDir', 'out', 'XTick', LONLIM(1):60:LONLIM(2), ...
        'YTick', LATLIM(1):30:LATLIM(2), 'FontSize', 10)
end
set(gcf, 'Renderer', 'painters')
figdisp([mfilename '_bnds'], [], [], 2, [], 'epstopdf')

%% THCK: layer thickness
% make a plot
figure(5)
set(gcf, 'Units', 'inches', 'Position', [0 1 12 6])
clf
for ii = 1:NUMLAYERS
    subplot('Position', FIGPOS(ii,:))
    if ii < NUMLAYERS
        bnd_top = flipud(reshape(bnds(:,ii), [NUMLONS NUMLATS])');
        bnd_bot = flipud(reshape(bnds(:,ii+1), [NUMLONS NUMLATS])');
    else
        bnd_top = flipud(reshape(bnds(:,1), [NUMLONS NUMLATS])');
        bnd_bot = flipud(reshape(bnds(:,NUMLAYERS), [NUMLONS NUMLATS])');
    end
    th = bnd_top-bnd_bot;
    % make zero thickness NaN
    th = th .* (th ./ th);

    im = imagesc(LONLIM, LATLIM, th);
    cb = colorbar;
    colormap('kelicol')
    set(im, 'AlphaData', ~isnan(th))
    set(cb, 'TickDirection', 'out')
    set(get(cb, 'Label'), 'String', 'thickness (km)')
    grid on
    axis xy
    axis tight
    axis equal
    if ii < NUMLAYERS
        title(sprintf('layer %d: %s', ii, LAYERNAMES{ii}))
    else
        title(sprintf('layer %d to %d: surface to Moho depth', 1, NUMLAYERS-1))
    end
    set(gca, 'TickDir', 'out', 'XTick', LONLIM(1):60:LONLIM(2), ...
        'YTick', LATLIM(1):30:LATLIM(2), 'FontSize', 10)
end
set(gcf, 'Renderer', 'painters')
figdisp([mfilename '_thck'], [], [], 2, [], 'epstopdf')

%% Extra
% 1. global sediment thickness
% 2. ocean sediment thickness
% 3. sediment+crustal thickness
% 4. crustal thickness
% 5. compensation depth pressure (assumes homogeneous Earth gravity)
% 6. 
end