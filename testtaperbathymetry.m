function varargout = testtaperbathymetry(s2, nu, rho, dx, big, small, taper, plt)
% [bath_big, bath_small, bath_small_old] = ...
%     TESTTAPERBATHYMETRY(s2, nu, rho, dx, big, small, taper)
%
% Generates the tapered bathymetry grids, both the cosine tapered and the
% modified version, as well as the big bathymetry grids serving as the
% ground truth to test the effects of tapered bahtymetry on waveform
% modling.
%
% INPUT:
% s2            The first Matern parameter, aka sigma^2, in field units^2 
% nu            The second Matern parameter (smoothness)
% rho           The third Matern parameter, in the units of the grid
%               (length scale)
% dx            element size
% big           the number of sample points of the big bathymetry grid
% small         the number of sample points of the small bathymetry grid
% taper         struct with a following member variables
%       trim        width of the 4 edges that will be set to be flat
%       cosine      width of the cosine taper
%
% OUTPUT:
% bath_big          big bathymetry grid
% bath_small        small bathymetry grid (new method, modified 2D taper)
% bath_small_old    small bathymetry grid (old method, 2D cosine taper)
%
% SEE ALSO:
% TAPERBATHYMETRY, SIMULOSL
% 
% Last modified by sirawich-at-princeton.edu, 04/24/2026

defval('s2', 10000)
defval('nu', 0.5)
defval('rho', 2000)
defval('dx', 250)
defval('big', 257)
defval('small', 97)
defval('plt', true)
defstruct('taper', 'trim', 8)
defstruct('taper', 'cosine', 16)

% generate a random field
p.NyNx = [big big];
p.dydx = [dx dx];
p.blurs = Inf;
p.taper = 1;
th0 = [s2 nu rho];
topo_big = v2s(simulosl(th0, p));

% taper the big random field
topo_big_trim = topo_big(1+taper.trim:big-taper.trim, ...
    1+taper.trim:big-taper.trim);
topo_big_trim_taper = taperbathymetry(topo_big_trim, ...
    taper.cosine / (big - 2*taper.trim));
topo_big_taper = ones(size(topo_big)) * topo_big_trim_taper(1,1);
topo_big_taper(1+taper.trim:big-taper.trim, ...
    1+taper.trim:big-taper.trim) = topo_big_trim_taper;

% taper the small random field
halfdiff = (big - small) / 2;
topo_small = topo_big(1+halfdiff:big-halfdiff, 1+halfdiff:big-halfdiff);
topo_small_trim = topo_small(1+taper.trim:small-taper.trim, ...
    1+taper.trim:small-taper.trim);
topo_small_trim_taper = taperbathymetry(topo_small_trim, ...
    taper.cosine / (small - 2*taper.trim));
topo_small_taper = ones(size(topo_small)) * topo_small_trim_taper(1,1);
topo_small_taper(1+taper.trim:small-taper.trim, ...
    1+taper.trim:small-taper.trim) = topo_small_trim_taper;

% taper the small random field in an old fashion
w = shanning2d(size(topo_small_trim), ...
    taper.cosine / (small - 2*taper.trim), 0, 0);
topo_small_trim_taper_old = topo_small_trim .* w;
topo_small_taper_old = ones(size(topo_small)) * ...
    topo_small_trim_taper_old(1,1);
topo_small_taper_old(1+taper.trim:small-taper.trim, ...
    1+taper.trim:small-taper.trim) = topo_small_trim_taper_old;

if plt
    axlim_big = [-1 1] * (big-1)/2 * dx / 1000;
    axlim_small = [-1 1] * (small-1)/2 * dx / 1000;
    figure
    set(gcf, 'Units', 'inches', 'Position', [0 1 10 6])
    
    % title
    subplot('Position', [0.13 0.94 0.74 0.01])
    titlestring = sprintf(['testtaperbathymetry: \\sigma^2 = %d' ...
                           ' | \\nu = %.2f | \\rho = %d | ' ...
                           '\\Deltax = %d'], s2, nu, rho, dx);
    title(titlestring)
    nolabels(gca, 3)
    set(gca, 'Box', 'off', 'CLim', [-800 1000], 'FontSize', 12, ...
        'Color', 'none')
    set(getfield(gca, 'XAxis'), 'Visible', 'off')
    set(getfield(gca, 'YAxis'), 'Visible', 'off')

    subplot('Position', [0.06 0.54 0.27 0.35])
    imagesc(axlim_big, axlim_big, topo_big_taper)
    axis xy
    axis equal
    axis tight
    grid on
    hold on
    plot([-1 1 1 -1 -1] * axlim_small(2), ...
        [-1 -1 1 1 -1] * axlim_small(2), ...
        'LineWidth', 1, 'Color', 'k')
    plot([-1 1 1 -1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        'LineWidth', 0.75, 'Color', 'k', 'LineStyle', ':')
    plot([-1 1 1 -1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        'LineWidth', 1.25, 'Color', 'k', 'LineStyle', ':')
    xlabel('easting (km)')
    ylabel('northing (km)')
    title('big bathymetry grid')
    colormap(gca, 'kelicol')
    clim([-3 3] * sqrt(s2))
    set(gca, 'FontSize', 12, 'Box', 'on', 'TickDir', 'out', ...
        'XTick', get(gca, 'YTick'))
    set(gca, 'YTick', get(gca, 'XTick'))

    subplot('Position', [0.35 0.54 0.27 0.35])
    imagesc(axlim_small, axlim_small, topo_small_taper)
    axis xy
    axis equal
    axis tight
    grid on
    hold on
    plot([-1 1 1 -1 -1] * axlim_small(2), ...
        [-1 -1 1 1 -1] * axlim_small(2), ...
        'LineWidth', 2, 'Color', 'k', 'LineStyle', '-')
    plot([-1 1 1 -1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        'LineWidth', 0.75, 'Color', 'k', 'LineStyle', ':')
    plot([-1 1 1 -1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        'LineWidth', 1.25, 'Color', 'k', 'LineStyle', ':')
    xlabel('easting (km)')
    ylabel('northing (km)')
    title('small bathymetry grid')
    colormap(gca, 'kelicol')
    clim([-3 3] * sqrt(s2))
    set(gca, 'FontSize', 12, 'Box', 'on', 'TickDir', 'out', ...
        'XTick', get(gca, 'YTick'))
    set(gca, 'YTick', get(gca, 'XTick'))

    subplot('Position', [0.64 0.54 0.34 0.35])
    imagesc(axlim_small, axlim_small, topo_small_taper_old)
    axis xy
    axis equal
    axis tight
    grid on
    hold on
    plot([-1 1 1 -1 -1] * axlim_small(2), ...
        [-1 -1 1 1 -1] * axlim_small(2), ...
        'LineWidth', 2, 'Color', 'k', 'LineStyle', '-')
    plot([-1 1 1 -1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        'LineWidth', 0.75, 'Color', 'k', 'LineStyle', ':')
    plot([-1 1 1 -1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        'LineWidth', 1.25, 'Color', 'k', 'LineStyle', ':')
    xlabel('easting (km)')
    ylabel('northing (km)')
    title('small bathymetry grid old')
    colormap(gca, 'kelicol')
    clim([-3 3] * sqrt(s2))
    set(gca, 'FontSize', 12, 'Box', 'on', 'TickDir', 'out', ...
        'XTick', get(gca, 'YTick'))
    set(gca, 'YTick', get(gca, 'XTick'))
    cb = colorbar;
    set(cb, 'FontSize', 12, 'TickDirection', 'both')
    set(get(cb, 'Label'), 'String', 'elevation (m)')

    % compute gradient magnitudes
    [gx_big, gy_big] = gradient(topo_big_taper, dx);
    G_big = sqrt(gx_big.^2 + gy_big.^2);
    [gx_small, gy_small] = gradient(topo_small_taper, dx);
    G_small = sqrt(gx_small.^2 + gy_small.^2);
    [gx_small_old, gy_small_old] = gradient(topo_small_taper_old, dx);
    G_small_old = sqrt(gx_small_old.^2 + gy_small_old.^2);

    maxgradient = max([max(G_big, [], 'all') ...
        max(G_small, [], 'all') ...
        max(G_small_old, [], 'all')]);
    
    subplot('Position', [0.06 0.09 0.27 0.35])
    imagesc(axlim_big, axlim_big, G_big)
    axis xy
    axis equal
    axis tight
    grid on
    hold on
    plot([-1 1 1 -1 -1] * axlim_small(2), ...
        [-1 -1 1 1 -1] * axlim_small(2), ...
        'LineWidth', 2, 'Color', [1 0.75 1], 'LineStyle', '-')
    plot([-1 1 1 -1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        'LineWidth', 0.75, 'Color', [1 0.75 1], 'LineStyle', ':')
    plot([-1 1 1 -1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        'LineWidth', 1.25, 'Color', [1 0.5 1], 'LineStyle', ':')
    xlabel('easting (km)')
    ylabel('northing (km)')
    colormap(gca, 'turbo')
    clim([0 maxgradient])
    set(gca, 'FontSize', 12, 'Box', 'on', 'TickDir', 'out', ...
        'XTick', get(gca, 'YTick'))
    set(gca, 'YTick', get(gca, 'XTick'))

    subplot('Position', [0.35 0.09 0.27 0.35])
    imagesc(axlim_small, axlim_small, G_small)
    axis xy
    axis equal
    axis tight
    grid on
    hold on
    plot([-1 1 1 -1 -1] * axlim_small(2), ...
        [-1 -1 1 1 -1] * axlim_small(2), ...
        'LineWidth', 2, 'Color', [1 0.75 1], 'LineStyle', '-')
    plot([-1 1 1 -1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        'LineWidth', 0.75, 'Color', [1 0.75 1], 'LineStyle', ':')
    plot([-1 1 1 -1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        'LineWidth', 1.25, 'Color', [1 0.5 1], 'LineStyle', ':')
    xlabel('easting (km)')
    ylabel('northing (km)')
    colormap(gca, 'turbo')
    clim([0 maxgradient])
    set(gca, 'FontSize', 12, 'Box', 'on', 'TickDir', 'out', ...
        'XTick', get(gca, 'YTick'))
    set(gca, 'YTick', get(gca, 'XTick'))

    subplot('Position', [0.64 0.09 0.34 0.35])
    imagesc(axlim_small, axlim_small, G_small_old)
    axis xy
    axis equal
    axis tight
    grid on
    hold on
    plot([-1 1 1 -1 -1] * axlim_small(2), ...
        [-1 -1 1 1 -1] * axlim_small(2), ...
        'LineWidth', 2, 'Color', [1 0.75 1], 'LineStyle', '-')
    plot([-1 1 1 -1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*taper.trim -1) / 2 * dx / 1000, ...
        'LineWidth', 0.75, 'Color', [1 0.75 1], 'LineStyle', ':')
    plot([-1 1 1 -1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        [-1 -1 1 1 -1] * (small - 2*(taper.trim+taper.cosine) -1) / 2 * dx / 1000, ...
        'LineWidth', 1.25, 'Color', [1 0.5 1], 'LineStyle', ':')
    xlabel('easting (km)')
    ylabel('northing (km)')
    colormap(gca, 'turbo')
    clim([0 maxgradient])
    set(gca, 'FontSize', 12, 'Box', 'on', 'TickDir', 'out', ...
        'XTick', get(gca, 'YTick'))
    set(gca, 'YTick', get(gca, 'XTick'))
    cb = colorbar;
    set(cb, 'FontSize', 12, 'TickDirection', 'both')
    set(get(cb, 'Label'), 'String', 'gradient magnitude')
end
    

% collect the outputs
varns = {topo_big_taper, topo_small_taper, topo_small_taper_old};
varargout = varns(1:nargout);
end