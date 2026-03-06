function varargout = taperbathymetry(z, r)
% [ztapered, zflat] = taperbathymetry(z, r)
%
% Tapers the bathymetry grid to the median elevation of the edge of
% untapered region. It applies 2D cosine taper (Hanning window) to the
% deviatory bathymetry.
%
% INPUT:
% z             bathymetry grid
% r             the fraction of the window that is tapered
%
% OUTPUT:
% ztapered      tapered bathymetry
% zflat         bathymetry where the soon-to-be-tapered region is flattened
%               to edge of the non-tapered region
%
% EXAMPLES:
% taperbathymetry('demo', 1)
% taperbathymetry('demo', 2)
% taperbathymetry('demo', 3)
%
% Last modified by sirawich@princeton.edu, 03/06/2026


if strcmp(z, 'demo')
    defval('r', 1)
    savename = sprintf('%s_demo%d', mfilename, r);

    % color axis limit for the gradient magnitude plot
    if r == 1
        topo_clim = [-800 1000];
        grad_clim = [0 0.25];
        s2 = 100000;
        nu = 1.5;
        rho = 4000;
        dx = 250;
    elseif r == 2
        topo_clim = [-1200 1200];
        grad_clim = [0 1];
        s2 = 100000;
        nu = 1;
        rho = 1000;
        dx = 250;
    elseif r == 3
        topo_clim = [0 1200];
        grad_clim = [0 0.25];
        s2 = 100000;
        nu = 1;
        rho = 8000;
        dx = 100;
    end

    % read a topography file
    fname = fullfile(getenv('MERMAID3'), 'examples', ...
        sprintf('taperbathymetry_topofile%d.dat', r));
    fid = fopen(fname, 'r');
    z = fscanf(fid, '%f');
    fclose(fid);
    z = reshape(z, 20000/dx + 1, 20000/dx + 1);
    r = 0.25;

    % title
    figure(2)
    set(gcf, 'Units', 'inches', 'Position', [0 1 11.5 6])
    clf
    subplot('Position', [0.13 0.92 0.74 0.001])
    titlestring = sprintf(['Taperbathymetry demo: \\sigma^2 = %d' ...
                           ' | \\nu = %.2f | \\rho = %d | ' ...
                           '\\Deltax = %d'], s2, nu, rho, dx);
    title(titlestring)
    nolabels(gca, 3)
    set(gca, 'Box', 'off', 'CLim', [-800 1000], 'FontSize', 12, ...
        'Color', 'none')
    set(getfield(gca, 'XAxis'), 'Visible', 'off')
    set(getfield(gca, 'YAxis'), 'Visible', 'off')

    % original, no taper
    axlim = [-10 10]; % in km
    subplot('Position', [0.05 0.51 0.19 0.35])
    imagesc(axlim, axlim, z')
    colormap(gca, "kelicol")
    grid on
    hold on
    plot([-5 5 5 -5 -5], [-5 -5 5 5 -5], 'Color', [0.3 0.3 0.3]);
    title('Original, no taper')
    ylabel('northing (km)')
    set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11, ...
        'CLim', topo_clim)
    nolabels(gca, 1);

    % taper to the edge
    zedge = [z(:,1); z(:,end); z(1,:)'; z(end,:)'];
    z0 = median(zedge);
    w = shanning2d(size(z), r);
    subplot('Position', [0.27 0.51 0.19 0.35])
    imagesc(axlim, axlim, w .* (z' - z0) + z0)
    colormap(gca, "kelicol")
    grid on
    hold on
    plot([-5 5 5 -5 -5], [-5 -5 5 5 -5], 'Color', [0.3 0.3 0.3]);
    title('Cosine taper to the outer edge')
    set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11, ...
        'CLim', topo_clim)
    nolabels(gca, 3);
    cb = colorbar('Location', 'eastoutside');
    set(gca, 'Position', [0.27 0.51 0.19 0.35])
    set(cb, 'FontSize', 8)
    cb.Position(3) = 0.01;
    set(cb.('Label'), 'String', 'elevation (m)', 'FontSize', 10)

    % gradient magnitude of the tapered bathymetry
    [gx0, gy0] = gradient(w .* (z - z0) + z0, dx, dx);
    G0 = sqrt(gx0.^2 + gy0.^2);
    subplot('Position', [0.54 0.51 0.19 0.35])
    imagesc(axlim, axlim, G0')
    colormap(gca, "turbo")
    grid on
    title('gradient magnitude')
    set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11, ...
        'CLim', grad_clim)
    nolabels(gca, 3);
    cb = colorbar('Location', 'eastoutside');
    set(gca, 'Position', [0.54 0.51 0.19 0.35])
    set(cb, 'FontSize', 8)
    cb.Position(3) = 0.01;
    set(cb.('Label'), 'String', 'gradient magnitude', 'FontSize', 10)

    % flatten the outer edge
    [ztapered, zflat] = taperbathymetry(z, r);
    subplot('Position', [0.05 0.08 0.19 0.35])
    imagesc(axlim, axlim, zflat')
    colormap(gca, "kelicol")
    grid on
    hold on
    plot([-5 5 5 -5 -5], [-5 -5 5 5 -5], 'Color', [0.3 0.3 0.3]);
    title('Flattened, no taper')
    xlabel('easting (km)')
    ylabel('northing (km)')
    set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11, ...
        'CLim', topo_clim)

    % taper to the inner boundary
    subplot('Position', [0.27 0.08 0.19 0.35])
    imagesc(axlim, axlim, ztapered')
    colormap(gca, "kelicol")
    grid on
    hold on
    plot([-5 5 5 -5 -5], [-5 -5 5 5 -5], 'Color', [0.3 0.3 0.3]);
    title('Flatten, taper to the inner edge')
    xlabel('easting (km)')
    set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11, ...
        'CLim', topo_clim)
    nolabels(gca, 2);

    % absolute gradient
    [gx, gy] = gradient(ztapered, dx, dx);
    G = sqrt(gx.^2 + gy.^2);
    subplot('Position', [0.54 0.08 0.19 0.35])
    imagesc(axlim, axlim, G')
    colormap(gca, "turbo")
    grid on
    title('gradient magnitude')
    xlabel('easting (km)')
    set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11, ...
        'CLim', grad_clim)
    nolabels(gca, 2);

    % plot where gradient magnitude is greater than that of other method
    subplot('Position', [0.80 0.51 0.19 0.35])
    cm = colormap(gca, 'turbo');
    cm(1,:) = [0 0 0];
    threshold = grad_clim(2)/size(cm,1);
    imagesc(axlim, axlim, G0' .* (G0'-G' > threshold));
    colormap(gca, cm)
    grid on
    title('worse gradient magnitude')
    set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11, ...
        'CLim', grad_clim)
    nolabels(gca, 3);

    subplot('Position', [0.80 0.08 0.19 0.35])
    imagesc(axlim, axlim, G' .* (G'-G0' > threshold));
    colormap(gca, cm)
    grid on
    title('worse gradient magnitude')
    xlabel('easting (km)')
    set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11, ...
        'CLim', grad_clim)
    nolabels(gca, 2);
    
    set(gcf, 'Renderer', 'painter')
    figdisp(savename, [], [], 2, [], 'epstopdf');
    return
end

defval('r', 0.5)
if r > 0.5
    error('r must be between 0 and 0.5.')
end

% size of the box
Nx = size(z, 1);
Ny = size(z, 2);

% width of tapered region
Wtx = round(r * Nx);
Wty = round(r * Ny);

% determine the median of the edge of the non-tapered region
% the corner points are not double counted
zedge = [z(1+Wtx, 1+Wty:Ny-Wty)'; ...
         z(Nx-Wtx, 1+Wty:Ny-Wty)'; ...
         z(2+Nx:Nx-Wtx-1, 1+Wty); ...
         z(2+Nx:Nx-Wtx-1, Ny-Wty)];
z0 = median(zedge);

% flatten the topography in the tapered region
zflat = z;
for ii = 1:Wtx
    zflat(ii,:) = zflat(Wtx+1,:);
end
for ii = Nx-Wtx+1:Nx
    zflat(ii,:) = zflat(Nx-Wtx,:);
end
for ii = 1:Wty
    zflat(:,ii) = zflat(:,Wty+1);
end
for ii = Ny-Wty+1:Ny
    zflat(:,ii) = zflat(:,Ny-Wty);
end
    
% 2D cosine taper
w = shanning2d([Nx Ny], r);

% tapering
ztapered = w .* (zflat - z0) + z0;

varns = {ztapered, zflat};
varargout = varns(1:nargout);
end