function drawmesh3d(ddir, direction, ii_plot, scaling)
% DRAWMESH3D(ddir, direction, ii_plot, scaling)
%
% Draw the cross-section of a SPECFEM3D simulation
%
% INPUT:
% ddir          directory to a FK-SPECFEM3D run
% direction     'longitudinal' (west->east) or 'transverse' (south->north)
% ii_plot       index of the grid of the cross section
%               [default: [] which is the middle of the grid
% scaling       the number to divide plotted values (e.g. 1 for meters,
%               1000 for kilometers) [default: 1]
%
% SEE ALSO:
% DRAWBACKGROUND (for SPECFEM2D simulation)
%
% Last modified by spipatprathanporn@ucsd.edu, 08/30/2026

defval('direction', 'longitudinal')
defval('ii_plot', [])
defval('scaling', 1)

% displayed colors in the drawing: ocean, sediment, crust, mantle
COLORS = csscolor({'cornflowerblue', 'tan', 'darkgray', 'salmon'});

% data file names
fkmodelfile = fullfile(ddir, 'DATA', 'FKMODEL');
meshparfile = fullfile(ddir, 'DATA', 'meshfem3D_files', 'Mesh_Par_file');
interffile = fullfile(ddir, 'DATA', 'meshfem3D_files', 'interfaces.dat');

% read the data files
fkmodel = loadfkmodel(fkmodelfile);
meshparams = loadmeshparfile3d(meshparfile);
[itfs, layers] = loadinterfacefiles3d(interffile);

% determine box dimension
height = meshparams.DEPTH_BLOCK_KM * 1000;

itfs_str = array2struct(flipud(itfs));
% resamples the mesh of the box bottom to match the number of elemements in
% each direction
itfs_str.Z{1} = zeros(meshparams.NEX_ETA+1, meshparams.NEX_XI+1);
for ii = 2:length(itfs_str.Z)
    if any(size(itfs_str.Z{ii}) ~= [meshparams.NEX_ETA+1 meshparams.NEX_XI+1])
        % interpolate
        [xxold, yyold] = meshgrid(...
            itfs_str.LON_MIN(ii) + (0:itfs_str.NXI(ii)-1) * itfs_str.SPACING_XI(ii), ...
            itfs_str.LAT_MIN(ii) + (0:itfs_str.NETA(ii)-1) * itfs_str.SPACING_ETA(ii) ...
            );
        [xxnew, yynew] = meshgrid(...
            linspace(meshparams.LONGITUDE_MIN, meshparams.LONGITUDE_MAX, meshparams.NEX_XI+1), ...
            linspace(meshparams.LATITUDE_MIN, meshparams.LATITUDE_MAX, meshparams.NEX_ETA+1) ...
            );
        itfs_str.Z{ii} = interp2(xxold, yyold, itfs_str.Z{ii}, ...
            xxnew, yynew, 'linear');
    end
end
itfs_str.ZBOT = [itfs_str.Z(2:end); ones(meshparams.NEX_ETA+1, meshparams.NEX_XI+1) * ...
        -height];
itfs_str.NZ = flipud(layers);

% apply scaling
height = height / scaling;
for ii = 1:length(itfs_str.NZ)
    itfs_str.Z{ii} = itfs_str.Z{ii} / scaling;
    itfs_str.ZBOT{ii} = itfs_str.ZBOT{ii} / scaling;
end

figure
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 9 5])
subplot('Position', [0.08 0.08 0.39 0.84])
imagesc([meshparams.LONGITUDE_MIN meshparams.LONGITUDE_MAX] / scaling, ...
    [meshparams.LATITUDE_MIN meshparams.LATITUDE_MAX] /scaling, ...
    itfs_str.ZBOT{1})
axis xy
axis tight
axis equal
cb = colorbar(gca);
set(cb, 'TickDirection', 'out')
if scaling == 1
    set(get(cb, 'Label'), 'String', 'elevation (m)')
elseif scaling == 1000
    set(get(cb, 'Label'), 'String', 'elevation (km)')
else
    set(get(cb, 'Label'), 'String', sprintf('elevation (x%g m)', scaling))
end
colormap(kelicol)
grid on
hold on
% draw the cross-ection line
if ~strcmpi(direction, 'longitudinal')
    if isempty(ii_plot)
        ii_plot = ceil((meshparams.NEX_XI+1)/2);
    end
    x = (meshparams.LONGITUDE_MIN + (ii_plot-1) / meshparams.NEX_XI * ...
        (meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN)) * [1 1];
    x = x / scaling;
    y = [meshparams.LATITUDE_MIN meshparams.LATITUDE_MAX];
    y = y / scaling;
    plot(x, y, 'Color', 'k', 'LineWidth', 1)
else
    if isempty(ii_plot)
        ii_plot = ceil((meshparams.NEX_ETA+1)/2);
    end
    x = [meshparams.LONGITUDE_MIN meshparams.LONGITUDE_MAX];
    x = x / scaling;
    y = (meshparams.LATITUDE_MIN + (ii_plot-1) / meshparams.NEX_ETA * ...
        (meshparams.LATITUDE_MAX - meshparams.LATITUDE_MIN)) * [1 1];
    y = y / scaling;
    plot(x, y, 'Color', 'k', 'LineWidth', 1)
end
if scaling == 1
    xlabel('longitudinal (m)')
    ylabel('transverse (m)')
elseif scaling == 1000
    xlabel('longitudinal (km)')
    ylabel('transverse (km)')
else
    xlabel(sprintf('longitudinal (x%g m)', scaling))
    ylabel(sprintf('transverse (x%g m)', scaling))
end
title('elevation map')
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12)

subplot('Position', [0.58 0.08 0.38 0.84])
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: roatate the mesh horizontally so that the longitudinal direction 
% points to the right
if mod(fkmodel.baz, 360) == 0
    % rotate 90 degrees counter clockwise

elseif mod(fkmodel.baz, 360) == 90
    % rotate 180 degrees counter clockwise


elseif mod(fkmodel.baz, 360) == 180
    % rotate 270 degrees counter clockwise

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmpi(direction, 'longitudinal')
    x = linspace(meshparams.LATITUDE_MIN, meshparams.LATITUDE_MAX, ...
        meshparams.NEX_ETA+1)' / scaling;
    if isempty(ii_plot)
        ii_plot = ceil((meshparams.NEX_XI+1)/2);
    end

    % draw the layer box
    pgons = {};
    for ii = 1:length(itfs_str.NZ)
        pgons{ii} = polyshape([x; flipud(x)], [itfs_str.Z{ii}(:,ii_plot); ...
            flipud(itfs_str.ZBOT{ii}(:,ii_plot))]);
        plot(pgons{ii}, 'FaceColor', COLORS(ii,:), 'LineWidth', 1, ...
            'EdgeColor', 'k', 'FaceAlpha', 1)
    end

    % draw the vertical grid lines
    for ii = 2:meshparams.NEX_ETA
        plot([1 1] * x(ii), [0 -height], 'LineWidth', 0.25, ...
            'Color', 'k')
    end

    % draw the horiztonal grid lines
    for ii = 1:length(itfs_str.NZ)
        for jj = 1:itfs_str.NZ(ii)-1
            plot(x, itfs_str.ZBOT{ii}(:,ii_plot) + (jj/itfs_str.NZ(ii)) * ...
                (itfs_str.Z{ii}(:,ii_plot) - itfs_str.ZBOT{ii}(:,ii_plot)), ...
                'LineWidth', 0.25, 'Color', 'k')
        end
    end

    % redraw the box outline
    for ii = 1:length(itfs_str.NZ)
        plot(pgons{ii}, 'FaceColor', 'none', 'LineWidth', 1, ...
                'EdgeColor', 'k', 'FaceAlpha', 1)
    end

    if scaling == 1
        xlabel('transverse (m)')
        ylabel('elevation (m)')
    elseif scaling == 1000
        xlabel('transverse (km)')
        ylabel('elevation (km)')
    else
        xlabel(sprintf('transverse (x%g m)', scaling))
        ylabel(sprintf('elevation (x%g m)', scaling))
    end
else
    x = linspace(meshparams.LONGITUDE_MIN, meshparams.LONGITUDE_MAX, ...
        meshparams.NEX_XI+1) / scaling;
    if isempty(ii_plot)
        ii_plot = ceil((meshparams.NEX_ETA+1)/2);
    end

    % draw the layer box
    pgons = {};
    for ii = 1:length(itfs_str.NZ)
        pgons{ii} = polyshape([x fliplr(x)], [itfs_str.Z{ii}(ii_plot,:) ...
            fliplr(itfs_str.ZBOT{ii}(ii_plot,:))]);
        plot(pgons{ii}, 'FaceColor', COLORS(ii,:), 'LineWidth', 1, ...
            'EdgeColor', 'k', 'FaceAlpha', 1)
    end

    % draw the vertical grid lines
    for ii = 2:meshparams.NEX_XI
        plot([1 1] * x(ii), [0 -height], 'LineWidth', 0.25, ...
            'Color', 'k')
    end

    % draw the horiztonal grid lines
    for ii = 1:length(itfs_str.NZ)
        for jj = 1:itfs_str.NZ(ii)-1
            plot(x, itfs_str.ZBOT{ii}(ii_plot,:) + (jj/itfs_str.NZ(ii)) * ...
                (itfs_str.Z{ii}(ii_plot,:) - itfs_str.ZBOT{ii}(ii_plot,:)), ...
                'LineWidth', 0.25, 'Color', 'k')
        end
    end

    % redraw the box outline
    for ii = 1:length(itfs_str.NZ)
        plot(pgons{ii}, 'FaceColor', 'none', 'LineWidth', 1, ...
                'EdgeColor', 'k', 'FaceAlpha', 1)
    end

    if scaling == 1
        xlabel('longitudinal (m)')
        ylabel('elevation (m)')
    elseif scaling == 1000
        xlabel('longitudinal (km)')
        ylabel('elevation (km)')
    else
        xlabel(sprintf('longitudinal (x%g m)', scaling))
        ylabel(sprintf('elevation (x%g m)', scaling))
    end
end
title('cross section')
set(gca, 'TickDir', 'out', 'Box', 'on', 'FontSize', 12, ...
    'XLim', [x(1) x(end)], 'YLim', [-height 0], 'DataAspectRatio', [1 1 1])
set(gcf, 'Renderer', 'painters')
end