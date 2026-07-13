function [itfs, layers] = loadinterfacefiles3d_demo
% [itfs, layers] = LOADINTERFACEFILES3D_DEMO
% Demo of LOADINTERFACEFILES3D
%
% It reads the LA basin example and plot the surface elevations. This will
% give you sense of how to read/write topography files and how to orient
% yourself.
%
% The SPECFEM3D example is available at
% https://github.com/SPECFEM/specfem3d/tree/master/EXAMPLES/applications/meshfem3D_examples/simple_model
%
% Read https://specfem3d.readthedocs.io/en/latest/03_mesh_generation on how
% the interfaces are stored in files
%
% SEE ALSO:
% LOADINTERFACEFILES3D, WRITEINTERFACEFILES3D
%
% Last modified by spipatprathanporn@ucsd.edu, 07/13/2026

% Local directory of the example
% Please set $SPECFEM3D to be where you clone your SPECFEM3D.
ddir = fullfile(getenv('SPECFEM3D'), 'EXAMPLES', 'applications', ...
    'meshfem3D_examples', 'simple_model');

% Reads the interface files
[itfs, layers] = loadinterfacefiles3d(fullfile(ddir, 'DATA', ...
    'meshfem3D_files', 'interfaces.dat'));

for ii = 1:2
    % plot interface
    x = itfs{ii}.LON_MIN + (0:itfs{ii}.NXI-1) * itfs{ii}.SPACING_XI;
    y = itfs{ii}.LAT_MIN + (0:itfs{ii}.NETA-1) * itfs{ii}.SPACING_ETA;
    
    figure(ii+10)
    clf
    set(gcf, 'Unit', 'inches', 'Position', [0 1 6 9])
    subplot(311)
    imagesc(x, y, itfs{ii}.Z)
    axis xy
    axis tight
    axis equal
    % add the colorbar
    if ii == 1        
        cb = colorbar;
        colormap('parula')
    else
        cb = cax2dem([-6000 4000]);
    end
    cm = colormap(gca);
    set(get(cb, 'Label'), 'String', 'elevation (m)')
    set(cb, 'TickDirection', 'out')
    grid on
    if ii == 1
        xlabel('easting (m)')
        ylabel('northing (m)')
        title('subsurface elevation: SUPPRESS\_UTM\_PROJECTION == .true.')
    else
        xlabel('lonitude')
        ylabel('latitude')
        title('Los Angeles elevation: SUPPRESS\_UTM\_PROJECTION == .false.')
    end
    subtitle(sprintf('imagesc(x, y, itfs\\{%d\\}.Z); axis xy', ii))
    
    subplot(312)
    imagesc(itfs{ii}.Z)
    axis tight
    axis equal
    cb = colorbar;
    colormap(gca, 'parula')
    set(get(cb, 'Label'), 'String', 'Z(i,j) = Z(i\_ETA, i\_XI)')
    set(cb, 'TickDirection', 'out')
    grid on
    xlabel('j: XI index')
    ylabel('i: ETA index')
    title(sprintf('Elevation grid (Z): NXI = %d, NETA = %d', itfs{ii}.NXI, ...
        itfs{ii}.NETA))
    subtitle(sprintf('imagesc(itfs\\{%d\\}.Z)', ii))
    
    subplot(313)
    [xx, yy] = meshgrid(x, y);
    surface(xx, yy, itfs{ii}.Z, 'EdgeColor', 'none')
    colormap(gca, cm)
    if ii == 2
        clim([-6000 4000])
    end
    cb = colorbar;
    set(get(cb, 'Label'), 'String', 'elevation (m)')
    set(cb, 'TickDirection', 'out')
    grid on
    box on
    if ii == 1
        xlabel('easting (m)')
        ylabel('northing (m)')
    else
        xlabel('longitude')
        ylabel('latitude')
    end
    title('Surface plot')
    subtitle(sprintf('[xx,yy] = meshgrid(x,y); surface(xx,yy,itfs\\{%d\\}.Z, "EdgeColor", "none"); ', ii))
    if ii == 1
        set(gca, 'CameraPosition', [-0.5055 -1.7020 0.0179] * 1e7, ...
            'CameraTarget', [2000000 5085000 -30000])
    else
        set(gca, 'CameraPosition', [-128.4078 -2.8601 3.3444e+04], ...
            'CameraTarget', [-118 34.5000 0])
    end
end
end