function [itfs, layers] = loadinterfacefiles3d(fname, plt)
% [itfs, layers] = LOADINTERFACEFILES3D(fname, plt)
%
% Reads interfaces from an interface file and other accompanied files for a
% SPECFEM3D_Cartesian simulation. Please make sure that other accompanied
% files are in the same directory as the interface file.
%
% Read https://specfem3d.readthedocs.io/en/latest/03_mesh_generation on how
% the interfaces are stored in files
%
% INPUT:
% fname         name of the interface file
% plt           whether to show examples of plots (how to plot)
%               0   - no plot [default]
%               n>0 - plot n-th interfaces from the top
%
% OUTPUT:
% itfs          interfaces, an array of struct with following fields
%       SUPPRESS_UTM_PROJECTION     whether to suppress UTM projection
%       NXI                         number of elements in x-direction
%       NETA                        number of elements in y-direction
%       LON_MIN                     minimum longitude (or x value)
%       LAT_MIN                     minimum latitude  (or y value)
%       SPACING_XI                  spacing in x-direction
%       SPACING_ETA                 spacing in y-direction
%       FILE                        elevation file name
%       Z                           elevation grid at (Y,X) or (LAT,LON)
% layers        number of vertical spectral elements for each layer
%
% SEE ALSO:
% MAKEINTERFACES3D, WRITEINTERFACEFILES3D, LOADINTERFACEFILES3D_DEMO
%
% Last modified by sirawich-at-princeton.edu, 07/13/2026

defval('plt', 0)

ddir = strcat(fileparts(fname), filesep);

%% open the file
fid = fopen(fname, 'r');
line = strip(fgetl(fid));
% skip comments / headers
while isempty(line) || strcmp(line(1), '#')
    line = strip(fgetl(fid));
end

%% read number of interfaces
numinterfaces = sscanf(line, '%d', 1);

itfs = cell(numinterfaces, 1);
layers = nan(numinterfaces, 1);

%% read the interfaces
for ii = 1:numinterfaces
    line = strip(fgetl(fid));
    % skip comments / headers
    while isempty(line) || strcmp(line(1), '#')
        line = strip(fgetl(fid));
    end
    [boolword, ~, ~, ni] = sscanf(line, '%s', 1);
    if strcmp(boolword, '.true.')
        itfs{ii}.SUPPRESS_UTM_PROJECTION = true;
    elseif strcmp(boolword, '.false.')
        itfs{ii}.SUPPRESS_UTM_PROJECTION = false;
    else
        warning(['SUPPRESS_UTM_PROJECTION flag could not be ' ...
            'interpreted, assumed to be true then.'])
        itfs{ii}.SUPPRESS_UTM_PROJECTION = true;
    end
    line = replace(line, 'd', 'e');
    nums = sscanf(line(ni:end), '%g');
    itfs{ii}.NXI = nums(1);
    itfs{ii}.NETA = nums(2);
    itfs{ii}.LON_MIN = nums(3);
    itfs{ii}.LAT_MIN = nums(4);
    itfs{ii}.SPACING_XI = nums(5);
    itfs{ii}.SPACING_ETA = nums(6);
    
    line = strip(fgetl(fid));
    itfs{ii}.FILE = strcat(ddir, sscanf(line, '%s', 1));
    % try to read the elevation from the file
    try
        fid_ii = fopen(itfs{ii}.FILE, 'r');
        z = fscanf(fid_ii, '%f');
        if itfs{ii}.SUPPRESS_UTM_PROJECTION
            itfs{ii}.Z = reshape(z, itfs{ii}.NETA, itfs{ii}.NXI);
        else
            itfs{ii}.Z = reshape(z, itfs{ii}.NXI, itfs{ii}.NETA)';
        end
        fclose(fid_ii);
    catch ME
        fprintf('Encounter problems while reading the file %s\n', ...
            itfs{ii}.FILE);
        getReport(ME)
        continue
    end
end

%% read the number of spectral elements in the vertical direction
for ii = 1:numinterfaces
    line = strip(fgetl(fid));
    % skip comments / headers
    while isempty(line) || strcmp(line(1), '#')
        line = strip(fgetl(fid));
    end
    layers(ii) = sscanf(line, '%d', 1);
end

%% close the file
fclose(fid);

%% plot
if plt
    n = length(itfs);
    i_plt = n + 1 - plt;
    x = itfs{2}.LON_MIN + (0:itfs{i_plt}.NXI-1) * itfs{i_plt}.SPACING_XI;
    y = itfs{2}.LAT_MIN + (0:itfs{i_plt}.NETA-1) * itfs{i_plt}.SPACING_ETA;
    x = x/1000;
    y = y/1000;

    figure(12)
    clf
    set(gcf, 'Unit', 'inches', 'Position', [0 1 6 9])
    subplot(311)
    imagesc(x, y, itfs{i_plt}.Z)
    axis xy
    axis tight
    axis equal
    cb = colorbar;
    colormap(kelicol)
    set(get(cb, 'Label'), 'String', 'elevation (m)')
    set(cb, 'TickDirection', 'out')
    grid on
    xlabel('easting (km)')
    ylabel('northing (km)')
    title('Image plot')
    subtitle(sprintf('imagesc(x, y, itfs\\{%d\\}.Z); axis xy', i_plt))
    
    subplot(312)
    imagesc(itfs{2}.Z)
    axis tight
    axis equal
    cb = colorbar;
    colormap(kelicol)
    set(get(cb, 'Label'), 'String', 'Z(i,j) = Z(i\_ETA, i\_XI)')
    set(cb, 'TickDirection', 'out')
    grid on
    xlabel('j: XI index')
    ylabel('i: ETA index')
    title(sprintf('Elevation grid (Z): NXI = %d, NETA = %d', itfs{i_plt}.NXI, ...
        itfs{i_plt}.NETA))
    subtitle(sprintf('imagesc(itfs\\{%d\\})', i_plt))
    
    subplot(313)
    [xx, yy] = meshgrid(x, y);
    surface(xx, yy, itfs{i_plt}.Z, 'EdgeColor', 'none')
    cb = colorbar;
    colormap(kelicol)
    set(get(cb, 'Label'), 'String', 'elevation (m)')
    set(cb, 'TickDirection', 'out')
    grid on
    box on
    xlabel('easting (km)')
    ylabel('northing (km)')
    zlabel('elevation (m)')
    title('Surface plot')
    subtitle(sprintf('[xx,yy] = meshgrid(x,y); surface(xx,yy,itfs\\{%d\\}.Z, "EdgeColor", "none");', i_plt))

end
end