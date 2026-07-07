function remesh3d(ddir, dx, dz)
% REMESH3D(ddir, dx, [dz1 dz2 dz3 ...])
%
% Updates the element size in each layer. You may leave dx or dz blank if
% you wish not to change them.
%
% INPUTS:
% ddir              directory to a FK-SPECFEM3D run
% dx                element size in x and y directions
% [dz1 dz2 dz3 ...] element size in z direction for each layer from top to
%                   bottom
%
% Last modified by spipatprathanporn@ucsd.edu, 07/06/2026

defval('ddir', [])
defval('dx', [])
defval('dz', [])

% If you leave ddir blank, you do not wish to remesh.
if isempty(ddir)
    return
end

% make sure that dz is a column vector
if size(dz, 2) > 1
    dz = dz';
end

% data file names
fkfile = fullfile(ddir, 'DATA', 'FKMODEL');
parfile = fullfile(ddir, 'DATA', 'Par_file');
meshparfile = fullfile(ddir, 'DATA', 'meshfem3D_files', 'Mesh_Par_file');
interffile = fullfile(ddir, 'DATA', 'meshfem3D_files', 'interfaces.dat');

% read the data files
fkmodel = loadfkmodel(fkfile);
params = loadparfile3d(parfile);
meshparams = loadmeshparfile3d(meshparfile);
[itfs, layers] = loadinterfacefiles3d(interffile);

% determine box dimension
lon_width = meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN;
lat_width = meshparams.LATITUDE_MAX - meshparams.LATITUDE_MIN;
height = meshparams.DEPTH_BLOCK_KM * 1000;

% determine thickness of the layer outside of the box
fklayers = array2struct(fkmodel.layers);
fklayers.zbot = [fklayers.ztop(2:end); -meshparams.DEPTH_BLOCK_KM*1000];
fklayers.thickness = fklayers.ztop - fklayers.zbot;

% list of elastic/acoustic wave speed in the media
vp = fklayers.vp;
vs = fklayers.vs;
vp(vp <= 0) = nan;
vs(vs <= 0) = nan;
cmin = min(vp, vs);
cmax = max(vp, vs);

% determine minimum and maximum recommended element size
dxmin = (7 + sqrt(21)) * max(cmax) * params.DT;
dxmax = 4/5 * min(cmin) / fkmodel.fmax;

dzmin = (7 + sqrt(21)) * cmax * params.DT;
dzmax = 4/5 * cmin / fkmodel.fmax;

% report the minimum and maximum allowed values:
fprintf('Minimum dx = %.3g\n', dxmin)
fprintf('Maximum dx = %.3g\n', dxmax)
fprintf('Minimum dz -- maximum dz in each layer:\n')
for ii = 1:length(dzmin)
    fprintf('Layer %d: %.2f -- %.2f\n', ii, dzmin(ii), dzmax(ii))
end

% verify the input
if dxmin > dxmax
    error(sprintf('%s:impossibleCondition', upper(mfilename)), ...
        ['Impossible condition: the minimum element size is greater ' ...
        'than the maximum element size. Consider decreasing the time ' ...
        'step size (DT in Par_file) or decreasing the maximum resolved ' ...
        'freqeuncy (fmax in FKMODEL).'])
end
if floor(lon_width / dxmin) < ceil(lon_width / dxmax)
    error(sprintf('%s:tooNarrowDXWindow', upper(mfilename)), ...
        ['dxmin and dxmax are so close to each other that there is no ' ...
        'number of elements that satisfy this criterion. Consider ' ...
        'decreasing the time step size (DT in Par_file) or decreasing ' ...
        'the maximum resolved freqeuncy (fmax in FKMODEL).'])
end
if floor(lat_width / dxmin) < ceil(lat_width / dxmax)
    error(sprintf('%s:tooNarrowDXWindow', upper(mfilename)), ...
        ['dxmin and dxmax are so close to each other that there is no ' ...
        'number of elements that satisfy this criterion. Consider ' ...
        'decreasing the time step size (DT in Par_file) or decreasing ' ...
        'the maximum resolved freqeuncy (fmax in FKMODEL).'])
end
if ~isempty(dx) && dx < dxmin
    warning(sprintf('%s:tooSmallDX', upper(mfilename)), ...
        'dx is smaller than recommended element size of %.3g', dxmin)
    dx = dxmin;
end
if ~isempty(dx) && dx > dxmax
    warning(sprintf('%s:tooLargeDX', upper(mfilename)), ...
        'dx is larger than recommended element size of %.3g', dxmax)
    dx = dxmax;
end
if ~isempty(dz) && any(dz < dzmin)
    wh = (dz < dzmin);
    for ii = 1:length(wh)
        if wh(ii)
            warning(sprintf('%s:tooSmallDZ', upper(mfilename)), ...
                ['dz(%d) is smaller than recommended vertical element ' ...
                'size in layer %d of %.3g'], ii, ii, dzmin(ii));
        end
    end
    dz = max(dz, dzmin);
end
if ~isempty(dz) && any(dz > dzmax)
    wh = (dz > dzmax);
    for ii = 1:length(wh)
        if wh(ii)
            warning(sprintf('%s:tooLargeDZ', upper(mfilename)), ...
                ['dz(%d) is larger than recommended vertical element ' ...
                'size in layer %d of %.3g'], ii, ii, dzmax(ii));
        end
    end
    dz = min(dz, dzmax);
end

if ~isempty(dx)
    % determine the number of elements in x (XI/longitude) direction
    NEX_XI = ceil(lon_width / dx);
    if lon_width / NEX_XI < dxmin
        fprintf('Rounding up NXI resulted dx < dxmin. Rounding down\n');
        NEX_XI = NEX_XI - 1;
        if lon_width / NEX_XI > dxmax
            fprintf('Rounding down NEX_XI resulted dx > dxmax.\n');
            error(sprintf('%s:tooNarrowDXWindow', upper(mfilename)), ...
                ['dxmin and dxmax are so close to each other that there ' ...
                'is no number of elements that satisfy this criterion.'])
        end
    end
    
    % determine the number of elements in y (ETA/latitude) direction
    NEX_ETA = ceil(lat_width / dx);
    if lat_width / NEX_ETA < dxmin
        fprintf('Rounding up NXI resulted dx < dxmin. Rounding down\n');
        NEX_ETA = NEX_ETA - 1;
        if lat_width / NEX_ETA > dxmax
            fprintf('Rounding down NEX_ETA resulted dx > dxmax.\n');
            error(sprintf('%s:tooNarrowDXWindow', upper(mfilename)), ...
                ['dxmin and dxmax are so close to each other that there ' ...
                'is no number of elements that satisfy this criterion.'])
        end
    end

    % interpolate the interface grid
    for ii = 1:length(itfs)
        % only interplate if the grid is not a 2x2 matrix
        % It is pointless to interpolate a plane with no curvature.
        if any(size(itfs{ii}.Z) ~= [2 2])
            [xxold, yyold] = meshgrid(...
                linspace(0, lon_width, itfs{ii}.NXI), ...
                linspace(0, lat_width, itfs{ii}.NETA) ...
                );
            [xxnew, yynew] = meshgrid(...
                linspace(0, lon_width, NEX_XI+1), ...
                linspace(0, lat_width, NEX_ETA+1) ...
                );
            % update the interface array
            itfs{ii}.Z = interp2(xxold, yyold, itfs{ii}.Z, xxnew, ...
                yynew, 'linear');
            itfs{ii}.NXI = NEX_XI + 1;
            itfs{ii}.NETA = NEX_ETA + 1;
            itfs{ii}.SPACING_XI = lon_width / NEX_XI;
            itfs{ii}.SPACING_ETA = lat_width / NEX_ETA;
        end
    end
else
    NEX_XI = meshparams.NEX_XI;
    NEX_ETA = meshparams.NEX_ETA;
end

if ~isempty(dz)
    % Reminder: itfs and layers variables are ascending in elevation while
    % fkmodel.layers is descending in elevation. Here I want everything to
    % be descending in elevation.
    itfs_str = array2struct(flipud(itfs));
    for ii = 1:length(itfs_str.Z)
        if all(size(itfs_str.Z{ii}) == [2 2])
            [xxold, yyold] = meshgrid([0 lon_width], [0 lat_width]);
            [xxnew, yynew] = meshgrid(...
                linspace(0, lon_width, NEX_XI+1), ...
                linspace(0, lat_width, NEX_ETA+1) ...
                );
            itfs_str.Z{ii} = interp2(xxold, yyold, ...
                itfs_str.Z{ii}, xxnew, yynew);
        end
    end
    itfs_str.ZBOT = [itfs_str.Z(2:end); ones(NEX_XI+1, NEX_ETA+1) * ...
        -height];
    itfs_str.NZ = nan(size(itfs_str.Z));
    fprintf('Recommended NZ for each layer: min - (use) - max\n')
    for ii = 1:length(itfs_str.Z)
        thickness = itfs_str.Z{ii} - itfs_str.ZBOT{ii};
        % check if the interfaces allow any possible NZ.
        nzmin = ceil(max(thickness, [], 'all') / dzmax(ii));
        nzmax = floor(min(thickness, [], 'all') / dzmin(ii));
        if nzmin > nzmax
            error(sprintf('%s:impossibleThickness', upper(mfilename)), ...
                ['Layer %d has so much thickness variation that ' ...
                'appropriate number of elements does not exist. ' ...
                'Consider modifying the interfaces, decreasing ' ...
                'time step size, and/or maximum resolved frequency'], ii);
        end
        itfs_str.NZ(ii) = max(nzmin, ...
            min(ceil(fklayers.thickness(ii) / dz(ii)), nzmax));
        fprintf('Layer %d: %d - (%d) - %d layers\n', ii, nzmin, ...
            itfs_str.NZ(ii), nzmax);
    end

    % update variable "layers"
    layers = flipud(itfs_str.NZ);
end

% update mesh parameter
meshparams.NEX_XI = NEX_XI;
meshparams.NEX_ETA = NEX_ETA;
for ii = meshparams.NREGIONS:-1:1
    meshparams.REGIONS{ii}.NEX_XI_END = NEX_XI;
    meshparams.REGIONS{ii}.NEX_ETA_END = NEX_ETA;
    if ii < meshparams.NREGIONS
        meshparams.REGIONS{ii}.NZ_BEGIN = ...
            meshparams.REGIONS{ii+1}.NZ_END + 1;
    end
    meshparams.REGIONS{ii}.NZ_END = meshparams.REGIONS{ii}.NZ_BEGIN + ...
        layers(ii);
end

% write the mesh parameter and interface to the file
writemeshparfile3d(meshparams, meshparfile);
writeinterfacefiles3d(itfs, layers, interffile);
end