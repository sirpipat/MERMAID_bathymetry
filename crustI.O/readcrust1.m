function varargout = readcrust1(ddir, varname, layer, lons, lats)
% val = READCRUST1(ddir, varname, layer, lons, lats)
%
% Returns values of a given variable at given longitudes, latitudes, and
% layers.
%
% INPUT:
% ddir              name of the CRUST1.0 directory
%                   [Default: $IFILES/EARTHMODELS/PHYSICAL/crust1.0]
% varname           name of the variable names
%   "rho"               density
%   "vp"                P-wave speed
%   "vs"                S-wave speed
%   "bnds", "z"         boundary surface elevation (<0 for below sealevel)
%   "thick", "dz"       layer thickness
%   "iswater"           is it water
%   "isice"             is it ice
%   "lsbnds"            list boundaries
%   "lslayers"          list layers
% layer             layer/boundary number
% lons              longitude(s) which can be either
%   - []                entire Earth (Default)
%   - single value      single point
%   - two values        [lonmin lonmax] with 1 degree spacing
%   - >2 values         list of longitudes of the grid
% lats              latitude(s) which can be either
%   - []                entire earth (Default)
%   - single value      single point
%   - two values        [latmin latmax] with 1 degree spacing
%   - >2 values         list of latitudes of the grid
%
% OUTPUT:
% val               values of the variable at the lon/lat grid
%
% EXAMPLES:
% readcrust1([], 'rho', 9)
% readcrust1([], 'bnds', 3, [160 260], [-40 10])
% readcrust1('list')
% readcrust1([], 'lslayers')
% readcryst1([], 'lsbnds')
%
% Last modified by sirawich@spipatprathanporn@ucsd.edu, 06/24/2026

% list of constants
NUMLAYERS = 9;
NUMLONS   = 360;
NUMLATS   = 180;
LONLIM    = [-180 180];
LATLIM    = [-90 90];
LONS      = -179.5:179.5;
LATS      =  -89.5:89.5;

% list of default input values
defval('ddir', fullfile(getenv('IFILES'), 'EARTHMODELS', 'PHYSICAL', 'crust1.0'))
defval('lons', [])
defval('lats', [])

varnamelists = {{'rho', 'density'},                             ... %1
                {'vp', 'p-wave speed'},                         ... %2
                {'vs', 's-wave speed'},                         ... %3
                {'bnds', 'z', 'boundary', 'elev'},              ... %4
                {'thick', 'dz', 'thickness'},                   ... %5
                {'iswater'},                                    ... %6
                {'isice'},                                      ... %7
                {'lsbnds', 'list boundaries'},                  ... %8
                {'lslayers', 'list layers'}};                       %9

varnamedescriptions = {'density of the layer',                  ... %1
                       'P-wave speed of the layer',             ... %2
                       'S-wave speed of the layer',             ... %3
                       'Boundary surface elevation',            ... %4
                       'Layer thickness',                       ... %5
                       'Whether it is water',                   ... %6
                       'Whether it is ice',                     ... %7
                       'List boundary names and numbers',       ... %8
                       'List layer names and numbers'};             %9

listboundaries = {{'Earth surface', 'top of water'}, ...
                  {'ETOPO1-ice', 'Bathymetry', 'bottom of water'}, ...
                  {'ETOPO1-bedrock', 'bottom of ice'}, ...
                  {'bottom of sediments 1'}, ...
                  {'bottom of sediments 2'}, ...
                  {'bedrock', 'bottom of sediments 3'}, ...
                  {'bottom of crystalline crust 1'}, ...
                  {'bottom of crystalline crust 2'}, ...
                  {'Moho', 'bottom of crystalline crust 3'}};

listlayers = {{'water'}, ...
              {'ice'}, ...
              {'sediments 1 (upper)'}, ...
              {'sediments 2 (middle)'}, ...
              {'sediments 3 (lower)'}, ...
              {'crystalline crust 1 (upper)'}, ...
              {'crystalline crust 2 (middle)'}, ...
              {'crystalilne crust 3 (lower)'}, ...
              {'mantle (layer thickness not available)'}};

listlayers_extra = {{'all sediments (layer 3-5)'}, ...
                    {'all crystalline crust (layer 6-8)'}, ...
                    {'all crust (layer 2-8)'}};

if strcmpi(ddir, 'list') || strcmpi(varname, 'list')
    % print the listed variables
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                     varname(s) : descriptions               \n');
    fprintf('-----------------------------------------------------------------\n');
    for ii = 1:length(varnamelists)
        fprintf('%31s : %s\n', cell2commasepstr(varnamelists{ii}), ...
            varnamedescriptions{ii});
    end
    fprintf('-----------------------------------------------------------------\n');
    outputlist = {[]};
    varargout = outputlist(1:nargout);
    return
elseif any(strcmpi(varname, varnamelists{1}))
    % return empty if the layer doesn't have a density
    if layer <= 0 || layer > NUMLAYERS
        fprintf('Density is avilable in layer 1-%d only\n', NUMLAYERS);
        outputlist = {[]};
        varargout = outputlist(1:nargout);
        return
    end
    % identify which file to read
    fname = fullfile(ddir, 'crust1.rho');

    % read, slice for the layer, and reshape
    fid = fopen(fname, 'r');
    val = double(fscanf(fid, '%f', [NUMLAYERS Inf]));
    fclose(fid);

    val = val';
    val = flipud(reshape(val(:,layer), [NUMLONS NUMLATS])');

    % handle longitudes from 0 to 360
    if any(lons > 180)
        % split and join the grid, so 180E is in the middle
        val = val(:,[181:360 1:180]);
        LONS = mod(LONS(:,[181:360 1:180]), 360);
    end

elseif any(strcmpi(varname, varnamelists{2}))
    % return empty if the layer doesn't have a P-wave speed
    if layer <= 0 || layer > NUMLAYERS
        fprintf('P-wave speed is avilable in layer 1-%d only\n', NUMLAYERS);
        outputlist = {[]};
        varargout = outputlist(1:nargout);
        return
    end
    % identify which file to read
    fname = fullfile(ddir, 'crust1.vp');

    % read, slice for the layer, and reshape
    fid = fopen(fname, 'r');
    val = double(fscanf(fid, '%f', [NUMLAYERS Inf]));
    fclose(fid);

    val = val';
    val = flipud(reshape(val(:,layer), [NUMLONS NUMLATS])');

    % handle longitudes from 0 to 360
    if any(lons > 180)
        % split and join the grid, so 180E is in the middle
        val = val(:,[181:360 1:180]);
        LONS = mod(LONS(:,[181:360 1:180]), 360);
    end
elseif any(strcmpi(varname, varnamelists{3}))
    % return empty if the layer doesn't have a S-wave speed
    if layer <= 0 || layer > NUMLAYERS
        fprintf('S-wave speed is avilable in layer 1-%d only\n', NUMLAYERS);
        outputlist = {[]};
        varargout = outputlist(1:nargout);
        return
    end
    % identify which file to read
    fname = fullfile(ddir, 'crust1.vs');

    % read, slice for the layer, and reshape
    fid = fopen(fname, 'r');
    val = double(fscanf(fid, '%f', [NUMLAYERS Inf]));
    fclose(fid);

    val = val';
    val = flipud(reshape(val(:,layer), [NUMLONS NUMLATS])');

    % handle longitudes from 0 to 360
    if any(lons > 180)
        % split and join the grid, so 180E is in the middle
        val = val(:,[181:360 1:180]);
        LONS = mod(LONS(:,[181:360 1:180]), 360);
    end

elseif any(strcmpi(varname, varnamelists{4}))
    % return empty if the layer doesn't have a boundary surface elevation
    if layer <= 0 || layer > NUMLAYERS
        fprintf('Boundary surface elevation is avilable in layer 1-%d only\n', NUMLAYERS);
        outputlist = {[]};
        varargout = outputlist(1:nargout);
        return
    end
    % identify which file to read
    fname = fullfile(ddir, 'crust1.bnds');

    % read, slice for the layer, and reshape
    fid = fopen(fname, 'r');
    val = double(fscanf(fid, '%f', [NUMLAYERS Inf]));
    fclose(fid);

    val = val';
    val = flipud(reshape(val(:,layer), [NUMLONS NUMLATS])');

    % handle longitudes from 0 to 360
    if any(lons > 180)
        % split and join the grid, so 180E is in the middle
        val = val(:,[181:360 1:180]);
        LONS = mod(LONS(:,[181:360 1:180]), 360);
    end

elseif any(strcmpi(varname, varnamelists{5}))
    % return empty if the layer doesn't have a boundary surface elevation
    if (layer <= 0 || layer > NUMLAYERS + length(listlayers_extra)) || ...
            layer == 9
        fprintf('Layer thickness is not available for this layer.\n');
        fprintf('Try readcrust1([], "lslayers") to see list of layers.\n')
        outputlist = {[]};
        varargout = outputlist(1:nargout);
        return
    end
    
    % identify which file to read
    fname = fullfile(ddir, 'crust1.bnds');
    fid = fopen(fname, 'r');
    val = double(fscanf(fid, '%f', [NUMLAYERS Inf]));
    fclose(fid);

    val = val';

    if layer >= 1 && layer <= 8
        val_top = flipud(reshape(val(:,layer), [NUMLONS NUMLATS])');
        val_bot = flipud(reshape(val(:,layer+1), [NUMLONS NUMLATS])');
    elseif layer == 10
        val_top = flipud(reshape(val(:,3), [NUMLONS NUMLATS])');
        val_bot = flipud(reshape(val(:,6), [NUMLONS NUMLATS])');
    elseif layer == 11
        val_top = flipud(reshape(val(:,6), [NUMLONS NUMLATS])');
        val_bot = flipud(reshape(val(:,9), [NUMLONS NUMLATS])');
    elseif layer == 12
        val_top = flipud(reshape(val(:,2), [NUMLONS NUMLATS])');
        val_bot = flipud(reshape(val(:,9), [NUMLONS NUMLATS])');
    end
    val = val_top - val_bot;

    % handle longitudes from 0 to 360
    if any(lons > 180)
        % split and join the grid, so 180E is in the middle
        val = val(:,[181:360 1:180]);
        LONS = mod(LONS(:,[181:360 1:180]), 360);
    end
elseif any(strcmpi(varname, varnamelists{6}))
    % water: water layer is not zero thick
    waterthickness = readcrust1(ddir, 'thick', 1, lons, lats);
    val = (waterthickness > 0);

    outputlist = {val};
    varargout = outputlist(1:nargout);
    return
elseif any(strcmpi(varname, varnamelists{7}))
    % ice: ice layer is not zero thick
    icethickness = readcrust1(ddir, 'thick', 2, lons, lats);
    val = (icethickness > 0);

    outputlist = {val};
    varargout = outputlist(1:nargout);
    return
elseif any(strcmpi(varname, varnamelists{8}))
    % print the listed boundary numbers
    fprintf('-----------------------------------------------------------------\n');
    fprintf(' No. : boundary descriptions               \n');
    fprintf('-----------------------------------------------------------------\n');
    for ii = 1:length(listboundaries)
        fprintf('%4d : %s\n', ii, strjoin(listboundaries{ii}, ', '));
    end
    fprintf('-----------------------------------------------------------------\n');
    outputlist = {[]};
    varargout = outputlist(1:nargout);
    return
elseif any(strcmpi(varname, varnamelists{9}))
    % print the listed layer numbers
    fprintf('-----------------------------------------------------------------\n');
    fprintf(' No. : layer descriptions               \n');
    fprintf('-----------------------------------------------------------------\n');
    for ii = 1:length(listlayers)
        fprintf('%4d : %s\n', ii, strjoin(listlayers{ii}, ', '));
    end
    fprintf('-----------------------------------------------------------------\n');
    fprintf(' Only layer thickness is available in these layers.\n');
    fprintf('-----------------------------------------------------------------\n');
    for ii = 1:length(listlayers_extra)
        fprintf('%4d : %s\n', ii+length(listlayers), ...
            strjoin(listlayers_extra{ii}, ', '));
    end
    fprintf('-----------------------------------------------------------------\n');
    outputlist = {[]};
    varargout = outputlist(1:nargout);
    return
else
    outputlist = {[]};
    varargout = outputlist(1:nargout);
    return
end

% process the input longitudes and latitudes
if length(lons) == 2
    lons = lons(1):lons(2);
end
if length(lats) == 2
    lats = lats(1):lats(2);
end

% linearly interpolate for the value at the requested grid
if ~isempty(lats) && ~isempty(lons)
    [xx, yy] = meshgrid(lons, lats);
    val = interp2(LONS, LATS, val, xx, yy);
end
    
outputlist = {val};
varargout = outputlist(1:nargout);
return
end