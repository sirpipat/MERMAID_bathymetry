function p = crust1profile(lon, lat, method, ddir)
% p = CRUST1PROFILE(lon, lat, method)
% p = CRUST1PROFILE(lon, lat, method, ddir)
%
% Returns a crustal profile at (lon,lat) according to CRUST1.0 model.
%
% INPUT:
% ddif          name of the CRUST1.0 directory
%               [Default: $IFILES/EARTHMODELS/PHYSICAL/crust1.0]
% lon           longitude of the profile
% lat           latitude of the profile
% method        interpolation method [Default: 'nearest']
%
% OUTPUT:
% p             crustal profile: a cell array of structs (layer) with a 
%               following member variables:
%       - rho       density of medium of the layer     (kg/m^3)
%       - vp        P-wave speed                       (m/s)
%       - vs        S-wave speed                       (m/s)
%       - ztop      Elevation of the top of the layer  (m, above sea level)
%
% SEE ALSO:
% MAKEFKMODEL
%
% Last modified by spipatprathanporn@ucsd.edu, 06/29/2026

% list of constants
NUMLAYERS = 9;
NUMLONS   = 360;
NUMLATS   = 180;
LONLIM    = [-180 180];
LATLIM    = [-90 90];
LONS      = -179.5:179.5;
LATS      =  -89.5:89.5;

% list of default input values
defval('method', 'linear')
defval('ddir', fullfile(getenv('IFILES'), 'EARTHMODELS', 'PHYSICAL', 'crust1.0'))

% ensure that -180 <= lon <= 180
lon = mod(lon + 180, 360) - 180;

% interpolated to the nearest grid point
if strcmpi(method, 'nearest')
    lon = interp1(LONS, LONS, lon, 'nearest');
    lat = interp1(LATS, LATS, lat, 'nearest');
end

% read the crust1.0 files
fid = fopen(fullfile(ddir, 'crust1.rho'));
rho = fscanf(fid, '%f', [NUMLAYERS Inf]);
rho = rho';
fclose(fid);

fid = fopen(fullfile(ddir, 'crust1.vp'));
vp = fscanf(fid, '%f', [NUMLAYERS Inf]);
vp = vp';
fclose(fid);

fid = fopen(fullfile(ddir, 'crust1.vs'));
vs = fscanf(fid, '%f', [NUMLAYERS Inf]);
vs = vs';
fclose(fid);

fid = fopen(fullfile(ddir, 'crust1.bnds'));
bnds = fscanf(fid, '%f', [NUMLAYERS Inf]);
bnds = bnds';
fclose(fid);

% crustal profile
p = repmat(struct('rho', nan, 'vp', nan, 'vs', nan, 'ztop', nan),...
    NUMLAYERS, 1);

for ii = 1:NUMLAYERS
    rholayer = flipud(reshape(rho(:,ii), [NUMLONS NUMLATS])');
    p(ii).rho = interp2(LONS, LATS, rholayer, lon, lat, method) * 1000;

    vplayer = flipud(reshape(vp(:,ii), [NUMLONS NUMLATS])');
    p(ii).vp = interp2(LONS, LATS, vplayer, lon, lat, method) * 1000;

    vslayer = flipud(reshape(vs(:,ii), [NUMLONS NUMLATS])');
    p(ii).vs = interp2(LONS, LATS, vslayer, lon, lat, method) * 1000;

    bndslayer = flipud(reshape(bnds(:,ii), [NUMLONS NUMLATS])');
    p(ii).ztop = interp2(LONS, LATS, bndslayer, lon, lat, method) * 1000;
end

% remove zero-thick layers
ii = 1;
while ii <= length(p)
    if p(ii).rho == 0
        p(ii) = [];
    elseif ii < length(p) && p(ii).ztop - p(ii+1).ztop == 0
        p(ii) = [];
    else
        ii = ii+1;
    end
end

% convert to a cell array to match the data type of the FKMODEL
p = mat2cell(p, ones(size(p)));
end