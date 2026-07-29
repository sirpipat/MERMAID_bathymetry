function writestations3d(stations, fname)
% WRITESTATIONS3D(stations, fname)
%
% Writes a STATION file of SPECFEM3D_Cartesian.
% 
% INPUT:
% stations          struct containing following fields
%       name            station name
%       network         network name
%       lat or y        latitude  or y-coordinate
%       lon or x        longitude or x-coordinate
%       elev            elevation
%       z               burial or z-coorndiate or depth 
% fname             filename of a STATION file
%
% SEE ALSO:
% LOADSTATIONS3D, READSTATIONS3D, MAKESTATIONS3D
%
% Last modified by sirawich-at-princeton.edu, 07/29/2026

if isempty(fname)
    % standard output aka console output
    fid = 1;
else
    fid = fopen(fname, 'w');
end

if ~isfield(stations, 'lat')
    stations.lat = stations.y;
end
if ~isfield(stations, 'lon')
    stations.lon = stations.x;
end

for ii = 1:length(stations.name)
    fprintf(fid, '%5s %2s %11.4f %11.4f %11.4f %11.4f\n', ...
        stations.name{ii}, stations.network{ii}, stations.lat(ii), ...
        stations.lon(ii), stations.elev(ii), stations.z(ii));
end

% close the file
if fid >= 3
    fclose(fid);
end
end