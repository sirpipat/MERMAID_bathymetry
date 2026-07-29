function stations = loadstations3d(fname)
% stations = LOADSTATIONS3D(fname)
% 
% Reads a STATION file for SPECFEM3D simulation and stores as a struct,
% fully compatitble with WRITESTAIONS3D
%
% INPUT
% fname         full filename of a STATION file
%
% OUTPUT
% stations      a struct containing following fields
%     network       station network names
%     name          station names
%     x             x-coordinates or longitudes
%     y             y-coordinates or latitudes
%     elev          elevations
%     z             z-coordinates
%
% SEE ALSO:
% READSTATIONS3D, WRITESTATIONS3D, MAKESTATIONS3D, READ_STATIONS
%
% Last modified by Sirawich Pipatprathanporn, 07/29/2026

% read the station file as a table
opts = detectImportOptions(fname, 'FileType', 'text');
stations = readtable(fname, opts);

% if the table is empty, assume the table has only one entry
if isempty(stations)
    fid = fopen(fname, 'r');
    line = fgetl(fid);
    words = split(line);
    stations = struct('network', {words(2)}, 'name', {words(1)}, ...
        'x', str2double(words{4}), ...
        'y', str2double(words{3}), ...
        'elev', str2double(words{5}), ...
        'z', str2double(words{6}));
else
    stations = struct('network', {stations.Var2}, 'name', {stations.Var1}, ...
        'x', stations.Var4, ...
        'y', stations.Var3, ...
        'elev', stations.Var5, ...
        'z', stations.Var6);
end
end