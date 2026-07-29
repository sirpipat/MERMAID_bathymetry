function [n, name, network, x, y, z] = readstations3d(fname)
% [n, name, network, x, y, z] = READSTATIONS3D(fname)
% 
% Reads a STATION file for SPECFEM3D simulation and returns as vectors of
% column names. Consider using LOADSTATIONS3D if you want a vector of
% struct output.
%
% INPUT
% fname         full filename of a STATION file
%
% OUTPUT
% n             the number of stations
% name          station names
% network       station network names
% x             x-coordinates or longitude
% y             y-coordinates or latitude
% z             z-coordinates
%
% SEE ALSO:
% LOADSTATIONS3D, WRITESTATIONS3D, MAKESTATIONS3D, READ_STATIONS
%
% Last modified by Sirawich Pipatprathanporn, 07/29/2026

% read the station file as a table
opts = detectImportOptions(fname, 'FileType', 'text');
T = readtable(fname, opts);

% if the table is empty, assume the table has only one entry
if isempty(T)
    fid = fopen(fname, 'r');
    line = fgetl(fid);
    words = split(line);
    T = struct('Var1', {words(1)}, 'Var2', {words(2)}, ...
        'Var3', str2double(words{3}), ...
        'Var4', str2double(words{4}), ...
        'Var6', str2double(words{6}));
end

% Format
% name     network   y/lat    x/lon   elev    z  
name = T.Var1;
network = T.Var2;
x = T.Var4;
y = T.Var3;
z = T.Var6;
n = size(T, 1);
end