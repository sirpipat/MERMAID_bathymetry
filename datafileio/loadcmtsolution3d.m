function cmt = loadcmtsolution3d(fname)
% cmt = LOADCMTSOLUTION3D(fname)
% cmt = LOADCMTSOLUTION3D(ddir)
%
% Loads a CMTSOLUTION file of SPECFEM3D_Cartesian.
%
% INPUT
% fname         name of the CMTSOLUTION file
% ddir          directory the FK-SPECFEM3D run where the file is
%               ddir/DATA/CMTSOLUTION
%
% OUTPUT
% cmt           CMT solutions, array of struct(s) with following fields
%               - HDR
%                   - TYPE
%                   - YEAR
%                   - MONTH
%                   - DAY
%                   - HOUR
%                   - MINUTE
%                   - SECOND
%                   - LAT
%                   - LON
%                   - DEPTH
%                   - mb
%                   - Ms
%                   - NAME
%               - NAME              event name
%               - TSHIFT            time shift
%               - HALFDUR           half duration
%               - LAT               latitude
%               - LON               longitude
%               - DEPTH             depth
%               - Mrr
%               - Mtt
%               - Mpp
%               - Mrt
%               - Mrp
%               - Mtp
%
% The CMT solution format can be found at
% https://specfem3d.readthedocs.io/en/latest/05_running_the_solver/
%
% SEE ALSO:
% LOADSOURCE, MAKECMTSOLUTION3D, WRITECMTSOLUTION3D
%
% Last modified by sirawich-at-princeton.edu, 10/13/2026

cmt = {};
num_sources = 0;

fid = fopen(fname, 'r');
if fid == -1
    fprintf('It is not a file name. Assume it is a directory:\n');
    fname = fullfile(fname, 'DATA', 'CMTSOLUTION');
    fprintf('Read from %s\n', fname);
    fid = fopen(fname, 'r');
end
line = fgetl(fid);

while ~isempty(line)
    % skip comments / headers
    while isempty(line) || strcmp(line(1), '#')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    % read the header line
    header.TYPE = strip(line(1:4));
    words = split(strip(line(5:end)));
    header.YEAR = sscanf(words{1}, '%d', 1);
    header.MONTH = sscanf(words{2}, '%d', 1);
    header.DAY = sscanf(words{3}, '%d', 1);
    header.HOUR = sscanf(words{4}, '%d', 1);
    header.MINUTE = sscanf(words{5}, '%d', 1);
    header.SECOND = sscanf(words{6}, '%f', 1);
    header.LAT = sscanf(words{7}, '%f', 1);
    header.LON = sscanf(words{8}, '%f', 1);
    header.DEPTH = sscanf(words{9}, '%f', 1);
    header.mb = sscanf(words{10}, '%f', 1);
    header.Ms = sscanf(words{11}, '%f', 1);
    header.NAME = cindeks(join(words(12:end)), 1);
    source.HDR = header;
    
    % read the rest of the CMT solution
    line = fgetl(fid);
    source.NAME = readstring(line);
    
    line = fgetl(fid);
    source.TSHIFT = readfloat(line);
    
    line = fgetl(fid);
    source.HALFDUR = readfloat(line);
    
    line = fgetl(fid);
    source.LAT = readfloat(line);
    
    line = fgetl(fid);
    source.LON = readfloat(line);
    
    line = fgetl(fid);
    source.DEPTH = readfloat(line);
    
    line = fgetl(fid);
    source.Mrr = readfloat(line);
    
    line = fgetl(fid);
    source.Mtt = readfloat(line);
    
    line = fgetl(fid);
    source.Mpp = readfloat(line);
    
    line = fgetl(fid);
    source.Mrt = readfloat(line);
    
    line = fgetl(fid);
    source.Mrp = readfloat(line);
    
    line = fgetl(fid);
    source.Mtp = readfloat(line);
    
    % add the source to cmt array
    num_sources = num_sources + 1;
    cmt{num_sources} = source;
    
    % read next source
    line = fgetl(fid);
end

fclose(fid);
end

function value = readstring(line)
% find equal sign
where_start = strfind(line, ':');

% find # where the comment starts
where_end = strfind(line, '#');
if isempty(where_end)
    where_end = length(line) + 1;
end

% read the value
value = strip(sscanf(line((where_start+1):(where_end-1)), '%c'));
end


function value =  readfloat(line)
% find equal sign
where = strfind(line, ':');
% change the exponent notation syntax from 'd' to 'e'
line = replace(line, 'd', 'e');
% read the value
value = sscanf(line((where+1):end), '%f', 1);
end