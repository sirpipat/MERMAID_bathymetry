function writeinterfacefiles3d(itfs, layers, fname)
% WRITEINTERFACEFILES3D(itfs, layers, fname)
%
% Writes interfaces and layers to an interface file and the elevation files
% used by this interface file. These elevation files will be in the same
% directory as the interface file as specified by FNAME variable. The
% directories of ITFS{ii}.FILE are not used here.
%
% Read https://specfem3d.readthedocs.io/en/latest/03_mesh_generation on how
% the interfaces are stored in files
%
% DISCLAIMER: This is not the official way to read/write an interface 
% file. I just go through comments and parameters in an instant of 
% interface file and read/write accordingly.
%
% INPUT:
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
% fname         full-path name of the interface file you want to create
%               Default: [] -- everything is written to standard output
%
% SEE ALSO:
% LOADINTERFACEFILES3D
%
% Last modified by sirawich-at-princeton.edu, 07/13/2026

defval('fname', [])

%% open the file
if isempty(fname)
    % standard output aka console output
    fid = 1;
    ddir = [];
else
    fid = fopen(fname, 'w');
    ddir = strcat(fileparts(fname), filesep);
end

%% write the interfaces information
fprintf(fid, '# number of interfaces\n');
fprintf(fid, '%d\n', length(itfs));
fprintf(fid, '#\n');
fprintf(fid, ['# We describe each interface below, structured as a ' ...
    '2D-grid, with several parameters :\n']);
fprintf(fid, ['# number of points along XI and ETA, minimal XI ETA ' ...
    'coordinates\n']);
fprintf(fid, '# and spacing between points which must be constant.\n');
fprintf(fid, ['# Then the records contain the Z coordinates of the ' ...
    'NXI x NETA points.\n']);
fprintf(fid, '#\n');
for ii = 1:length(itfs)
    fprintf(fid, '# interface number %d\n', ii);
    fprintf(fid, ['# SUPPRESS_UTM_PROJECTION  NXI  NETA LONG_MIN   ' ...
        'LAT_MIN    SPACING_XI SPACING_ETA\n']);
    if itfs{ii}.SUPPRESS_UTM_PROJECTION
        bool_word = '.true.';
    else
        bool_word = '.false.';
    end
    fprintf(fid, '%-22s %d %d %g %g %g %g\n', bool_word, ...
        itfs{ii}.NXI, itfs{ii}.NETA, itfs{ii}.LON_MIN, ...
        itfs{ii}.LAT_MIN, itfs{ii}.SPACING_XI, itfs{ii}.SPACING_ETA);
    fprintf(fid, '%s\n', removepath(itfs{ii}.FILE));
    fprintf(fid, '\n');
    
    % write the topography file
    if isempty(ddir)
        fid_ii = 1;
    else
        fid_ii = fopen(strcat(ddir, removepath(itfs{ii}.FILE)), 'w');
    end
    if itfs{ii}.SUPPRESS_UTM_PROJECTION
        fprintf(fid_ii, '%g\n', reshape(itfs{ii}.Z, ...
            itfs{ii}.NXI * itfs{ii}.NETA, 1));
    else
        fprintf(fid_ii, '%g\n', reshape(itfs{ii}.Z', ...
            itfs{ii}.NXI * itfs{ii}.NETA, 1));
    end
    if fid_ii >= 3
        fclose(fid_ii);
    end
end

%% write the layers information
fprintf(fid, '\n');
fprintf(fid, ['# for each layer, we give the number of spectral ' ...
    'elements in the vertical direction\n']);
for ii = 1:length(layers)
    fprintf(fid, '# layer number %d', ii);
    if ii == 1
        fprintf(fid, ' (bottom layer)\n');
    elseif ii == length(layers)
        fprintf(fid, ' (top layer)\n');
    else
        fprintf(fid, '\n');
    end
    fprintf(fid, '%d\n', layers(ii));
end

%% close the file
if fid >= 3
    fclose(fid);
end
end