function meshparams = loadmeshparfile3d(fname)
% meshparams = LOADMESHPARFILE3D(fname)
% meshparams = LOADMESHPARFILE3D(ddir)
%
% Reads mesh parameters from a Mesh_Par_file of SPECFEM3D_Cartesian.
%
% INPUT:
% fname             name of a Mesh_Par_file
% ddir              directory the FK-SPECFEM3D run where the Mesh par file
%                   is: ddir/DATA/meshfem3D_files/Mesh_Par_file
%
% OUTPUT:
% meshparams        mesh parameters
%
% SEE ALSO:
% MAKEMESHPARAMS3D, WRITEMESHPARFILE3D, LOADPARFILE3D
%
% Last modified by sirawich-at-princeton.edu, 10/13/2026

% open the file
fid = fopen(fname, 'r');
if fid == -1
    fprintf('It is not a file name. Assume it is a directory:\n');
    fprintf('Read from %s\n', ...
        fullfile(fname, 'DATA', 'meshfem3D_files', 'Mesh_Par_file'))
    fid = fopen(fullfile(fname, 'DATA', 'meshfem3D_files', ...
        'Mesh_Par_file'), 'r');
end

names = {};
values = {};
numvar = 0;

% Read the parameters
line = fgetl(fid);
while ischar(line)
    % skip comments / headers
    while isempty(line) || strcmp(line(1), '#') || ~contains(line, '=')
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
    end
    if ~ischar(line)
        break;
    end
    
    numvar = numvar + 1;
    [name, value] = readgeneric(line);
    names{numvar} = name;
    values{numvar} = value;
    
    % handle special parameter input formats (anything but NAME = value)
    if strcmp(name, 'NDOUBLINGS')
        % skip the comment (instruction)
        line = fgetl(fid);
        while isempty(line) || strcmp(line(1), '#')
            line = fgetl(fid);
        end
        if ~ischar(line)
            break;
        end
        
        % there must be at least 2 NZ_DOUBLINGS_n entries regardless of the
        % value of NDOUBLINGS
        NDOUBLINGS = max(2, value);
        name = 'NZ_DOUBLINGS';
        NZ_DOUBLINGS = nan(NDOUBLINGS, 1);
        % read the doubling layers
        for ii = 1:NDOUBLINGS
            NZ_DOUBLINGS(ii) = readint(line);
            if ii < NDOUBLINGS
                line = fgetl(fid);
            end
        end
        
        numvar = numvar + 1;
        names{numvar} = name;
        values{numvar} = NZ_DOUBLINGS;
    elseif strcmp(name, 'NMATERIALS')
        % skip the comment (instruction)
        line = fgetl(fid);
        while isempty(line) || strcmp(line(1), '#')
            line = fgetl(fid);
        end
        if ~ischar(line)
            break;
        end
        
        NMATERIALS = value;
        name = 'MATERIALS';
        MATERIALS = cell(NMATERIALS, 1);
        % read the material
        for ii = 1:NMATERIALS
            nums = sscanf(replace(line, 'd', 'e'), '%f');
            MATERIALS{ii}.rho               = nums(2);
            MATERIALS{ii}.vp                = nums(3);
            MATERIALS{ii}.vs                = nums(4);
            MATERIALS{ii}.Q_Kappa           = nums(5);
            MATERIALS{ii}.Q_mu              = nums(6);
            MATERIALS{ii}.anisotropy_flag   = nums(7);
            MATERIALS{ii}.domain_id         = nums(8);
            if ii < NMATERIALS
                line = fgetl(fid);
            end
        end
        
        numvar = numvar + 1;
        names{numvar} = name;
        values{numvar} = MATERIALS;
    elseif strcmp(name, 'NREGIONS')
        % skip the comment (instruction)
        line = fgetl(fid);
        while isempty(line) || strcmp(line(1), '#')
            line = fgetl(fid);
        end
        if ~ischar(line)
            break;
        end
        
        NREGIONS = value;
        name = 'REGIONS';
        REGIONS = cell(NREGIONS, 1);
        % read the domain regions
        for ii = 1:NREGIONS
            nums = sscanf(line, '%d');
            REGIONS{ii}.NEX_XI_BEGIN    = nums(1);
            REGIONS{ii}.NEX_XI_END      = nums(2);
            REGIONS{ii}.NEX_ETA_BEGIN   = nums(3);
            REGIONS{ii}.NEX_ETA_END     = nums(4);
            REGIONS{ii}.NZ_BEGIN        = nums(5);
            REGIONS{ii}.NZ_END          = nums(6);
            REGIONS{ii}.material_id     = nums(7);
            if ii < NREGIONS
                line = fgetl(fid);
            end
        end
        
        numvar = numvar + 1;
        names{numvar} = name;
        values{numvar} = REGIONS;
    end
    
    % read the next line
    line = fgetl(fid);
end

% close the files
fclose(fid);

meshparams = cell2struct(values, names, 2);
end

function [name, value] = readgeneric(line)
% list of variables with known data type
bool_var   = {'SUPPRESS_UTM_PROJECTION', ...
              'USE_REGULAR_MESH', ...
              'CREATE_ABAQUS_FILES', ...
              'CREATE_DX_FILES', ...
              'CREATE_VTK_FILES', ...
              'SAVE_MESH_AS_CUBIT'
              };
int_var    = {'UTM_PROJECTION_ZONE', ...
              'NEX_XI', ...
              'NEX_ETA', ...
              'NPROC_XI', ...
              'NPROC_ETA', ...
              'NDOUBLINGS', ...
              'NMATERIALS', ...
              'NREGIONS'
              };
float_var  = {'LATITUDE_MIN', ...
              'LATITUDE_MAX', ...
              'LONGITUDE_MIN', ...
              'LONGITUDE_MAX', ...
              'DEPTH_BLOCK_KM', ...
              'THICKNESS_OF_X_PML', ...
              'THICKNESS_OF_Y_PML', ...
              'THICKNESS_OF_Z_PML'
              };
string_var = {'INTERFACES_FILE', ...
              'CAVITY_FILE', ...
              'LOCAL_PATH'
              };

name = sscanf(line, '%s', 1);

if any(strcmp(bool_var, name))
    value = readbool(line);
elseif any(strcmp(int_var, name))
    value = readint(line);
elseif any(strcmp(float_var, name))
    value = readfloat(line);
elseif any(strcmp(string_var, name))
    value = readstring(line);
else
    fprintf('unable to determine type of "%s" variable. Read as string\n', ...
        name);
    keyboard;
    value = readstring(line);
end
end

function value = readstring(line)
% find equal sign
where_start = strfind(line, '=');

% find # where the comment starts
where_end = strfind(line, '#');
if isempty(where_end)
    where_end = length(line) + 1;
end

% read the value
value = strip(sscanf(line((where_start+1):(where_end-1)), '%c'));
end

function value = readbool(line)
value = readstring(line);
if strcmp(value, '.true.')
    value = true;
elseif strcmp(value, '.false.')
    value = false;
else
    % do not know what to do
    error(strcat('ValueError: cannot read a boolean\n', ...
                 sprintf('line >> %s\n', line)));
end
end

function value = readint(line)
% find equal sign
where = strfind(line, '=');
% read the value
value = sscanf(line((where+1):end), '%d', 1);
end

function value =  readfloat(line)
% find equal sign
where = strfind(line, '=');
% change the exponent notation syntax from 'd' to 'e'
line = replace(line, 'd', 'e');
% read the value
value = sscanf(line((where+1):end), '%f', 1);
end