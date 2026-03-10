function fkmodel = loadfkmodel(fname)
% fkmodel = LOADFKMODEL(fname)
%
% Loads a FKMODEL file of SPECFEM3D_Cartesian.
%
% List of parameters is taken from:
% src/specfem3D/couple_with_injection.f90 in SPECFEM3D_Cartesian code
%
% INPUT:
% fname         name of the FKMODEL file
%
% OUTPUT:
% fkmodel       struct containing elastic properties of the media in
%               FK-domain and the properties of the incoming wave
%
% SEE ALSO:
% MAKEFKMODEL, WRITEFKMODEL
%
% Last modified by sirawich-at-princeton.edu, 03/09/2026

% set default values
fkmodel = makefkmodel;

fid = fopen(fname, 'r');
while true
    line = fgetl(fid);
    % skip comments
    while strcmp(line(1), '#')
        line = fgetl(fid);
    end
    % exit at the end of file
    if ~isempty(line) && ~ischar(line)
        break
    end
    
    % check variable names in FK model file
    words = split(line);
    switch words{1}
        case 'NLAYER'
            fkmodel.nlayers = str2double(words{2});
        case 'LAYER'
            fkmodel.layers = cell(fkmodel.nlayers, 1);
            fkmodel.layers{1,1}.rho  = str2double(words{3});
            fkmodel.layers{1,1}.vp   = str2double(words{4});
            fkmodel.layers{1,1}.vs   = str2double(words{5});
            fkmodel.layers{1,1}.ztop = str2double(words{6});
            for ii = 2:fkmodel.nlayers
                line = fgetl(fid);
                words = split(line);
                if strcmp(words{1}, 'LAYER')
                    fkmodel.layers{ii}.rho  = str2double(words{3});
                    fkmodel.layers{ii}.vp   = str2double(words{4});
                    fkmodel.layers{ii}.vs   = str2double(words{5});
                    fkmodel.layers{ii}.ztop = str2double(words{6});
                else
                    warning('One or more layers are missing');
                    break
                end
            end
        case 'INCIDENT_WAVE'
            fkmodel.wave = words{2};
        case 'BACK_AZIMITH'
            fkmodel.baz = str2double(words{2});
        case 'AZIMUTH'
            fkmodel.baz = mod(str2double(words{2}) + 180, 360);
        case 'TAKE_OFF'
            fkmodel.theta = str2double(words{2});
        case 'ORIGIN_WAVEFRONT'
            fkmodel.origin_wavefront = [str2double(words{2}) ...
                                        str2double(words{3}) ...
                                        str2double(words{4})];
        case 'ORIGIN_TIME'
            fkmodel.origin_time = str2double(words{2});
        case 'FREQUENCY_MAX'
            fkmodel.fmax = str2double(words{2});
        case 'FREQUENCY_SAMPLING'
            fkmodel.fs = str2double(words{2});
        case 'TIME_WINDOW'
            fkmodel.twindow = str2double(words{2});
        case 'AMPLITUDE'
            fkmodel.amplitude = str2double(words{2});
        case 'TIME_FUNCTION_TYPE'
            fkmodel.stf_type = str2double(words{2});
        case 'NAME_OF_SOURCE_FILE'
            fkmodel.stf_file = words{2};
    end
end

fclose(fid);
end