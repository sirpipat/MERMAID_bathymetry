function fkmodel = makefkmodel(model, varargin)
% fkmodel = MAKEFKMODEL
% fkmodel = MAKEFKMODEL(model)
% fkmodel = MAKEFKMODEL(model, varargin)
%
% Makes a sensible fkmodel for FK-SPECFEM3D_Cartesian.
%
% INPUT:
% model         model name
%   - 'PS2024'    Ocean bottom from Pipatprathanporn & Simons 2024 (default)
%   - 'OCM'       Ocean + Crust + Mantle
%   - 'OSCM'      Ocean + Sediment + Crust + Mantle
%   - 'Crust1'    Example crust1.0 profile (entry #40720)
%   - [lat lon]   Nearest crust1.0 profile (not implemented yet)
% varargin      keyword argument pairs to modify fkmodel
%
% OUTPUT:
% fkmodel       model for embeded FK modelling
%
% SEE ALSO:
% LOADFKMODEL, WRITEFKMODEL
%
% Last modified by sirawich-at-princeton.edu, 10/29/2026

defval('model', 'PS2024')

if isstring(model) || ischar(model)
    switch lower(model)
        case 'ps2024'
            fkmodel.nlayers = 2;
            rho  = [1020  2500];
            vp   = [1500  3400];
            vs   = [   0  1963];
            ztop = [   0 -4400];
        case 'oc'
            fkmodel.nlayers = 2;
            rho  = [1020  2850];
            vp   = [1500  6500];
            vs   = [   0  3700];
            ztop = [   0 -4400];
        case 'ocm'
            fkmodel.nlayers = 3;
            rho  = [1020  2850   3330];
            vp   = [1500  6500   8080];
            vs   = [   0  3700   4490];
            ztop = [   0 -4400 -11370];
        case 'oscm'
            fkmodel.nlayers = 4;
            rho  = [1020  1820  2850   3330];
            vp   = [1500  1750  6500   8080];
            vs   = [   0   540  3700   4490];
            ztop = [   0 -4400 -4650 -11370];
        case 'crust1'
            fkmodel.nlayers = 6;
            rho  = [1020  1820  2550  2850  3050   3330];
            vp   = [1500  1750  5000  6500  7100   8080];
            vs   = [   0   340  2700  3700  4050   4490];
            ztop = [   0 -4400 -4650 -5130 -6660 -11370];
        otherwise
            fkmodel.nlayers = 3;
            rho = [1020 2677.5 3300];
            vp = [1500 5250 8031.15];
            vs = [0 2835 4576.15];
            ztop = [0 -5000 -7500];
    end
else
    % Nearest crust1.0 profile
    ii_lat = min(max(round(89.5 - model(1)), 0), 179);
    ii_lon = min(max(round(model(2) + 179.5), 0), 359);
    
    lineno = 360 * ii_lat + ii_lon + 1;
    
    crustdir = fullfile(getenv('IFILES'), 'EARTHMODELS', 'PHYSICAL', ...
        'crust1.0');
    
    [~, line] = system(sprintf('head -n %d %s%s%s | tail -n 1', lineno, ...
        crustdir, filesep, 'crust1.bnds'));
    ztop = str2num(line);
    
    [~, line] = system(sprintf('head -n %d %s%s%s | tail -n 1', lineno, ...
        crustdir, filesep, 'crust1.rho'));
    rho = str2num(line);
    
    [~, line] = system(sprintf('head -n %d %s%s%s | tail -n 1', lineno, ...
        crustdir, filesep, 'crust1.vp'));
    vp = str2num(line);
    
    [~, line] = system(sprintf('head -n %d %s%s%s | tail -n 1', lineno, ...
        crustdir, filesep, 'crust1.vs'));
    vs = str2num(line);
    
    % remove unused layers
    wh   = [(ztop(2:end) - ztop(1:end-1) ~= 0) true];
    ztop = ztop(wh) * 1000;
    rho  = rho(wh) * 1000;
    vp   = vp(wh) * 1000;
    vs   = vs(wh) * 1000;
    
    fkmodel.nlayers = length(ztop);
end

for ii = 1:fkmodel.nlayers
    fkmodel.layers{ii,1} = struct('rho', rho(ii), 'vp', vp(ii), ...
        'vs', vs(ii), 'ztop', ztop(ii));
end

% incident wave
fkmodel.wave = 'p';

% angle of incoming wave
fkmodel.baz = 270;
fkmodel.theta = 0;

% frequency max
fkmodel.fmax = 1;

% sampling frequency
fkmodel.fs = 10;

% time window
fkmodel.twindow = 50;

% amplitude
fkmodel.amplitude = 1;

% optional
fkmodel.origin_wavefront = nan(1,3);
fkmodel.origin_time = nan;

fkmodel.stf_type = nan;
fkmodel.stf_file = "n/a";

% modify fkmodel from keyword argument pairs
for ii = 1:2:length(varargin)
    if isfield(fkmodel, lower(varargin{ii}))
        fkmodel.(lower(varargin{ii})) = varargin{ii+1};
    else
        warning('%s is not a variable of fkmodel.', varargin{ii})
    end
end
end