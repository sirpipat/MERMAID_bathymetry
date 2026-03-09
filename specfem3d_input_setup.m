function specfem3d_input_setup(ddir, depth, topo, inject, solver, with2d)
% SPECFEM3D_INPUT_SETUP(ddir, depth, topo, inject, solver, with2d)
%
% Generates Par_file, cmtsolution file, station file, Mesh_Par_file, and 
% interface file for a fluid-solid simulation.
%
% INPUT:
% ddir          directory for the input files
% depth         depth of the hydrophone (deeper --> more negative)
%               [Default: -1500]
% topo          ocean bottom topography: struct containing following variables
% - type          'flat' [default] or
%                 'randomfield' (two-dimensional Matern covariance field)
% - depth         mean depth [Default: -4000]
% +++ ONLY FOR topo.type == 'randomfield' +++
% - s2            the first Matern parameter, aka sigma^2, in field units^2
%                 [Default: 100]
% - nu            the second Matern parameter [Default: 0.5]
% - rho           the third Matern parameter, in the units of the grid
%                 [Default: 1000]
% - dydx          sampling interval in the y and x directions [m m]
%                 [Default: [500 500]]
% - NyNx          number of samples in the y and x directions
%                 [Default: [57 57]]
% - blurs         0 No wavenumber blurring
%                 1 No wavenumber blurring, effectively
%                 N Fejer convolutional  BLUROS  on an N-times refined grid
%                -1 Fejer multiplicative BLUROSY using exact procedure
%               Inf Simulate using SGP invariant embedding [default]
% - taper         0 there is no taper near of far
%                 1 it's a unit taper, implicitly [default]
%                 OR an appropriately sized taper with proper values 
%                   (1 is yes and 0 is no and everything in between)
% inject        injection parameters: struct containing following variables
% - type          'Gaussian' [default] or 'stf'
% - theta         incident angle [Default: 0]
% - freq          scaling frequency of Gaussian waveform [Default: 1]
% - fs            sampling frequency of injected wave [Default: 10]
% - stf           source-time function [t d(t)]
%                 [Default: [(-100:0.01:100)', 1-Hz Gaussian]]
% - origin_time   origin time [Default: 0]
% solver        struct containing following variables
% - gpu_mode      whether to run on GPUs [Default: false]
% - nsteps        number of steps [Default: 10000]
% with2d        whether to create a SPECFEM2D run for comparison
%               The run will be the same as 3D except and the 2D box will
%               be a plane x=0 of the 3D box. [Default: false]
%
% SEE ALSO:
% SPECFEM3D_INPUT_SETUP_FLAT, SIMULOSL
%
% Last modified by sirawich-at-princeton.edu: 03/08/2026

% default parameters
defstruct('topo', 'type', 'flat')
    
defstruct('inject', 'type', 'Gaussian')
defstruct('inject', 'theta', 0)
defstruct('inject', 'freq', 1)
defstruct('inject', 'fs', 10)
defstruct('inject', 'origin_time', 0)

defstruct('solver', 'gpu_mode', false)
defstruct('solver', 'nsteps', 10000)

defval('with2d', false)

%% determine topography
if strcmp(topo.type, 'flat')
    defstruct('topo', 'depth', -4000)
    NXI =  2;
    NETA = 2;
    bath = ones(2, 2) * topo.depth;
elseif strcmp(topo.type, 'randomfield')
    defstruct('topo', 'depth', -4000)
    defstruct('topo', 's2', 100)
    defstruct('topo', 'nu', 0.5)
    defstruct('topo', 'rho', 1000);
    defstruct('topo', 'dydx', [500 500])
    defstruct('topo', 'NyNx', [57 57])
    defstruct('topo', 'blurs', Inf)
    defstruct('topo', 'taper', 1)
    defstruct('topo', 'shanning', [0.15 0.15])
    
    % generate a random field
    p.dydx = topo.dydx;
    p.NyNx = topo.NyNx - [16 16];
    p.blurs = topo.blurs;
    p.taper = topo.taper;
    Hx = simulosl([topo.s2 topo.nu topo.rho], p);
    Hxy = v2s(Hx, p);
    
    % apply hann filter
    wY = shanning(topo.NyNx(1)-16, topo.shanning(1));
    wX = shanning(topo.NyNx(2)-16, topo.shanning(2));
    Hxy = (wY * wX') .* Hxy;
    
    NXI = topo.NyNx(2);
    NETA = topo.NyNx(1);
    
    bath = ones(NETA, NXI) * topo.depth;
    bath(9:NETA-8, 9:NXI-8) = bath(9:NETA-8, 9:NXI-8) + Hxy;
    
elseif strcmp(topo.type, 'bathymetry')
    defstruct('topo', 'dydx', [500 500])
    defstruct('topo', 'NyNx', [57 57])
    defstruct('topo', 'shanning', [0.15 0.15])
    defstruct('topo', 'lon', -168.38)
    defstruct('topo', 'lat', -12.00)
    defstruct('topo', 'baz', 266.2356)
    
    NXI = topo.NyNx(2);
    NETA = topo.NyNx(1);
    
    [~, ~, zz] = bathymetryprofile2d(topo.dydx .* ...
        [topo.NyNx(1)-1 topo.NyNx(2)-1], [topo.NyNx(1) topo.NyNx(2)], ...
        [topo.lon topo.lat], topo.baz);
    
    % rotate counter clockwise 90 degrees
    % zz from bathymetryprofile2d, the wave enters from the left while
    % zz for this setup, the wave enters from the bottom
    zz = rot90(zz);
    
    % mean depth
    topo.depth = mean(zz, 'all');
    
    % apply hann filter
    wY = [zeros(4,1); shanning(topo.NyNx(1)-8, topo.shanning(1)); zeros(4,1)];
    wX = [zeros(4,1); shanning(topo.NyNx(2)-8, topo.shanning(2)); zeros(4,1)];
    bath = (wY * wX') .* (zz - topo.depth) + topo.depth;
    
% for direction test purpose for now
elseif strcmp(topo.type, 'slope')
    defstruct('topo', 'depth', -5000)
    defstruct('topo', 'slope', 0.1)
    defstruct('topo', 'dydx', [500 500])
    defstruct('topo', 'NyNx', [57 57])
    defstruct('topo', 'shanning', [0.15 0.15])
    
    NXI = topo.NyNx(2);
    NETA = topo.NyNx(1);
    
    % slope 
    Hxy = (1:NETA-16)' * topo.slope * topo.dydx(1) + (1:NXI-16) * 0;
    
    % apply hann filter
    wY = shanning(NETA-16, topo.shanning(1));
    wX = shanning(NXI-16, topo.shanning(2));
    Hxy = (wY * wX') .* Hxy;
    
    bath = ones(NETA, NXI) * topo.depth;
    bath(9:NETA-8, 9:NXI-8) = bath(9:NETA-8, 9:NXI-8) + Hxy;
end

% elevation of ocean bottom at the middle of the box
if strcmp(topo.type, 'flat')
    bottom = topo.depth;
else
    if mod(NETA, 2)
        if mod(NXI, 2)
            bottom = bath((NETA+1)/2, (NXI+1)/2);
            
        else
            bottom = mean(bath((NETA+1)/2, NXI/2:NXI/2+1));
        end
    else
        if mod(NXI, 2)
            bottom = mean(bath(NETA/2:NETA/2+1, (NXI+1)/2));
        else
            bottom = mean(bath(NETA/2:NETA/2+1, NXI/2:NXI/2+1), 'all');
        end
    end
end
%% process input parameters and write into files
% create the DATA/ folder for the files generated by this function
system(sprintf('mkdir -p %s', ddir));
system(sprintf('mkdir -p %s/DATA/', ddir));
system(sprintf('mkdir -p %s/DATA/meshfem3D_files/', ddir));

% xspecfem3d parameters
params = makeparams3d;
fparams = fullfile(ddir, 'DATA', 'Par_file');
params.NSTEP = solver.nsteps;
params.GPU_MODE = solver.gpu_mode;
writeparfile3d(params, fparams);

% fkmodel parameters
fkmodel = makefkmodel;
ffkmodel = fullfile(ddir, 'DATA', 'FKMODEL');
if strcmp(inject.type, 'Gaussian')
    inject.fs = max(inject.fs, 2*inject.freq);
elseif strcmp(inject.type, 'Ricker')
    inject.fs = max(inject.fs, 2*inject.freq);
    t_stf = (-100:0.01:100)';
    x_stf = (1 - 2 * inject.freq^2 * t_stf.^2) .* ...
        exp(-(inject.freq * t_stf).^2);
    inject.stf = [t_stf x_stf];
    stf_fname = 'stf_file.txt';
    writetimeseries(inject.stf(:,1), inject.stf(:,2), ...
        fullfile(ddir, stf_fname));
elseif strcmp(inject.type, 'stf')
    defstruct('inject', 'stf', [(-100:0.01:100)' ...
        exp(-((-100:0.01:100)' / 1).^2)])
    if ischar(inject.stf)
        stf_fname = removepath(inject.stf);
        system(sprintf('cp %s %s', inject.stf, fullfile(ddir, stf_fname)));
    elseif all(~isnan(inject.stf), 'all')
        stf_fname = 'stf_file.txt';
        writetimeseries(inject.stf(:,1), inject.stf(:,2), ...
            fullfile(ddir, stf_fname));
    end
end
fkmodel.nlayers = 2;
layer1 = struct('rho', 1020, 'vp', 1500, 'vs', 0, 'ztop', 0);
layer2 = struct('rho', 2500, 'vp', 3400, 'vs', 1963, 'ztop', topo.depth);
fkmodel.layers = {layer1; layer2};
fkmodel.baz = 180; % fixed for now
fkmodel.theta = inject.theta;
fkmodel.fmax = inject.freq;
fkmodel.fs = inject.fs;
fkmodel.twindow = params.NSTEP * params.DT;
fkmodel.origin_time = inject.origin_time;
fkmodel.origin_wavefront = [0 -14000 -10000];
if strcmp(inject.type, 'Gaussian')
    fkmodel.stf_type = 1;
    fkmodel.stf_file = "n/a";
else
    fkmodel.stf_type = 4;
    fkmodel.stf_file = stf_fname;
end
writefkmodel(fkmodel, ffkmodel);

% source parameters (unused)
cmt = makecmtsolution3d;
fcmt = fullfile(ddir, 'DATA', 'CMTSOLUTION');
writecmtsolution3d(cmt, fcmt);

% station
stations = struct('name', {{'OBS01', 'P0009'}}, ...
    'network', {{'AA', 'MH'}}, 'lat', [0 0], ...
    'lon', [0 0], 'elev', [bottom bottom], 'z', [bottom depth]);
fstations = fullfile(ddir, 'DATA', 'STATIONS');
writestations3d(stations, fstations);

% xmeshfem3d parameters
meshparams = makemeshparams3d;
fmesh = fullfile(ddir, 'DATA', 'meshfem3D_files', 'Mesh_Par_file');
% water
material1 = struct('rho', 1020, 'vp', 1500, 'vs', 0, 'Q_Kappa', 9999, ...
    'Q_mu', 9999, 'anisotropy_flag', 0, 'domain_id', 1);
% crust
material2 = struct('rho', 2500, 'vp', 3400, 'vs', 1963, ...
    'Q_Kappa', 9999, 'Q_mu', 9999, 'anisotropy_flag', 0, 'domain_id', 2);
materials = {material1; material2};
meshparams.NMATERIALS = 2;
meshparams.MATERIALS = materials;
% regions
nz_all = 20;
nz_top = round(nz_all*abs(topo.depth)/(meshparams.DEPTH_BLOCK_KM * 1000));
nz_bottom = nz_all - nz_top;
region1 = struct(...
    'NEX_XI_BEGIN'      , 1             , ...
    'NEX_XI_END'        , 56            , ...
    'NEX_ETA_BEGIN'     , 1             , ...
    'NEX_ETA_END'       , 56            , ...
    'NZ_BEGIN'          , nz_bottom+1   , ...
    'NZ_END'            , 20            , ...
    'material_id'       , 1               ...
);

region2 = struct(...
    'NEX_XI_BEGIN'      , 1             , ...
    'NEX_XI_END'        , 56            , ...
    'NEX_ETA_BEGIN'     , 1             , ...
    'NEX_ETA_END'       , 56            , ...
    'NZ_BEGIN'          , 1             , ...
    'NZ_END'            , nz_bottom     , ...
    'material_id'       , 2               ...
);
meshparams.NREGIONS = 2;
meshparams.REGIONS = {region1; region2};
writemeshparfile3d(meshparams, fmesh);

% interface file parameters
finterf = fullfile(ddir, 'DATA', 'meshfem3D_files', ...
    meshparams.INTERFACES_FILE);
% interfaces
% ocean bottom
itf1 = struct(...
    'SUPPRESS_UTM_PROJECTION'   , true  , ...
    'NXI'                       , NXI     , ...
    'NETA'                      , NETA     , ...
    'LON_MIN'   , meshparams.LONGITUDE_MIN, ...
    'LAT_MIN'   , meshparams.LATITUDE_MIN, ...
    'SPACING_XI', (meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN) / (NXI - 1), ...
    'SPACING_ETA', (meshparams.LATITUDE_MAX - meshparams.LATITUDE_MIN) / (NETA - 1), ...
    'FILE', fullfile(ddir, 'DATA', 'meshfem3D_files', 'interf_ocean.dat'), ...
    'Z', bath ...
);
% sea surface
itf2 = struct(...
    'SUPPRESS_UTM_PROJECTION'   , true  , ...
    'NXI'                       , 2   , ...
    'NETA'                      , 2  , ...
    'LON_MIN'   , meshparams.LONGITUDE_MIN, ...
    'LAT_MIN'   , meshparams.LATITUDE_MIN, ...
    'SPACING_XI', meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN, ...
    'SPACING_ETA', meshparams.LATITUDE_MAX - meshparams.LATITUDE_MIN, ...
    'FILE', fullfile(ddir, 'DATA', 'meshfem3D_files', 'interf_top.dat'), ...
    'Z', zeros(2,2) ...
);
itfs = {itf1; itf2};
layers = [nz_bottom; nz_top];
writeinterfacefiles3d(itfs, layers, finterf);

if with2d
    % output directory for specfem2d
    words = split(ddir, filesep);
    ddir2d = fullfile(getenv('REMOTE2D'), words{end-1}, words{end});
    ddir2d = [ddir2d filesep];
    
    system(sprintf('mkdir -p %s', ddir2d));
    
    % cross section of bathymetry profile at x=0
    tparams.X = linspace(meshparams.LATITUDE_MIN, meshparams.LATITUDE_MAX, NETA)';
    tparams.Z = bath(:,ceil(NXI/2));
    
    % convert (X,Z) coordinates to those in SPECFEM2D
    % (0,0) in SPECFEM3D == (10000,9600) in SPECFEM2D 
    tparams.X = tparams.X + 10000;
    tparams.Z = tparams.Z +  9600;
    
    if strcmp(inject.type, 'Ricker')
        time_function_type = 1;
    else
        time_function_type = 3;
    end
    
    specfem2d_input_setup_response('from3d', ...
                                   'custom', ...
                                   tparams, ...
                                   -depth, ...
                                   'homogeneous', ...
                                   'homogeneous', ...
                                   inject.freq / pi, ... % Gaussian in 2D: exp(-(pi * t * f).^2)
                                   inject.theta, ...
                                   time_function_type, ...
                                   [], ...
                                   ddir2d, ...
                                   false, ...
                                   'devel', ...
                                   solver.gpu_mode);
end
end