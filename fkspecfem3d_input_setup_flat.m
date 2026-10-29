function fkspecfem3d_input_setup_flat(ddir, depth, fkmodel, model, stf, solver)
% FKSPECFEM3D_INPUT_SETUP_FLAT(ddir, depth, fkmodel, model, stf, solver)
%
% Sets up a generic FK-SPECFEM3D run with homogenous flat layered media.
% Then, you can customize your model using SETBOXSIZE, ADDOCEANBOTTOM,
% REMESH3D etc.
%
% INPUT:
% ddir          directory for the FK-SPECFEM3D run
% depth         depth of the hydrophone float
% fkmodel       the fkmodel from FKMODEL file
% model         model name. See options below
% - 'OC', 'PS2024'      Ocean + Crust                         [default for 2 layers]
% - 'OCM'               Ocean + Crust + Mantle                [default for 3 layers]
% - 'OSCM'              Ocean + Sediment + Crust + Mantle     [default for 4 layers]
% - 'OS3C3M', 'Crust1'  Ocean + 3 Sediment + 3 Crust + Mantle [default for 8 layers]
% - 'OSmCnM'            Ocean + m Sediment + n Crust + Mantle
% - 'RIO(...)'          aiR + Ice + Ocean + ...
% - '(...)CnLA'         ... + n Crust + Lithosphere + Aesthenosphere
% stf           source-time function [t d(t)]
%               [Default: [(-100:0.01:100)', 5-Hz Gaussian]]
% solver        struct containing following variables
% - fs              sampling rate
% - nstep           number of time steps
% - NPROC_XI        number of processors in X-direction
% - NPROC_ETA       number of processors in Y-direction
% - gpu_mode        whether to use GPU mode
%
% SEE ALSO:
% SPECFEM3D_INPUT_SETUP_FLAT, MAKEFKMODEL,
% SETBOXSIZE, ADDOCEANBOTTOM, ADDLAYER, REMESH3D, SETLAYER,
% SETSEDIMENTLAYER
%
% Last modified by spipatprathanporn@ucsd.edu, 10/29/2026

defval('stf', [(-100:0.01:100)' exp(-((-100:0.01:100)' / 0.2).^2)])

numlayers = fkmodel.nlayers;

if ~exist('model', 'var') || isempty(model)
    switch numlayers
        case 2
            model = 'OC';
        case 3
            model = 'OCM';
        case 4
            model = 'OSCM';
        case 8
            model = 'OS3C3M';
        otherwise
            model = sprintf('OSC%dM', numlayers-3);
    end
end
if strcmpi(model, 'PS2024')
    model = 'OC';
elseif strcmpi(model, 'Crust1')
    model = 'OS3C3M';
end

% make a list of appropriate filenames
filenames = cell(1,numlayers);
filenames{1} = 'interf_top.dat';
numfiles = 1;
letters = split(model, '');
letters = letters(2:end-1);
for ii = 2:length(letters)
    switch upper(letters{ii})
        case 'R'
            numfiles = numfiles + 1;
            filenames{numfiles} = 'interf_air.dat';
        case 'I'
            numfiles = numfiles + 1;
            filenames{numfiles} = 'interf_ice.dat';
        case 'O'
            numfiles = numfiles + 1;
            filenames{numfiles} = 'interf_seasurface.dat';
        case 'S'
            numfiles = numfiles + 1;
            filenames{numfiles} = 'interf_ocean.dat';

            % check if there are more than 1 sediment layers
            m = '';
            iletter = ii;
            while iletter < length(letters) && ...
                ~isnan(str2double(letters{iletter+1}))
                iletter = iletter+1;
                m = [m letters{iletter}];
            end
            m = str2double(m);
            if ~isnan(m)
                for jj = 2:m
                    numfiles = numfiles + 1;
                    filenames{numfiles} = ...
                        sprintf('interf_sediment%d.dat', jj);
                end
            end
        case 'C'
            % check if there are more than 1 crustal layers
            n = '';
            iletter = ii;
            while iletter < length(letters) && ...
                    ~isnan(str2double(letters{iletter+1}))
                iletter = iletter+1;
                n = [n letters{iletter}];
            end
            n = str2double(n);
            if ~isnan(n)
                for jj = 1:n
                    numfiles = numfiles + 1;
                    if jj == 1 && strcmpi(letters{ii-1}, 'O')
                        filenames{numfiles} = 'interf_ocean.dat';
                    else
                        filenames{numfiles} = ...
                            sprintf('interf_crust%d.dat', jj);
                    end
                end
            else
                numfiles = numfiles + 1;
                if strcmpi(letters{ii-1}, 'O')
                    filenames{numfiles} = 'interf_ocean.dat';
                else
                    filenames{numfiles} = 'interf_crust.dat';
                end
            end
        case {'M', 'L'}
            numfiles = numfiles + 1;
            filenames{numfiles} = 'itnerf_moho.dat';
        case 'A'
            numfiles = numfiles + 1;
            filenames{numfiles} = 'interf_lab.dat';
        case {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}
        otherwise
            numfiles = numfiles + 1;
            filenames{numfiles} = sprintf('interf_layer%02d.dat', numfiles);
    end
end
if numfiles ~= numlayers
    error('model name does not match the number of layers in fkmodel');
end

% make solver parameters
params = makeparams3d;
params.DT = 1 / solver.fs;
params.NSTEP = solver.nstep;
params.NPROC = solver.NPROC_XI * solver.NPROC_ETA;


% list of materials
MATERIALS = cell(numlayers,1);
for ii = 1:numlayers
    MATERIALS{ii}.rho             = fkmodel.layers{ii}.rho;
    MATERIALS{ii}.vp              = fkmodel.layers{ii}.vp;
    MATERIALS{ii}.vs              = fkmodel.layers{ii}.vs;
    MATERIALS{ii}.Q_Kappa         = 9999;
    MATERIALS{ii}.Q_mu            = 9999;
    MATERIALS{ii}.anisotropy_flag = 0;
    if fkmodel.layers{ii}.vs > 0
        MATERIALS{ii}.domain_id = 2;
    else
        MATERIALS{ii}.domain_id = 1;
    end
end

% make mesh parameters
meshparams = struct(...
               'LATITUDE_MIN', -14000               , ...
               'LATITUDE_MAX', 14000                , ...
              'LONGITUDE_MIN', -14000               , ...
              'LONGITUDE_MAX', 14000                , ...
             'DEPTH_BLOCK_KM', -fkmodel.origin_wavefront(3) / 1000, ...
        'UTM_PROJECTION_ZONE', 11                   , ...
    'SUPPRESS_UTM_PROJECTION', 1                    , ...
            'INTERFACES_FILE', 'interfaces.dat'     , ...
                'CAVITY_FILE', 'no_cavity.dat'      , ...
                     'NEX_XI', 56                   , ...
                    'NEX_ETA', 56                   , ...
                   'NPROC_XI', solver.NPROC_XI      , ...
                  'NPROC_ETA', solver.NPROC_ETA     , ...
           'USE_REGULAR_MESH', 1                    , ...
                 'NDOUBLINGS', 0                    , ...
               'NZ_DOUBLINGS', [11; 0]              , ...
        'CREATE_ABAQUS_FILES', 0                    , ...
            'CREATE_DX_FILES', 0                    , ...
           'CREATE_VTK_FILES', 1                    , ...
         'SAVE_MESH_AS_CUBIT', 0                    , ...
                 'LOCAL_PATH', './DATABASES_MPI'    , ...
         'THICKNESS_OF_X_PML', 12.3000              , ...
         'THICKNESS_OF_Y_PML', 12.3000              , ...
         'THICKNESS_OF_Z_PML', 12.3000              , ...
                 'NMATERIALS', numlayers            , ...
                  'MATERIALS', {MATERIALS}          , ...
                   'NREGIONS', numlayers            , ...
                    'REGIONS', {cell(numlayers, 1)}   ...
    );
% modify the regions
for ii = 1:numlayers
    meshparams.REGIONS{ii}.NEX_XI_BEGIN  = 1;
    meshparams.REGIONS{ii}.NEX_XI_END    = meshparams.NEX_XI;
    meshparams.REGIONS{ii}.NEX_ETA_BEGIN = 1;
    meshparams.REGIONS{ii}.NEX_ETA_END   = meshparams.NEX_ETA;
    meshparams.REGIONS{ii}.NZ_BEGIN      = numlayers+1-ii;
    meshparams.REGIONS{ii}.NZ_END        = numlayers+1-ii;
    meshparams.REGIONS{ii}.material_id   = ii;
end

% make interface parameters
layers = ones(numlayers, 1);
itfs = cell(numlayers, 1);
for ii = 1:numlayers
    itfs{ii}.SUPPRESS_UTM_PROJECTION = 1;
    itfs{ii}.NXI = 2;
    itfs{ii}.NETA = 2;
    itfs{ii}.LON_MIN = meshparams.LONGITUDE_MIN;
    itfs{ii}.LAT_MIN = meshparams.LATITUDE_MIN;
    itfs{ii}.SPACING_XI = meshparams.LONGITUDE_MAX - ...
        meshparams.LONGITUDE_MIN;
    itfs{ii}.SPACING_ETA = meshparams.LATITUDE_MAX - ...
        meshparams.LATITUDE_MIN;
    itfs{ii}.FILE = fullfile(ddir, 'DATA', 'meshfem3D_files', ...
        filenames{numlayers+1-ii});
    itfs{ii}.Z = repmat(fkmodel.layers{numlayers+1-ii}.ztop, [2 2]);
end

%% save the parameters
% create the DATA/ folder for the files generated by this function
system(sprintf('mkdir -p %s', ddir));
system(sprintf('mkdir -p %s/DATA/', ddir));
system(sprintf('mkdir -p %s/DATA/meshfem3D_files/', ddir));


parfile = fullfile(ddir, 'DATA', 'Par_file');
fkfile = fullfile(ddir, params.FKMODEL_FILE);
meshparfile = fullfile(ddir, 'DATA', 'meshfem3D_files', 'Mesh_Par_file');
interfacefile = fullfile(ddir, 'DATA', 'meshfem3D_files', ...
    meshparams.INTERFACES_FILE);
stf_fname = fullfile(ddir, fkmodel.stf_file);

writefkmodel(fkmodel, fkfile);
writeparfile3d(params, parfile);
writemeshparfile3d(meshparams, meshparfile);
writeinterfacefiles3d(itfs, layers, interfacefile);

% source parameters (unused)
cmt = makecmtsolution3d;
fcmt = fullfile(ddir, 'DATA', 'CMTSOLUTION');
writecmtsolution3d(cmt, fcmt);

% station
bottom = fkmodel.layers{strcmp(filenames, 'interf_ocean.dat')}.ztop;
stations = struct('name', {{'OBS01', 'P0009'}}, ...
    'network', {{'AA', 'MH'}}, 'lat', [0 0], ...
    'lon', [0 0], 'elev', [bottom -depth], ...
    'z', [bottom -depth]);
fstations = fullfile(ddir, 'DATA', 'STATIONS');
writestations3d(stations, fstations);

if fkmodel.stf_type == 4
    writetimeseries(stf(:,1), stf(:,2), stf_fname);
end
end