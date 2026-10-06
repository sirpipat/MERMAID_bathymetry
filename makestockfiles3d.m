function stockdir = makestockfiles3d(runmode)
% MAKESTOCKFILES3D(runmode)
%
% Create stock input files for SPECFEM3D runs.
%
% INPUT:
% runmode           see options
%       'default'   generate stock input files for the first ime
%       'notopo'    do not generate random field topography files
%       [s2 nu rho it]     generate random field topography files given a 
%                          variance, smoothness, correlation length, and
%                          the number of iterations
%
% Last modified by sirawich-at-princeton.edu, 10/06/2026

defval('runmode', 'default')

%% create stockfiles folder
stockdir = fullfile(getenv('REMOTE3D'), 'stockfiles');
system(sprintf('mkdir -p %s', stockdir));

fkdir = fullfile(stockdir, 'fkmodelfiles');
system(sprintf('mkdir -p %s', fkdir));

stfdir = fullfile(stockdir, 'stffiles');
system(sprintf('mkdir -p %s', stfdir));

pardir = fullfile(stockdir, 'parfiles');
system(sprintf('mkdir -p %s', pardir));

stadir = fullfile(stockdir, 'stationfiles');
system(sprintf('mkdir -p %s', stadir));

cmtdir = fullfile(stockdir, 'cmtsolutionfiles');
system(sprintf('mkdir -p %s', cmtdir));

meshdir = fullfile(stockdir, 'meshparfiles');
system(sprintf('mkdir -p %s', meshdir));

interfdir = fullfile(stockdir, 'interfacefiles');
system(sprintf('mkdir -p %s', interfdir));

topodir = fullfile(stockdir, 'topofiles');
system(sprintf('mkdir -p %s', topodir));

%% source-time function(s)
t_stf = (-100:0.01:100)';
freq = [0.1 0.5 1 2 4 10];

for ii = 1:length(freq)
    % Gaussian
    x_stf = exp(-(t_stf * freq(ii)).^2);
    fname = sprintf('stf_file_Gaussian_%s.txt', float2filestr(freq(ii)));
    writetimeseries(t_stf, x_stf, fullfile(stfdir, fname));
    
    % Ricker
    x_stf = (1 - 2 * freq(ii)^2 * t_stf.^2) .* exp(-(freq(ii) * t_stf).^2);
    fname = sprintf('stf_file_Ricker_%s.txt', float2filestr(freq(ii)));
    writetimeseries(t_stf, x_stf, fullfile(stfdir, fname));
end

%% FKMODEL(s)
theta = [0 10 40];

% Pipatprathanporn & Simons 2024 GJI
for ii = 1:length(freq)
    for jj = 1:length(theta)
        fkmodel = makefkmodel('PS2024', ...
            'twindow', 100, ...
            'stf_type', 4, ...
            'stf_file', 'stf_file.txt', ...
            'fmax', freq(ii), ...
            'theta', theta(jj), ...
            'origin_wavefront', [0 -10000 -10000], ...
            'origin_time', 0);
        
        fname = sprintf('FKMODEL_PS2024_THETA-%02d_FREQ-%s', theta(jj), ...
            float2filestr(freq(ii)));
        writefkmodel(fkmodel, fullfile(fkdir, fname));
        
        % Ocean + Crust
        fkmodel = makefkmodel('OC', ...
            'twindow', 100, ...
            'stf_type', 4, ...
            'stf_file', 'stf_file.txt', ...
            'fmax', freq(ii), ...
            'theta', theta(jj), ...
            'origin_wavefront', [0 -10000 -20000], ...
            'origin_time', 0);
        
        fname = sprintf('FKMODEL_OC_THETA-%02d_FREQ-%s', theta(jj), ...
            float2filestr(freq(ii)));
        writefkmodel(fkmodel, fullfile(fkdir, fname));
        
        % Ocean + Crust + Mantle
        fkmodel = makefkmodel('OCM', ...
            'twindow', 100, ...
            'stf_type', 4, ...
            'stf_file', 'stf_file.txt', ...
            'fmax', freq(ii), ...
            'theta', theta(jj), ...
            'origin_wavefront', [0 -10000 -20000], ...
            'origin_time', 0);

        fname = sprintf('FKMODEL_OCM_THETA-%02d_FREQ-%s', theta(jj), ...
            float2filestr(freq(ii)));
        writefkmodel(fkmodel, fullfile(fkdir, fname));

        % Ocean + Sediment + Crust + Mantle
        fkmodel = makefkmodel('OSCM', ...
            'twindow', 100, ...
            'stf_type', 4, ...
            'stf_file', 'stf_file.txt', ...
            'fmax', freq(ii), ...
            'theta', theta(jj), ...
            'origin_wavefront', [0 -10000 -20000], ...
            'origin_time', 0);

        fname = sprintf('FKMODEL_OSCM_THETA-%02d_FREQ-%s', theta(jj), ...
            float2filestr(freq(ii)));
        writefkmodel(fkmodel, fullfile(fkdir, fname));

        % Ocean + Sediment + 3 Crusts + Mantle
        fkmodel = makefkmodel('Crust1', ...
            'twindow', 100, ...
            'stf_type', 4, ...
            'stf_file', 'stf_file.txt', ...
            'fmax', freq(ii), ...
            'theta', theta(jj), ...
            'origin_wavefront', [0 -10000 -20000], ...
            'origin_time', 0);

        fname = sprintf('FKMODEL_Crust1_THETA-%02d_FREQ-%s', theta(jj), ...
            float2filestr(freq(ii)));
        writefkmodel(fkmodel, fullfile(fkdir, fname));
    end
end

%% Par_files
DTs = [0.01 0.005 0.0025 0.001];
for ii = 1:length(DTs)
    params = makeparams3d;
    params.DT = DTs(ii);
    params.NSTEP = round(70 / DTs(ii));
    params.GPU_MODE = true;
    fname = sprintf('Par_file_FS-%d', 1/DTs(ii));
    writeparfile3d(params, fullfile(pardir, fname));
end

%% STATION FILES
stations = struct('name', {{'OBS01', 'P0009'}}, ...
    'network', {{'AA', 'MH'}}, 'lat', [0 0], ...
    'lon', [0 0], 'elev', [-1500 -1500], 'z', [-4400 -1500]);
writestations3d(stations, fullfile(stadir, 'STATIONS'));

%% CMTSOLUTION FILES
cmt = makecmtsolution3d;
writecmtsolution3d(cmt, fullfile(cmtdir, 'CMTSOLUTION'));

%% MESH FILES
% size of the box in meters [x (lon), y (lat), z (depth)]
% The origin (0,0,0) is at the center of the top face
size_box = [28000 28000 10000; 20000 20000 20000];
size_box_name = {'original', 'thick'};

% element size
size_elem = [500 500 500; 250 250 250; 100 100 100];

% number of processors
NPROCS = [2 4; 4 4; 4 8];

% model name
model_name = {'PS2024', 'OC', 'OCM', 'OSCM', 'Crust1'};

for ii = 1:size(size_box, 1)
    for jj = 1:size(size_elem, 1)
        for kk = 1:size(NPROCS, 1)
            for ll = 1:length(model_name)
                % run only model 'PS2024' for original box size
                if ii == 1 && ll >= 2
                    continue
                end
                % skip model with sedimentary layer and elem size >= 500 m
                if jj == 1 && ll >= 4
                    continue
                end
                % materials: must consistent with layers in fkmodel
                meshparams3d = makemeshparams3d;
                meshparams3d.LATITUDE_MIN  = -size_box(ii, 2) / 2;
                meshparams3d.LATITUDE_MAX  =  size_box(ii, 2) / 2;
                meshparams3d.LONGITUDE_MIN = -size_box(ii, 1) / 2;
                meshparams3d.LONGITUDE_MAX =  size_box(ii, 1) / 2;
                meshparams3d.DEPTH_BLOCK_KM = size_box(ii, 3) / 1000;
                
                meshparams3d.NEX_XI = size_box(ii, 1) / size_elem(jj, 1);
                meshparams3d.NEX_ETA = size_box(ii, 2) / size_elem(jj, 2);
                meshparams3d.NPROC_XI = NPROCS(kk, 1);
                meshparams3d.NPROC_ETA = NPROCS(kk, 2);
                
                % materials + regions
                fkmodel = makefkmodel(model_name{ll});
                meshparams3d.NMATERIALS = fkmodel.nlayers;
                meshparams3d.NREGIONS = fkmodel.nlayers;
                for mm = fkmodel.nlayers:-1:1
                    layer = fkmodel.layers{mm};
                    meshparams3d.MATERIALS{mm} = struct('rho', layer.rho, ...
                        'vp', layer.vp, 'vs', layer.vs, 'Q_Kappa', 9999, ...
                        'Q_mu', 9999, 'anisotropy_flag', 0);
                    if layer.vs == 0
                        meshparams3d.MATERIALS{mm}.domain_id = 1;
                    else
                        meshparams3d.MATERIALS{mm}.domain_id = 2;
                    end
                    
                    meshparams3d.REGIONS{mm} = struct('NEX_XI_BEGIN', 1, ...
                        'NEX_XI_END', meshparams3d.NEX_XI, ...
                        'NEX_ETA_BEGIN', 1, ...
                        'NEX_ETA_END', meshparams3d.NEX_ETA);
                    if mm == fkmodel.nlayers
                        meshparams3d.REGIONS{mm}.NZ_BEGIN = 1;
                    else
                        meshparams3d.REGIONS{mm}.NZ_BEGIN = ...
                            meshparams3d.REGIONS{mm+1}.NZ_END + 1;
                    end
                    meshparams3d.REGIONS{mm}.NZ_END = round(...
                        (size_box(ii, 3) + layer.ztop) / size_elem(jj, 3));
                    meshparams3d.REGIONS{mm}.material_id = mm;
                end
                
                % name
                fname = sprintf(['Mesh_par_file_SIZE-%s_ELEM-%d_NPROCS' ...
                    '_X-%d_NPROCS_Y-%d_MODEL-%s'], size_box_name{ii}, ...
                    size_elem(jj, 1), NPROCS(kk, 1), NPROCS(kk, 2), ...
                    model_name{ll});
                
                % write mesh parameter file
                writemeshparfile3d(meshparams3d, fullfile(meshdir, fname));
            end
        end
    end
end


%% INTERFACES
% There is no makeinterface function here :(

% flat topography
subdir = fullfile(interfdir, 'FLAT');
system(sprintf('mkdir -p %s', subdir));
for kk = 1:length(model_name)
    fkmodel = makefkmodel(model_name{kk});
    for ii = 1:size(size_box, 1)
        for jj = 1:size(size_elem, 1)
            try
                fname_in = sprintf(['Mesh_par_file_SIZE-%s_ELEM-%d_NPROCS' ...
                    '_X-%d_NPROCS_Y-%d_MODEL-%s'], size_box_name{ii}, ...
                    size_elem(jj, 1), 2, 4, model_name{kk});
                meshparams3d = loadmeshparfile3d(fullfile(meshdir, fname_in));
            catch
                continue
            end
            layers = zeros(meshparams3d.NREGIONS, 1);
            for mm = fkmodel.nlayers:-1:1
                ii_itfs = fkmodel.nlayers - mm + 1;
                itfs{ii_itfs} = struct(...
                    'SUPPRESS_UTM_PROJECTION'   , true  , ...
                    'NXI'                       , 2   , ...
                    'NETA'                      , 2  , ...
                    'LON_MIN'   , meshparams3d.LONGITUDE_MIN, ...
                    'LAT_MIN'   , meshparams3d.LATITUDE_MIN, ...
                    'SPACING_XI', meshparams3d.LONGITUDE_MAX - meshparams3d.LONGITUDE_MIN, ...
                    'SPACING_ETA', meshparams3d.LATITUDE_MAX - meshparams3d.LATITUDE_MIN, ...
                    'FILE', fullfile(subdir, ...
                        sprintf(['interf_SIZE-%s_ELEM-%d_MODEL' ...
                        '-%s_%02d.dat'], size_box_name{ii}, ...
                        size_elem(jj, 1), model_name{kk}, ii_itfs)), ...
                    'Z', ones(2,2) * fkmodel.layers{mm}.ztop ...
                    );
                layers(ii_itfs) = meshparams3d.REGIONS{mm}.NZ_END - ...
                    meshparams3d.REGIONS{mm}.NZ_BEGIN + 1;
            end
            
            % write the interface file and its dependents
            fname = sprintf('interfaces_SIZE-%s_ELEM-%d_MODEL-%s.dat', ...
                size_box_name{ii}, size_elem(jj, 1), model_name{kk});
            writeinterfacefiles3d(itfs, layers, fullfile(subdir, fname));
            
            % reset the interfaces struct
            clear itfs
        end
    end
end



% NOTE: I separate the topography files (files with a bunch of elevations)
% from interface files to reduce the number of files and storage needed.
% For random field bathymetry or custom bathymetry, you need to combine an
% interface file and topography files together.
% create a directory for the inferface files and "topo" file
subdir = fullfile(interfdir, 'RANDOM');
system(sprintf('mkdir -p %s', subdir));

% ii,jj,kk loop over model, box size, and element size
for kk = 1:length(model_name)
    fkmodel = makefkmodel(model_name{kk});
    for ii = 1:size(size_box, 1)
        for jj = 1:size(size_elem, 1)
            try
                fname_in = sprintf(['Mesh_par_file_SIZE-%s_ELEM-%d_NPROCS' ...
                    '_X-%d_NPROCS_Y-%d_MODEL-%s'], size_box_name{ii}, ...
                    size_elem(jj, 1), 2, 4, model_name{kk});
                meshparams3d = loadmeshparfile3d(fullfile(meshdir, fname_in));
            catch
                continue
            end
            layers = zeros(meshparams3d.NREGIONS, 1);

            % assumes flat interfaces except the ocean bottom
            for mm = fkmodel.nlayers:-1:1
                ii_itfs = fkmodel.nlayers - mm + 1;
                if ii_itfs ~= fkmodel.nlayers - 1
                    itfs{ii_itfs} = struct(...
                        'SUPPRESS_UTM_PROJECTION'   , true  , ...
                        'NXI'                       , 2   , ...
                        'NETA'                      , 2  , ...
                        'LON_MIN'   , meshparams3d.LONGITUDE_MIN, ...
                        'LAT_MIN'   , meshparams3d.LATITUDE_MIN, ...
                        'SPACING_XI', meshparams3d.LONGITUDE_MAX - meshparams3d.LONGITUDE_MIN, ...
                        'SPACING_ETA', meshparams3d.LATITUDE_MAX - meshparams3d.LATITUDE_MIN, ...
                        'FILE', fullfile(subdir, ...
                            sprintf(['interf_SIZE-%s_ELEM-%d_MODEL' ...
                            '-%s_%02d.dat'], size_box_name{ii}, ...
                            size_elem(jj, 1), model_name{kk}, ii_itfs)), ...
                        'Z', ones(2,2) * fkmodel.layers{mm}.ztop ...
                        );
                else
                    NXI = meshparams3d.NEX_XI + 1;
                    NETA = meshparams3d.NEX_ETA + 1;
                    itfs{ii_itfs} = struct(...
                        'SUPPRESS_UTM_PROJECTION'   , true  , ...
                        'NXI'                       , NXI   , ...
                        'NETA'                      , NETA  , ...
                        'LON_MIN'   , meshparams3d.LONGITUDE_MIN, ...
                        'LAT_MIN'   , meshparams3d.LATITUDE_MIN, ...
                        'SPACING_XI'                , size_elem(jj, 1), ...
                        'SPACING_ETA'               , size_elem(jj, 2), ...
                        'FILE', fullfile(subdir, ...
                            sprintf(['interf_SIZE-%s_ELEM-%d_MODEL' ...
                            '-%s_%02d.dat'], size_box_name{ii}, ...
                            size_elem(jj, 1), model_name{kk}, ii_itfs)), ...
                        'Z', ones(NETA,NXI) * fkmodel.layers{mm}.ztop ...
                        );
                end
                layers(ii_itfs) = meshparams3d.REGIONS{mm}.NZ_END - ...
                    meshparams3d.REGIONS{mm}.NZ_BEGIN + 1;
                    
            end
            
            % write the interface file and its dependents
            fname = sprintf('interfaces_SIZE-%s_ELEM-%d_MODEL-%s.dat', ...
                size_box_name{ii}, size_elem(jj, 1), model_name{kk});
            writeinterfacefiles3d(itfs, layers, fullfile(subdir, fname));

            % reset the interfaces struct
            clear itfs
        end
    end
end

%% TOPOGRAPHY
% random topography
if strcmp(runmode, 'default')
    s2 = [100 1000 10000 100000];
    nu = [0.5 1.0 1.5];
    rho = [1000 2000 4000 8000];
    it = 1:10;
elseif strcmp(runmode, 'notopo')
    return
elseif isnumeric(runmode) && length(runmode) == 4
    s2 = runmode(1);
    nu = runmode(2);
    rho = runmode(3);
    it = 1:runmode(4);
end

% xx,yy,zz loops over random field's parameters
for xx = 1:length(s2)
    for yy = 1:length(nu)
        for zz = 1:length(rho)
            % nn loops over iteration
            for nn = it
                % Generate random field for the ocean bottom
                % Note that I make the grid of the field much larger than what 
                % I use, so that I can have larger rho
                MIN_SIZE_ELEM = 50;
                NEX_MAX = max(size_box(:,1)) / MIN_SIZE_ELEM;
                p.dydx = MIN_SIZE_ELEM * [1 1];
                p.NyNx = (1.125 * (rho(zz) / 1000) * NEX_MAX + 1) * [1 1];
                p.blurs = Inf;
                p.taper = 1;
                Hx = simulosl([s2(xx) nu(yy) rho(zz)], p);
                Hxy = v2s(Hx, p);

                % ii,jj,kk loop over model, box size, and element size
                for ii = 1:size(size_box, 1)
                    for jj = 1:size(size_elem, 1)
                        try
                            fname_in = sprintf(['Mesh_par_file_SIZE-%s_ELEM-%d_NPROCS' ...
                                '_X-%d_NPROCS_Y-%d_MODEL-%s'], size_box_name{ii}, ...
                                size_elem(jj, 1), 2, 4, 'PS2024');
                            meshparams3d = loadmeshparfile3d(fullfile(meshdir, fname_in));
                        catch
                            fname_in = sprintf(['Mesh_par_file_SIZE-%s_ELEM-%d_NPROCS' ...
                                '_X-%d_NPROCS_Y-%d_MODEL-%s'], size_box_name{ii}, ...
                                size_elem(jj, 1), 2, 4, 'OCM');
                            meshparams3d = loadmeshparfile3d(fullfile(meshdir, fname_in));
                        end
            
                        % apply the random field to the ocean bottom
                        % downsampling the random field
                        index_used = size_box(ii, 1) / MIN_SIZE_ELEM / 2;
                        downsampling_factor = size_elem(jj, 1) / MIN_SIZE_ELEM;
                        index_center = ceil(size(Hxy, 1) / 2);
                        index_used = ...
                            (-index_used:downsampling_factor:index_used)...
                            + index_center;
                        bath = Hxy(index_used, index_used);
                            

                        % apply Hann Filter
                        NXI = meshparams3d.NEX_XI + 1;
                        NETA = meshparams3d.NEX_ETA + 1;

                        % write the interface file and its dependents
                        fname = sprintf(['topo_RANDOM_S2-%d_' ...
                            'NU-%.2f_RHO-%d_SIZE-%d_ELEM-%d_IT-%03d.dat'], ...
                            s2(xx), nu(yy), rho(zz), ...
                            size_box(ii, 1), size_elem(jj, 1), nn);
                        fid = fopen(fullfile(topodir, fname), 'w');
                        fprintf(fid, '%g\n', reshape(bath', ...
                            NXI * NETA, 1));
                        fclose(fid);

                        % draw the bathymetry for reference
                        if ii > 1
                            figure(1)
                            clf
                            set(gcf, 'Units', 'inches', ...
                                'Position', [0 1 7 6])
                            imagesc([-10000 10000], [-10000 10000], ...
                                bath);
                            colormap(kelicol);
                            cb = colorbar( 'TickDirection', 'out', ...
                                'FontSize', 12, 'Ticks', ...
                                (-4:4) * round(sqrt(s2(xx)), 1, ...
                                'significant'));
                            cb.Label.String = 'elevation (m)';
                            grid on
                            titlestring = sprintf(['\\sigma^2 = %d' ...
                                ' | \\nu = %.2f | \\rho = %d | ' ...
                                '\\Deltax = %d'], s2(xx), nu(yy), ...
                                rho(zz), size_elem(jj, 1));
                            title(titlestring)
                            xlabel('easting (m)')
                            ylabel('northing (m)')
                            set(gca, 'TickDir', 'out', ...
                                'DataAspectRatio', [1 1 1], ...
                                'XTick', (-10:2:10) * 1e3, ...
                                'YTick', (-10:2:10) * 1e3, ...
                                'CLim', [-4 4] * round(sqrt(s2(xx)) ...
                                , 1, 'significant'), ...
                                'FontSize', 12)
                            set(gcf, 'Renderer', 'painters')

                            figname = sprintf(['%s_RANDOM_S2-%d_' ...
                                'NU-%.2f_RHO-%d_ELEM-%d_IT-%03d'], ...
                                mfilename, s2(xx), nu(yy), rho(zz), ...
                                size_elem(jj, 1), nn);
                            figname = replace(figname, '.', 'p');
                            figdisp(figname, [], [], 2, [], 'epstopdf');
                        end
                    end
                end
            end
        end
    end
end

%% COMMENTED OUT
% for xx = 1:length(s2)
%     for yy = 1:length(nu)
%         for zz = 1:length(rho)
%             % nn loops over iteration
%             for nn = 1:100
%                 % create a directory for the inferface files and "topo" file
%                 dirname = sprintf('RANDOM_S2-%d_NU-%.2f_RHO-%d_IT-%03d', ...
%                     s2(xx), nu(yy), rho(zz), nn);
%                 subdir = fullfile(interfdir, dirname);
%                 system(sprintf('mkdir -p %s', subdir));
% 
%                 % Generate random field for the ocean bottom
%                 % Note that I make the grid of the field much larger than what 
%                 % I use, so that I can have larger rho
%                 MIN_SIZE_ELEM = 50;
%                 NEX_MAX = max(size_box(:,1)) / MIN_SIZE_ELEM;
%                 p.dydx = MIN_SIZE_ELEM * [1 1];
%                 p.NyNx = (1.125 * (rho(zz) / 1000) * NEX_MAX + 1) * [1 1];
%                 p.blurs = Inf;
%                 p.taper = 1;
%                 Hx = simulosl([s2(xx) nu(yy) rho(zz)], p);
%                 Hxy = v2s(Hx, p);
% 
%                 % ii,jj,kk loop over model, box size, and element size
%                 for ii = 1:size(size_box, 1)
%                     for jj = 1:size(size_elem, 1)
%                         % apply the random field to the ocean bottom
%                         % downsampling the random field
%                         index_used = size_box(ii, 1) / MIN_SIZE_ELEM / 2;
%                         downsampling_factor = size_elem(jj, 1) / MIN_SIZE_ELEM;
%                         index_center = ceil(size(Hxy, 1) / 2);
%                         index_used = ...
%                             (-index_used:downsampling_factor:index_used)...
%                             + index_center;
%                         bath = Hxy(index_used, index_used);
% 
%                         for kk = 1:length(model_name)
%                             fkmodel = makefkmodel(model_name{kk});
%                             try
%                                 fname_in = sprintf(['Mesh_par_file_SIZE-%s_ELEM-%d_NPROCS' ...
%                                     '_X-%d_NPROCS_Y-%d_MODEL-%s'], size_box_name{ii}, ...
%                                     size_elem(jj, 1), 2, 4, model_name{kk});
%                                 meshparams3d = loadmeshparfile3d(fullfile(meshdir, fname_in));
%                             catch
%                                 continue
%                             end
%                             layers = zeros(meshparams3d.NREGIONS, 1);
%                             
%                             % assumes flat interfaces except the ocean bottom
%                             for mm = fkmodel.nlayers:-1:1
%                                 ii_itfs = fkmodel.nlayers - mm + 1;
%                                 itfs{ii_itfs} = struct(...
%                                     'SUPPRESS_UTM_PROJECTION'   , true  , ...
%                                     'NXI'                       , 2   , ...
%                                     'NETA'                      , 2  , ...
%                                     'LON_MIN'   , meshparams3d.LONGITUDE_MIN, ...
%                                     'LAT_MIN'   , meshparams3d.LATITUDE_MIN, ...
%                                     'SPACING_XI', meshparams3d.LONGITUDE_MAX - meshparams3d.LONGITUDE_MIN, ...
%                                     'SPACING_ETA', meshparams3d.LATITUDE_MAX - meshparams3d.LATITUDE_MIN, ...
%                                     'FILE', fullfile(subdir, ...
%                                         sprintf(['interf_SIZE-%s_ELEM-%d_MODEL' ...
%                                         '-%s_%02d.dat'], size_box_name{ii}, ...
%                                         size_elem(jj, 1), model_name{kk}, ii_itfs)), ...
%                                     'Z', ones(2,2) * fkmodel.layers{mm}.ztop ...
%                                     );
%                                 layers(ii_itfs) = meshparams3d.REGIONS{mm}.NZ_END - ...
%                                     meshparams3d.REGIONS{mm}.NZ_BEGIN + 1;
%                             end
% 
%                             % apply Hann Filter
%                             NXI = meshparams3d.NEX_XI + 1;
%                             NETA = meshparams3d.NEX_ETA + 1;
% 
%                             wY = shanning(NETA, 0.15);
%                             wX = shanning(NXI, 0.15);
%                             bath_taper = (wY * wX') .* bath;
%                             
%                             % apply random field to ocean bottom
%                             itfs{fkmodel.nlayers - 1}.NXI = NXI;
%                             itfs{fkmodel.nlayers - 1}.NETA = NETA;
%                             itfs{fkmodel.nlayers - 1}.SPACING_XI = size_elem(jj, 1);
%                             itfs{fkmodel.nlayers - 1}.SPACING_ETA = size_elem(jj, 2);
%                             itfs{fkmodel.nlayers - 1}.Z = bath_taper + ...
%                                 fkmodel.layers{2}.ztop;
%                             
%                             % lower the interfaces underneath to prevent
%                             % interface intersection or the layer too thin
%                             for mm = (fkmodel.nlayers - 2):-1:1
%                                 mm_fkmodel_top = fkmodel.nlayers - mm;
%                                 mm_fkmodel_bottom = fkmodel.nlayers - mm + 1;
%                                 layer_thickness = fkmodel.layers{mm_fkmodel_top}.ztop - ...
%                                     fkmodel.layers{mm_fkmodel_bottom}.ztop;
%                                 min_layer_thickness = 0.6 * layer_thickness;
%                                 if min(itfs{mm+1}.Z, [], 'all') - itfs{mm}.Z(1,1) <= min_layer_thickness
%                                     itfs{mm}.NXI = NXI;
%                                     itfs{mm}.NETA = NETA;
%                                     itfs{mm}.SPACING_XI = size_elem(jj, 1);
%                                     itfs{mm}.SPACING_ETA = size_elem(jj, 2);
%                                     itfs{mm}.Z = min(itfs{mm+1}.Z - ...
%                                         min_layer_thickness, ...
%                                         itfs{mm}.Z(1,1));
%                                 end
%                             end
% 
%                             % write the interface file and its dependents
%                             fname = sprintf('interfaces_SIZE-%s_ELEM-%d_MODEL-%s.dat', ...
%                                 size_box_name{ii}, size_elem(jj, 1), model_name{kk});
%                             writeinterfacefiles3d(itfs, layers, fullfile(subdir, fname));
%                             
%                             % reset the interfaces struct
%                             clear itfs
%                             
%                             % draw the bathymetry for reference
%                             if kk == 2
%                                 figure(1)
%                                 clf
%                                 set(gcf, 'Units', 'inches', ...
%                                     'Position', [0 1 7 6])
%                                 imagesc([-10000 10000], [-10000 10000], ...
%                                     bath_taper);
%                                 colormap(kelicol);
%                                 cb = colorbar( 'TickDirection', 'out', ...
%                                     'FontSize', 12, 'Ticks', ...
%                                     (-4:4) * round(sqrt(s2(xx)), 1, ...
%                                     'significant'));
%                                 cb.Label.String = 'elevation (m)';
%                                 grid on
%                                 titlestring = sprintf(['\\sigma^2 = %d' ...
%                                     ' | \\nu = %.1f | \\rho = %d | ' ...
%                                     '\\Deltax = %d'], s2(xx), nu(yy), ...
%                                     rho(zz), size_elem(jj, 1));
%                                 title(titlestring)
%                                 xlabel('easting (m)')
%                                 ylabel('northing (m)')
%                                 set(gca, 'TickDir', 'out', ...
%                                     'DataAspectRatio', [1 1 1], ...
%                                     'XTick', (-10:2:10) * 1e3, ...
%                                     'YTick', (-10:2:10) * 1e3, ...
%                                     'CLim', [-4 4] * round(sqrt(s2(xx)) ...
%                                     , 1, 'significant'), ...
%                                     'FontSize', 12)
%                                 set(gcf, 'Renderer', 'painters')
%                                 
%                                 figname = sprintf(['%s_RANDOM_S2-%d_' ...
%                                     'NU-%.2f_RHO-%d_ELEM-%d_IT-%03d'], ...
%                                     mfilename, s2(xx), nu(yy), rho(zz), ...
%                                     size_elem(jj, 1), nn);
%                                 figname = replace(figname, '.', 'p');
%                                 figdisp(figname, [], [], 2, [], 'epstopdf');
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
end
