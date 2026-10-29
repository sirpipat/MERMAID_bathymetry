function addlayer(ddir, n, thick, prop, fname, dz)
% ADDLAYER(ddir, n, thick, prop, fname, dz)
%
% Adds a layer to the FK-SPECFEM3D_Cartesian model. It assumes the mesh has
% the same layers as the outside layers assumed by FKMODEL file.
%
% INPUT:
% ddir      directory to a FK-SPECFEM3D run
% n         underneath which layer from top to insert a layer
% thick
% prop      material properties of the added layer
% - rho         density
% - vp          P-wave velocity
% - vs          S-wave velocity
% - Q_Kappa     bulk quality factor
% - Q_mu        shear quality factor
% fname     name of the topography file
% dz        vertical element size of the added layer
%
% Last modified by spipatprathanporn@ucsd.edu, 10/29/2026

% default topography file name for the added layer
defval('fname', sprintf('interf_newlayer_%s.dat', ...
    indeks(hash(now,'MD5'), 27:32)))

% determine the file names and read parameters
parfile = fullfile(ddir, 'DATA', 'Par_file');
params = loadparfile3d(parfile);

fkfile = fullfile(ddir, params.FKMODEL_FILE);
fkmodel = loadfkmodel(fkfile);

meshparfile = fullfile(ddir, 'DATA', 'meshfem3D_files', 'Mesh_Par_file');
meshparams = loadmeshparfile3d(meshparfile);

interfacefile = fullfile(ddir, 'DATA', 'meshfem3D_files', ...
    meshparams.INTERFACES_FILE);
[itfs, layers] = loadinterfacefiles3d(interfacefile);

% verify the number of layers above
if n < 0 || n > fkmodel.nlayers
    error('Cannot add a new layer. n must be between 0 and %d inclusive\n', ...
        fkmodel.nlayers);
end

% default vertical element size of the inserted layer
defval('dz', (meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN) / ...
    meshparams.NEX_XI);

% define the number of layers above and below the insertion
numlayers_old = fkmodel.nlayers;
numlayers_below = fkmodel.nlayers - n;
numlayers_above = n;

% the number of elements in Z-axis for the inserted layer
nznew = max(1, floor(thick / dz));

%% insert layer to fkmodel
% update bottom of the box if the shifted interfaces go below the bottom of
% the block
if -(min(itfs{1}.Z,[],'all') - thick) >= meshparams.DEPTH_BLOCK_KM * 1000
    to_remesh = false;
    meshparams.DEPTH_BLOCK_KM = meshparams.DEPTH_BLOCK_KM + thick / 1000;
    fkmodel.origin_wavefront(3) = -meshparams.DEPTH_BLOCK_KM * 1000;
else
    to_remesh = true;
    thickness_old = nan(numlayers_old, 1);
    for ii = 1:numlayers_old-1
        thickness_old(ii) = fkmodel.layers{ii}.ztop - ...
            fkmodel.layers{ii+1}.ztop;
    end
    thickness_old(numlayers_old) = fkmodel.layers{numlayers_old}.ztop - ...
        (-meshparams.DEPTH_BLOCK_KM * 1000);
    dzs_old = thickness_old ./ flipud(layers);
    dzs = [dzs_old(1:numlayers_above); dz; dzs_old(numlayers_above+1:end)];
end

% make a new layer of fkmodel
fknewlayer = struct('rho', prop.rho, 'vp', prop.vp, 'vs', prop.vs, ...
    'ztop', fkmodel.layers{n+1}.ztop);
% insert new layer
fkmodel.layers = [fkmodel.layers(1:numlayers_above); fknewlayer; ...
    fkmodel.layers(numlayers_above+1:end)];
fkmodel.nlayers = fkmodel.nlayers + 1;
% update ztop
for ii = (numlayers_above+1) + (1:numlayers_below)
    fkmodel.layers{ii}.ztop = fkmodel.layers{ii}.ztop - thick;
end


%% insert interface
layers = [layers(1:numlayers_below); nznew; layers(numlayers_below+1:end)];
itfn = struct( ...
    'SUPPRESS_UTM_PROJECTION', 1, ...
    'NXI', 2, ...
    'NETA', 2, ...
    'LON_MIN', meshparams.LONGITUDE_MIN, ...
    'LAT_MIN', meshparams.LATITUDE_MIN, ...
    'SPACING_XI', meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN, ...
    'SPACING_ETA', meshparams.LATITUDE_MAX - meshparams.LATITUDE_MIN, ...
    'FILE', fullfile(ddir, 'DATA', 'meshfem3D_files', fname), ...
    'Z', repmat(fkmodel.layers{numlayers_above+1}.ztop, [2 2]));
itfs = [itfs(1:numlayers_below); itfn; itfs(numlayers_below+1:end)];
% shift the interfaces down by the thickness of the inserted layer
for ii = 1:numlayers_below
    itfs{ii}.Z = itfs{ii}.Z - thick;
end

%% modify the mesh parametrs
% fill in missing fields required for MATERIAL used in mesh parameters
if ~isfield(prop, 'Q_Kappa')
    prop.Q_Kappa = 9999;
end
if ~isfield(prop, 'Q_mu')
    prop.Q_mu = 9999;
end
if ~isfield(prop, 'anisotropy_flag')
    prop.anisotropy_flag = 0;
end
if prop.vs > 0
    prop.domain_id = 2;
else
    prop.domain_id = 1;
end
% insert the material of the new layer to the list of materials
meshparams.NMATERIALS = meshparams.NMATERIALS + 1;
meshparams.MATERIALS = [meshparams.MATERIALS(1:numlayers_above); ...
    prop; meshparams.MATERIALS(numlayers_above+1:end)];

% make a new region
region = struct('NEX_XI_BEGIN', 1, ...
                'NEX_XI_END', meshparams.NEX_XI, ...
                'NEX_ETA_BEGIN', 1, ...
                'NEX_ETA_END', meshparams.NEX_ETA, ...
                'NZ_BEGIN', sum(layers(1:numlayers_below)) + 1, ...
                'NZ_END', sum(layers(1:numlayers_below)) + nznew, ...
                'material_id', numlayers_above+1);
% insert a new region
% assume the regions are ordered from top to bottom and span all XI,ETA
meshparams.NREGIONS = meshparams.NREGIONS + 1;
meshparams.REGIONS = [meshparams.REGIONS(1:numlayers_above); region; ...
    meshparams.REGIONS(numlayers_above+1:end)];
% shift NZ_BEGIN and NZ_End of the layers above
for ii = 1:numlayers_above
    meshparams.REGIONS{ii}.NZ_BEGIN = meshparams.REGIONS{ii}.NZ_BEGIN + ...
        nznew;
    meshparams.REGIONS{ii}.NZ_END = meshparams.REGIONS{ii}.NZ_END + nznew;
end
% update material ID
for ii = (numlayers_above+1) + (1:numlayers_below)
    meshparams.REGIONS{ii}.material_id = ii;
end

%% save the parameters
writefkmodel(fkmodel, fkfile);
writemeshparfile3d(meshparams, meshparfile);
writeinterfacefiles3d(itfs, layers, interfacefile);

%% remesh if necessary
if to_remesh
    dx = (meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN) / ...
        meshparams.NEX_XI;
    remesh3d(ddir, dx, dzs);
end
end