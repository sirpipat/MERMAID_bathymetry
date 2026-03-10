function varargout = setboxsize(ddir, boxsize, elemsize)
% [elemsize, numelems] = SETBOXSIZE(ddir, [x y z], [dx dy dz])
%
% Sets the box size of an FK-SPECFEM3D run. It also adjusts the number of
% elements, element sizes, and origin wavefront. It may removes layers and
% interfaces if the interfaces are below or intesect the bottom of the box.
%
% INPUT:
% ddir          directory to an FK-SPECFEM3D run
% [x y z]       box dimension [lon lat depth] in meters
% [dx dy dz]    element size [lon lat depth] in meters
%
% OUTPUT:
% elemsize      adjusted element size in meters
% numelems      number of elements in [lon lat depth] dimension
%
% SEE ALSO:
% SPECFEM3D_INPUT_SETUP, SPECFEM3D_INPUT_SETUP_FLAT
%
% Last modified by sirawich-at-princeton.edu, 03/09/2026

% what to change
% v meshparams.LATITUDE_MIN/MAX
% v meshparams.LONGITUDE_MIN/MAX
% v meshparams.DEPTH_BLOCK_KM
% v meshparams.NEX_XI/ETA
% v meshparams.NEX_XI/ETA_END
% v meshparams.REGIONS{mm}.NZ_BEGIN/END
% v itfs{ii}.LON_MIN
% v ifts{ii}.LAT_MIN
% v itfs{ii}.SPACING_XI
% v itfs{ii}.SPACING_ETA
% v layers(ii)
% v fkmodel.origin_wavefront --> [0 meshparams.LATITIUDE_MIN -boxsize(3)]

defval('boxsize', [20000 20000 20000])
defval('elemsize', [-1 -1 -1])

if isscalar(boxsize)
    boxsize = [1 1 1] * boxsize;
end
if isscalar(elemsize)
    elemsize = [1 1 1] * elemsize;
end

%% list of files to read and write
meshfile = fullfile(ddir, 'DATA', 'meshfem3D_files', 'Mesh_Par_file');
interfacefile = fullfile(ddir, 'DATA', 'meshfem3D_files', 'interfaces.dat');
fkmodelfile = fullfile(ddir, 'DATA', 'FKMODEL');

%% read the parameters
meshparams = loadmeshparfile3d(meshfile);
[itfs, layers] = loadinterfacefiles3d(interfacefile);
fkmodel = loadfkmodel(fkmodelfile);

%% check the input conditions
% elemsize
if elemsize(1) <= 0
    elemsize(1) = (meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN) / meshparams.NEX_ETA;
end
if elemsize(2) <= 0
    elemsize(2) = (meshparams.LATITUDE_MAX - meshparams.LATITUDE_MIN) / meshparams.NEX_XI;
end
if elemsize(3) <= 0
    elemsize(3) = (meshparams.DEPTH_BLOCK_KM * 1000) / sum(layers);
end

% Round the number of elements to be 8 times the multiple of the number of 
% processors and adjust the element size accordingly.
% See https://specfem3d.readthedocs.io/en/latest/03_mesh_generation
nex_eta_exact = boxsize(1) / elemsize(1);
nex_eta_round = round(nex_eta_exact / (8 * meshparams.NPROC_ETA)) * ...
    (8 * meshparams.NPROC_ETA);
elemsize(1) = boxsize(1) / nex_eta_round;

nex_xi_exact = boxsize(2) / elemsize(2);
nex_xi_round = round(nex_xi_exact / (8 * meshparams.NPROC_XI)) * ...
    (8 * meshparams.NPROC_XI);
elemsize(2) = boxsize(2) / nex_xi_round;

nez_exact = boxsize(3) / elemsize(3);
nez_round = round(nez_exact);
elemsize(3) = boxsize(3) / nez_round;

numelems = [nex_eta_round nex_xi_round nez_round];

fprintf('The element size is adjusted to %.2f x %.2f x %.2f m.\n', ...
    elemsize(1), elemsize(2), elemsize(3));
str = sprintf('%d', nex_eta_round * nex_xi_round * nez_round);
fprintf('The number of elements is %d x %d x %d = %s elements.\n', ...
    nex_eta_round, nex_xi_round, nez_round, ...
    insert(str, ',', flip(length(str)-2:-3:2)));

% fkmodel layers ztop
validlayers = false([fkmodel.nlayers 1]);
for ii = 1:fkmodel.nlayers
    validlayers(ii) = (fkmodel.layers{ii}.ztop > -boxsize(3));
end

% discard any interfaces that intersect the bottom of the box
for ii = 1:fkmodel.nlayers
    if min(itfs{fkmodel.nlayers+1-ii}.Z, [], 'all') <= -boxsize(3)
        validlayers(ii) = false;
    end
end
fkmodel.nlayers = sum(validlayers);
fkmodel.layers = fkmodel.layers(validlayers);
meshparams.NREGIONS = sum(validlayers);
meshparams.REGIONS = meshparams.REGIONS(validlayers);
layers = layers(flip(validlayers));
itfs = itfs(flip(validlayers));

%% modify values
meshparams.LATITUDE_MIN   = -boxsize(2) / 2;
meshparams.LATITUDE_MAX   =  boxsize(2) / 2;
meshparams.LONGITUDE_MIN  = -boxsize(1) / 2;
meshparams.LONGITUDE_MAX  =  boxsize(1) / 2;
meshparams.DEPTH_BLOCK_KM =  boxsize(3) / 1000;
meshparams.NEX_XI         =  nex_xi_round;
meshparams.NEX_ETA        =  nex_eta_round;

% compute the number of elements in z-dimension for each remaining layer
layers(1) = round((boxsize(3) + fkmodel.layers{fkmodel.nlayers}.ztop) / ...
    elemsize(3));
for ii = 2:fkmodel.nlayers
    layers(ii) = round((boxsize(3) + ...
        fkmodel.layers{fkmodel.nlayers+1-ii}.ztop) / elemsize(3)) - ...
        layers(ii-1);
end

% update the begin/end of element number in x, y, and z directions
for ii = length(meshparams.REGIONS):-1:1
    meshparams.REGIONS{ii}.NEX_XI_END = meshparams.NEX_XI;
    meshparams.REGIONS{ii}.NEX_ETA_END = meshparams.NEX_ETA;
    if ii == length(meshparams.REGIONS)
        meshparams.REGIONS{ii}.NZ_BEGIN = 1;
        meshparams.REGIONS{ii}.NZ_END = layers(length(meshparams.REGIONS)+1-ii);
    else
        meshparams.REGIONS{ii}.NZ_BEGIN = meshparams.REGIONS{ii+1}.NZ_END + 1;
        meshparams.REGIONS{ii}.NZ_END = meshparams.REGIONS{ii+1}.NZ_END + layers(length(meshparams.REGIONS)+1-ii);
    end
end

% update the interfaces
for ii = 1:length(itfs)
    itfs{ii}.LAT_MIN = meshparams.LATITUDE_MIN;
    itfs{ii}.LON_MIN = meshparams.LONGITUDE_MIN;
    itfs{ii}.SPACING_XI = boxsize(2) / (itfs{ii}.NXI - 1);
    itfs{ii}.SPACING_ETA = boxsize(1) / (itfs{ii}.NETA - 1);
end

% move the origin wavefront
fkmodel.origin_wavefront = [0 meshparams.LATITUDE_MIN -boxsize(3)];

%% write the parameters to the file
writemeshparfile3d(meshparams, meshfile);
writeinterfacefiles3d(itfs, layers, interfacefile);
writefkmodel(fkmodel, fkmodelfile);

%% function return values
varns = {elemsize, numelems};
varargout = varns(1:nargout);
end