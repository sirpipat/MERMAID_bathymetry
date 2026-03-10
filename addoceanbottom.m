function addoceanbottom(ddir, topo, z0, taper)
% ADDOCEANBOTTOM(ddir, topo, z0, taper)
%
% Modifies the interface file and a topography file of the ocean bottom for
% a FK-SPECFEM3D run. It updates the depth of the virtual OBS station in 
% the station file. It may modify the topgraphy files of deeper
% interfaces to keep the layer thickness reasonable (not too thin or
% negative).
%
% INPUT:
% ddir          directory to a FK-SPECFEM3D run
% topo          2D matrix describing the deviation of bathymetry from zero
% z0            mean elevation of the topography if topo has zero mean. If
%               z0 is empty, it will use elevation from FKMODEL file by
%               default.
% taper         whether to taper to z0 elevation or not
%
% Last modified by sirawich-at-princeton.edu, 03/09/2026

fkmodel = loadfkmodel(fullfile(ddir, 'DATA', 'FKMODEL'));

defval('z0', fkmodel.layers{2}.ztop)

% process elevation
if taper
    wn = shanning2d(size(topo), 0.15, 0, 0);
    topo = topo .* wn;
end

% add mean elevation
topo = topo + z0;

% modify the interface and topography files
[itfs, layers] = loadinterfacefiles3d(fullfile(ddir, 'DATA', ...
    'meshfem3D_files', 'interfaces.dat'));

% collect flat elevation for each interface
elev = zeros(size(layers));
for ii = 1:length(layers)
    elev(ii) = itfs{ii}.Z(1);
end

size_xi = itfs{end}.SPACING_XI;
size_eta = itfs{end}.SPACING_ETA;
spacing_xi = size_xi / (size(topo,1) - 1);
spacing_eta = size_eta / (size(topo,2) - 1);
itfs{end-1}.NXI = size(topo, 1);
itfs{end-1}.NETA = size(topo, 2);
itfs{end-1}.SPACING_XI = spacing_xi;
itfs{end-1}.SPACING_ETA = spacing_eta;
itfs{end-1}.Z = topo;

% modifies the topgraphy files of deeper interfaces to keep the layer 
% thickness reasonable (not too thin or negative).
for ii = (length(layers)-2):-1:1
    layer_thickness = elev(ii+1) - elev(ii);
    min_layer_thickness = 0.6 * layer_thickness;
    if min(itfs{ii+1}.Z, [], 'all') - elev(ii) <= min_layer_thickness
        itfs{ii}.NXI = size(topo, 1);
        itfs{ii}.NETA = size(topo, 2);
        itfs{ii}.SPACING_XI = spacing_xi;
        itfs{ii}.SPACING_ETA = spacing_eta;
        itfs{ii}.Z = min(itfs{ii+1}.Z - min_layer_thickness, ...
            itfs{ii}.Z(1,1));
    else
        break
    end
end

% adjust parameter files if z0 is different from ztop of
% uppermost crustal layer (i.e. ocean bottom)
if z0 ~= fkmodel.layers{2}.ztop
    % adjust FK model file
    fkmodel.layers{2}.ztop = z0;
    writefkmodel(fkmodel, fullfile(ddir, 'DATA', 'FKMODEL'))

    % adjust Mesh parameter file
    meshparams = loadmeshparfile3d(fullfile(ddir, 'DATA', ...
        'meshfem3D_files', 'Mesh_Par_file'));
    box_depth_m = meshparams.DEPTH_BLOCK_KM * 1000;
    NZ_MAX = meshparams.REGIONS{1}.NZ_END;
    meshparams.REGIONS{2}.NZ_END = round((box_depth_m + z0) / box_depth_m  * NZ_MAX);
    meshparams.REGIONS{1}.NZ_BEGIN = meshparams.REGIONS{2}.NZ_END + 1;
    writemeshparfile3d(meshparams, fullfile(ddir, 'DATA', ...
        'meshfem3D_files', 'Mesh_Par_file'));
    
    % adjust layers in interface file
    layers(end) = meshparams.REGIONS{1}.NZ_END - meshparams.REGIONS{2}.NZ_END;
    layers(end-1) = meshparams.REGIONS{2}.NZ_END - ...
        meshparams.REGIONS{2}.NZ_BEGIN + 1;
end

writeinterfacefiles3d(itfs, layers, fullfile(ddir, 'DATA', ...
    'meshfem3D_files', 'interfaces.dat'))

% fix the elevation of the station
[~, name, network, x, y, z] = readstations3d(fullfile(ddir, 'DATA', ...
    'STATIONS'));
z(1) = topo(floor((size(topo,1)+1)/2), floor((size(topo,2)+1)/2));
stations = struct('name', {name}, 'network', {network}, 'lat', y, ...
    'lon', x, 'elev', z, 'z', z);
writestations3d(stations, fullfile(ddir, 'DATA', 'STATIONS'))
end