function [dx, dy, dzs] = getelementsize(ddir)
% [dx, dy, dzs] = getelementsize(ddir)
%
% Gets element sizes in X/Y and Z direction for each layer.
%
% INPUT:
% ddir          directory the FK-SPECFEM3D run
%
% OUTPUT:
% dx            element size in X direction
% dy            element size in Y direction
% dzs           element sizes in Z direction:
%               dzs(1)      = top layer
%               ...
%               dzs(end)    = bottom layer
%
% Last modified by spipatprathanporn@ucsd.edu, 10/18/2026

% determine the file names and read parameters
parfile = fullfile(ddir, 'DATA', 'Par_file');
params = loadparfile3d(parfile);

fkfile = fullfile(ddir, params.FKMODEL_FILE);
fkmodel = loadfkmodel(fkfile);

meshparfile = fullfile(ddir, 'DATA', 'meshfem3D_files', 'Mesh_Par_file');
meshparams = loadmeshparfile3d(meshparfile);

interfacefile = fullfile(ddir, 'DATA', 'meshfem3D_files', ...
    meshparams.INTERFACES_FILE);
[~, layers] = loadinterfacefiles3d(interfacefile);

% compute the element size in X and Y directions
dx = (meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN) / ...
    meshparams.NEX_XI;
dy = (meshparams.LATITUDE_MAX - meshparams.LATITUDE_MIN) / ...
    meshparams.NEX_ETA;

% compute the element sizes in Z direction for each layer
thickness = nan(fkmodel.nlayers,1);
for ii = 1:fkmodel.nlayers-1
    thickness(ii) = fkmodel.layers{ii}.ztop - ...
        fkmodel.layers{ii+1}.ztop;
end
thickness(fkmodel.nlayers) = fkmodel.layers{fkmodel.nlayers}.ztop - ...
    (-meshparams.DEPTH_BLOCK_KM * 1000);
dzs = thickness ./ flipud(layers);
end