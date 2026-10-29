function setinferface(ddir, n, zz)
% SETINTERFACE(ddir, n, zz)
%
% Sets the interface at the top of n-th layer, counting from top. It will
% updates the interface file and the FKMODEL file. Please make sure you
% have tapered the interface before setting the interface. 
%
% INPUT:
% ddir          directory of the FK SPECFEM3D run
% n             n-th layer from the top
% zz            elevation grid of the interface
%
% Last modified by spipatprathanporn@ucsd.edu, 10/29/2026

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

% old element sizes
[dx, ~, dzs] = getelementsize(ddir);

% flip them to make indexing easier
itfs = flipud(itfs);
layers = flipud(layers);

zmin = min(zz, [], 'all');
zmax = max(zz, [], 'all');
if n > 1 && zz(1,1) >= fkmodel.layers{n-1}.ztop
    error('Your interface will overlap with the interface above.')
end
if n < length(layers) && zz(1,1) <= fkmodel.layers{n+1}.ztop
    error('Your interface will overlap with the interface below.')
end
if zmin <= -meshparams.DEPTH_BLOCK_KM * 1000
    error('Your interface will overlap with the bottom of the box.')
end

% edit the FK model
fkmodel.layers{n}.ztop = zz(1,1);

% edit the inteface
itfs{n}.NXI = size(zz,2);
itfs{n}.NETA = size(zz,1);
itfs{n}.SPACING_XI = (meshparams.LONGITUDE_MAX - ...
    meshparams.LONGITUDE_MIN) / (itfs{n}.NXI-1);
itfs{n}.SPACING_ETA = (meshparams.LATITUDE_MAX - ...
    meshparams.LATITUDE_MIN) / (itfs{n}.NETA-1);
itfs{n}.Z = zz;

% save
writefkmodel(fkmodel, fkfile);
writeinterfacefiles3d(flipud(itfs), flipud(layers), interfacefile);

% remesh
remesh3d(ddir, dx, dzs);
end