function setlayer(ddir, n, thick, prop, dz)
% SETLAYER(ddir, n, thick, prop, dz)
%
% Modifies the properties of the n-th layer (from top) of the
% FK-SPECFEM3D_Cartesian model. It assumes the mesh has the same layers as
% the outside layers in the FKMODEL file.
%
% INPUT:
% ddir      directory the FK-SPECFEM3D run
% n         which layer from the top to be modified
% thick     layer thickness
% prop      material properties
% - rho         density
% - vp          P-wave velocity
% - vs          S-wave velocity
% - Q_Kappa     bulk quality factor
% - Q_mu        shear quality factor
% dz        target thickness per element
%
% SEE ALSO:
% ADDLAYER, SETSEDIMENTLAYER
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

% old element size thickness
[dx, ~, dzs] = getelementsize(ddir);

% set the fkmodel
fn = fieldnames(prop);
for ii = 1:length(fn)
    switch lower(fn{ii})
        case 'rho'
            fkmodel.layers{n}.rho = prop.(fn{ii});
        case 'vp'
            fkmodel.layers{n}.vp = prop.(fn{ii});
        case 'vs'
            fkmodel.layers{n}.vs = prop.(fn{ii});
    end
end
% adjust ztop of the layers below
if n < fkmodel.nlayers
    zbot = fkmodel.layers{n+1}.ztop;
else
    zbot = -meshparams.DEPTH_BLOCK_KM * 1000;
end
thick_old = fkmodel.layers{n} - zbot;
if ~isempty(thick) || thick_old ~= thick
    for ii = (n+1):fkmodel.nlayers
        fkmodel.layers{ii}.ztop = fkmodel.layers{ii}.ztop - ...
            (thick - thick_old);
    end
end
% adjust the elevations of the interfaces
if ~isempty(thick)
    for ii = 1:(fkmodel.nlayers - n)
        itfs{ii}.Z = itfs{ii}.Z - (thick - thick_old);
    end
    meshparams.DEPTH_BLOCK_KM = meshparams.DEPTH_BLOCK_KM + ...
        (thick - thick_old);
end

% save the parameters
writefkmodel(fkmodel, fkfile);
writemeshparfile3d(meshparams, meshparfile);
writeinterfacefiles3d(itfs, layers, interfacefile);

% remesh
if ~isempty(dz)
    dzs(n) = dz;
end
remesh3d(ddir, dx, dzs);
end
