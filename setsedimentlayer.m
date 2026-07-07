function setsedimentlayer(ddir, thick, rho, vp, vs, isadding)
% SETSEDIMENTLAYER(ddir, thickness, density, isadding)
%
% Set or add the sediment layer to a FK-SPECFEM3D run.
%
% INPUT:
% ddir          directory to a FK-SPECFEM3D run
% thick         uniform sediment layer thickness (scalar, m)
% rho           uniform sediment density         (scalar, kg/m^3)
% vp            uniform sediment P-wave speed    (scalar, m/s)
% vs            uniform seidment S-wave speed)   (scalar, m/s)
% isadding      whether to add a sediment layer instead of modify the
%               existing sediment layer in the model [default: false]
%
% Last modified by spipatprathanporn@ucsd.edu, 07/06/2026

defval('thick', [])
defval('rho', [])
defval('vp', [])
defval('vs', [])
defval('isadding', false)
isadding = false;

% data file names
fkfile = fullfile(ddir, 'DATA', 'FKMODEL');
meshparfile = fullfile(ddir, 'DATA', 'meshfem3D_files', 'Mesh_Par_file');
interffile = fullfile(ddir, 'DATA', 'meshfem3D_files', 'interfaces.dat');

% read the data files
fkmodel = loadfkmodel(fkfile);
meshparams = loadmeshparfile3d(meshparfile);
[itfs, layers] = loadinterfacefiles3d(interffile);

% determine thickness of the layer outside of the box
fklayers = array2struct(fkmodel.layers);
fklayers.zbot = [fklayers.ztop(2:end); -meshparams.DEPTH_BLOCK_KM*1000];
fklayers.thickness = fklayers.ztop - fklayers.zbot;

% determine target element thickness
% dx = (meshparams.LONGITUDE_MAX - meshparams.LONGITUDE_MIN) / meshparams.NEX_XI;
% dz = fklayers.thickness ./ flipud(layers);

%% Part I: directly modify the fkmodel and interface
if isadding
    % TODO: Implement
else
    if ~isempty(rho)
        if rho <= 0
            error('rho must be greater than zero.')
        end
        fkmodel.layers{2}.rho = rho;
    end
    if ~isempty(vp)
        if vp < 0
            error('vp cannot be smaller than zero')
        end
        fkmodel.layers{2}.vp = vp;
    end
    if ~isempty(vs)
        if vs < 0
            error('vs cannot be smaller than zero')
        end
        fkmodel.layers{2}.vs = vs;
    end
    if ~isempty(thick)
        if thick <= 0
            error('dz must be greater than zero.')
        end
        fkmodel.layers{3}.ztop = fkmodel.layers{2}.ztop - thick;
        itfs{end-2}.Z = itfs{end-1}.Z - thick;
    end
end

writefkmodel(fkmodel, fkfile);
writeinterfacefiles3d(itfs, layers, interffile);

%% Part II: remesh
if isadding || ~isempty(thick)
    remesh3d(ddir, [], [250 200 320 350*ones(1, length(layers)-3)])
end
end