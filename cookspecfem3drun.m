function cookspecfem3drun(ddir, model, boxsize, elemsize, stftype, freq, theta, fs, nprocs, bath, s2, nu, rho, it)
% COOKSPECFEM3DRUN(ddir, model, boxsize, elemsize, stftype, freq, ...
%     theta, fs, nprocs, bath, s2, nu, rho, it)
%
% Makes a FK-SPECFEM3D run from stockfiles.
%
% SEE ALSO:
% MAKESTOCKFILE3D
%
% Last modified by sirawich-at-princeton.edu, 08/19/2026

if strcmp(ddir, 'options')
    % print out options
    return
end

defval('model', 'PS2024')
defval('boxsize', 'thick')
defval('elemsize', 250)
defval('stftype', 'Gaussian')
defval('freq', 4)
defval('theta', 0)
defval('fs', 200)
defval('nprocs', [2 4])
defval('bath', 'FLAT')
defval('s2', 100)
defval('nu', 0.5)
defval('rho', 1000)
defval('it', 1)

% input validation
if ~(strcmp(model, 'PS2024') || strcmp(model, 'OC') || ...
        strcmp(model, 'OCM') || strcmp(model, 'OSCM') || ...
        strcmp(model, 'Crust1'))
    error("model must be either 'PS2024', 'OC', 'OCM', 'OSCM', or 'Crust1'")
end
if ~(strcmp(boxsize, 'original') || strcmp(boxsize, 'thick'))
    error("boxsize must be either 'original' or 'thick'")
end
if ~any(elemsize == [100 250 500])
    error('elemsize must be eigther 100, 250, or 500')
end
if ~(strcmp(stftype, 'Gaussian') || strcmp(stftype, 'Ricker'))
    error("stftype must be either 'Gaussian' or 'Ricker'")
end
if ~any(freq == [0.1 0.5 1 2 4 10])
    error('freq must be from [0.1 0.5 1 2 4 10]')
end
if ~any(theta == [0 10 40])
    error('theta must be either 0, 10, or 40')
end
if ~any(fs == [100 200 400 1000])
    error('fs must be from [100 200 400 1000]')
end
if ~(all(nprocs == [2 4]) || all(nprocs == [4 4]) || all(nprocs == [4 8]))
    error('nprocs must be either [2 4], [4 4], or [4 8]')
end
if ~(strcmp(bath, 'FLAT') || strcmp(bath, 'RANDOM'))
    error("bath must be either 'FLAT' or 'RANDOM'")
end
if strcmp(bath, 'RANDOM') && ~any(s2 == [100 1000 10000 100000])
    error("s2 must be from [100 1000 10000 100000]")
end
if strcmp(bath, 'RANDOM') && ~any(nu == [0.5 1.0 1.5])
    error('nu must be from [0.5 1.0 1.5]')
end
if strcmp(bath, 'RANDOM') && ~any(rho == [1000 2000 4000 8000])
    error('rho must be from [1000 2000 4000 8000]')
end
if strcmp(bath, 'RANDOM') && (it < 1 || it > 100)
    error('1 <= it <= 100')
end

%% STOCKFILES directories
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

%% COPY
system(sprintf('mkdir -p %s', ddir));
dataddir = fullfile(ddir, 'DATA');
system(sprintf('mkdir -p %s', dataddir));
meshddir = fullfile(ddir, 'DATA', 'meshfem3D_files');
system(sprintf('mkdir -p %s', meshddir));

% FKMODEL
source = sprintf('FKMODEL_%s_THETA-%02d_FREQ-%.2f', model, theta, freq);
source = fullfile(fkdir, replace(source, '.', 'p'));
target = fullfile(dataddir, 'FKMODEL');
system(sprintf('cp %s %s', source, target));

% STF File
source = replace(sprintf('stf_file_%s_%.2f', stftype, freq), '.', 'p');
source = fullfile(stfdir, [source '.txt']);
target = fullfile(ddir, 'stf_file.txt');
system(sprintf('cp %s %s', source, target));

% CMTSOLUTION file
system(sprintf('cp %s %s', fullfile(cmtdir, 'CMTSOLUTION'), dataddir));

% Par file
source = fullfile(pardir, sprintf('Par_file_FS-%d', fs));
target = fullfile(dataddir, 'Par_file');
system(sprintf('cp %s %s', source, target));

% Mesh par file
source = sprintf(['Mesh_par_file_SIZE-%s_ELEM-%d_NPROCS_X-%d_NPROCS_Y' ...
    '-%d_MODEL-%s'], boxsize, elemsize, nprocs(1), nprocs(2), model);
source = fullfile(meshdir, source);
target = fullfile(meshddir, 'Mesh_Par_file');
system(sprintf('cp %s %s', source, target));

% Interface file
if strcmp(bath, 'FLAT')
    source = sprintf('interfaces_SIZE-%s_ELEM-%d_MODEL-%s.dat', ...
        boxsize, elemsize, model);
    source = fullfile(interfdir, 'FLAT', source);
    target = fullfile(meshddir, 'interfaces.dat');
    system(sprintf('cp %s %s', source, target));
    
    % topo files
    source = sprintf('interf_SIZE-%s_ELEM-%d_MODEL-%s_*.dat', ...
        boxsize, elemsize, model);
    source = fullfile(interfdir, 'FLAT', source);
    system(sprintf('cp %s %s', source, meshddir));
elseif strcmp(bath, 'RANDOM')
    source = sprintf('interfaces_SIZE-%s_ELEM-%d_MODEL-%s.dat', ...
        boxsize, elemsize, model);
    source = fullfile(interfdir, 'RANDOM', source);
    target = fullfile(meshddir, 'interfaces.dat');
    system(sprintf('cp %s %s', source, target));
    
    % topo files
    source = sprintf('interf_SIZE-%s_ELEM-%d_MODEL-%s_*.dat', ...
        boxsize, elemsize, model);
    source = fullfile(interfdir, 'RANDOM', source);
    system(sprintf('cp %s %s', source, meshddir));
    
    % read interfaces
    [itfs, layers] = loadinterfacefiles3d(target);
    % collect flat elevation for each interface
    elev = zeros(size(layers));
    for ii = 1:length(layers)
        elev(ii) = itfs{ii}.Z(1);
    end
    
    % read random field
    if strcmp(boxsize, 'original')
        boxsizenum = 28000;
    else
        boxsizenum = 20000;
    end
    fname = sprintf(['topo_RANDOM_S2-%d_NU-%.2f_RHO-%d_SIZE-%d_ELEM-%d_' ...
        'IT-%03d.dat'], s2, nu, rho, boxsizenum, elemsize, it);
    fname = fullfile(topodir, fname);
    fid = fopen(fname, 'r');
    topo = fscanf(fid, '%f');
    fclose(fid);
    topo = reshape(topo, sqrt(length(topo)), sqrt(length(topo)))';
    
    % assign topo to the ocean bottom interface
    NETA = itfs{length(layers)-1}.NETA;
    NXI = itfs{length(layers)-1}.NXI;
    wY = shanning(NETA, 0.15);
    wX = shanning(NXI, 0.15);
    topo_taper = (wY * wX') .* topo;
    itfs{length(layers)-1}.Z = itfs{length(layers)-1}.Z(1) +  topo_taper;
    
    for ii = (length(layers)-2):-1:1
        layer_thickness = elev(ii+1) - elev(ii);
        min_layer_thickness = 0.6 * layer_thickness;
        if min(itfs{ii+1}.Z, [], 'all') - elev(ii) <= min_layer_thickness
            itfs{ii}.NXI = NXI;
            itfs{ii}.NETA = NETA;
            itfs{ii}.SPACING_XI = elemsize;
            itfs{ii}.SPACING_ETA = elemsize;
            itfs{ii}.Z = min(itfs{ii+1}.Z - min_layer_thickness, ...
                itfs{ii}.Z(1,1));
        else
            break
        end
    end
    writeinterfacefiles3d(itfs, layers, target)
    
end

% STATION
source = fullfile(stadir, 'STATIONS');
[n, name, network, x, y, z] = readstations3d(source);
if strcmp(bath, 'RANDOM')
    z(1) = itfs{length(layers)-1}.Z(ceil(NETA/2), ceil(NXI/2));
end
stations.name = name;
stations.network = network;
stations.lat = y;
stations.lon = x;
stations.elev = z;
stations.z = z;
target = fullfile(dataddir, 'STATIONS');
writestations3d(stations, target);
end
