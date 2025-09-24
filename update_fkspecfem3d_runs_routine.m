function update_fkspecfem3d_runs_routine(obs_struct, synmasterdir, specmasterdir, n, min_snr, min_gcarc, min_depth)
% UPDATE_FKSPECFEM3D_RUNS_ROUTINE(obs_struct, synmasterdir, specmasterdir, n, min_snr, min_gcarc)
%
% Updates the FK-SPECFEM3D runs: added ocean bottom bathymetry, adjust the
% incident angle, and prepares stf_file for FK-SPECFEM3D runs from 
% Instaseis seimogram. This will override the existing stf_file.txt in the
% directories.
%
% INPUT:
% obs_struct        a struct containing
%   - snr               signal-to-noise ratio
%   - fcorners          [lower upper] corner frequencies
%   - CCmaxs            maximum correlation coefficients for
%                       [flat bath] cases
%   - metadata          SAC Headers associated to the obsfile
%   - presiduals        InstaSeis arrival - TauP prediction for first P
%                       arrival
% synmasterdir      the master directory to the synthetic files sorted into
%                   IRIS event ID folders
% specmasterdir     the master directory to the input folders from
%                   FK-SPECFEM3D run for fluid-solid setting. Must have the
%                   following subdirectories
%   - LAYERED_OC_MERMAID_HIGHCC_**     ** is 01 to n  (see below)
%   - LAYERED_OC_MERMAID_LOWCC_**      ** is 01 to n  (see below)
% n                 number of runs for HIGHCC or LOWCC  [default: 15]
% min_snr           Signal-to-noise ratio cut-off       [default: 0]
% min_gcarc         Epicentral distance cut-off         [default: 0]
% min_depth         Minimum event depth                 [default: 0]
% 
% SEE ALSO:
% ADDOCEANBOTTOM, STFMAKER
% 
% Last modified by sirawich-at-princeton.edu, 07/16/2025

defval('specmasterdir', ...
    fullfile(getenv('REMOTE3D'), '20250717_MERMAID_INSTASEIS'))
defval('n', 15)
defval('min_snr', 0)
defval('min_gcarc', 0)
defval('min_depth', 0)

wh_valid = and(and(obs_struct.snr(:,2) > min_snr, ...
    obs_struct.metadata.GCARC > min_gcarc), ...
    obs_struct.metadata.EVDP > min_depth);
ic_list = indeks((1:length(obs_struct.CCmaxs(:,2)))', wh_valid);

% sort the observations by decreasing correlation coefficient
[~, ic] = sort(obs_struct.CCmaxs(wh_valid,2), 'descend');
ic = ic_list(ic);
% loop through the runs
for ii = 1:n
    % directory
    ddir = fullfile(specmasterdir, ...
        sprintf('LAYERED_OC_MERMAID_HIGHCC_%02d', ii));

    % get ocean bathymetry
    [~, ~, zz] = bathymetryprofile2d([20000 20000], [81 81], ...
            [obs_struct.metadata.STLO(ic(ii)) ...
            obs_struct.metadata.STLA(ic(ii))], ...
            obs_struct.metadata.BAZ(ic(ii))+180);
    zz = double(zz');

    % median elevation
    z0 = median([zz(:,1); zz(:,end); zz(1,:)'; zz(end,:)']);

    % update ocean bottom
    addoceanbottom(ddir, zz - z0, z0, true);

    % fix the elevation of the mermaid
    % the elevation of the virtual OBS has been adjusted in addoceanbottom
    [~, name, network, x, y, z] = readstations3d(fullfile(ddir, 'DATA', ...
        'STATIONS'));
    z(2) = -obs_struct.metadata.STDP(ic(ii));
    stations = struct('name', {name}, 'network', {network}, 'lat', y, ...
        'lon', x, 'elev', z, 'z', z);
    writestations3d(stations, fullfile(ddir, 'DATA', 'STATIONS'))

    % update incident angle FKMODEL file
    fkmodel = loadfkmodel(fullfile(ddir, 'DATA', 'FKMODEL'));
    fkmodel.theta = asin(obs_struct.metadata.USER9(ic(ii)) * ...
        fkmodel.layers{fkmodel.nlayers}.vp / (6371000 - 20000)) * 180/pi;
    fkmodel.twindow = 100;
    writefkmodel(fkmodel, fullfile(ddir, 'DATA', 'FKMODEL'));
    
    % identify which Instaseis syntheteic seismogram to read
    stationid = indeks(obs_struct.metadata.KSTNM{ic(ii)}, '2:5');
    eventid = obs_struct.metadata.USER7(ic(ii));
    synfile = cindeks(ls2cell(sprintf('%s%d/*_%s_*.sac', synmasterdir, ...
        eventid, stationid), 1), 1);
    
    % read the Instaseis synthetic vertical displacement seismogram
    [seis_s, hdr_s, ~, ~, tims_s] = readsac(synfile);
    tims_s = tims_s';

    % adjust the time based on the ray-theory arrival time of first phase
    t = tims_s - hdr_s.T0;

    % compute the tapering window for FK-SPECFEM3D injection
    seis_s = seis_s .* stf_window(t);

    % write the source-time function from 200 s before to 200 s after the
    % ray-theory arrival time
    wh = and(t >= -200, t<= 200);
    writetimeseries(t(wh), seis_s(wh), fullfile(ddir, 'stf_file.txt'));
end

% sort the observations by increasing correlation coefficient
[~, ic] = sort(obs_struct.CCmaxs(wh_valid,2), 'ascend');
ic = ic_list(ic);
for ii = 1:n
    % directory
    ddir = fullfile(specmasterdir, ...
        sprintf('LAYERED_OC_MERMAID_LOWCC_%02d', ii));

    % get ocean bathymetry
    [~, ~, zz] = bathymetryprofile2d([20000 20000], [81 81], ...
            [obs_struct.metadata.STLO(ic(ii)) ...
            obs_struct.metadata.STLA(ic(ii))], ...
            obs_struct.metadata.BAZ(ic(ii))+180);
    zz = double(zz');

    % median elevation
    z0 = median([zz(:,1); zz(:,end); zz(1,:)'; zz(end,:)']);

    % update ocean bottom
    addoceanbottom(ddir, zz - z0, z0, true);

    % fix the elevation of the mermaid
    % the elevation of the virtual OBS has been adjusted in addoceanbottom
    [~, name, network, x, y, z] = readstations3d(fullfile(ddir, 'DATA', ...
        'STATIONS'));
    z(2) = -obs_struct.metadata.STDP(ic(ii));
    stations = struct('name', {name}, 'network', {network}, 'lat', y, ...
        'lon', x, 'elev', z, 'z', z);
    writestations3d(stations, fullfile(ddir, 'DATA', 'STATIONS'))

    % update incident angle FKMODEL file
    fkmodel = loadfkmodel(fullfile(ddir, 'DATA', 'FKMODEL'));
    fkmodel.theta = asin(obs_struct.metadata.USER9(ic(ii)) * ...
        fkmodel.layers{fkmodel.nlayers}.vp / (6371000 - 20000)) * 180/pi;
    fkmodel.twindow = 100;
    writefkmodel(fkmodel, fullfile(ddir, 'DATA', 'FKMODEL'));
    
    % identify which Instaseis syntheteic seismogram to read
    stationid = indeks(obs_struct.metadata.KSTNM{ic(ii)}, '2:5');
    eventid = obs_struct.metadata.USER7(ic(ii));
    synfile = cindeks(ls2cell(sprintf('%s%d/*_%s_*.sac', synmasterdir, ...
        eventid, stationid), 1), 1);
    
    % read the Instaseis synthetic vertical displacement seismogram
    [seis_s, hdr_s, ~, ~, tims_s] = readsac(synfile);
    tims_s = tims_s';

    % adjust the time based on the ray-theory arrival time of first phase
    t = tims_s - hdr_s.T0 - hdr_s.B;

    % compute the tapering window for FK-SPECFEM3D injection
    seis_s = seis_s .* stf_window(t);

    % write the source-time function from 200 s before to 200 s after the
    % ray-theory arrival time
    wh = and(t >= -200, t<= 200);
    writetimeseries(t(wh), seis_s(wh), fullfile(ddir, 'stf_file.txt'));
end
end