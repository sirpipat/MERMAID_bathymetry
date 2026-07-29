function stations = makestations3d(s)
% stations = MAKESTATIONS           % makes an example of generic structure
% stations = MAKESTATIONS3D(s)
%
% Makes a struct of stations for SPECFEM3D cartesian.
%
% INPUT: 
% s                 a struct containing the following fields
%       n               number of stations in a network
%       network         network code
%       xbeg            x-coordinate of the first station
%       xend            x-coordinate of the last station
%       ybeg            y-coordinate of the first station
%       yend            y-coordinate of the last station
%       zbeg            z-coordinate of the first station
%       zend            z-coordinate of the last station
%       ebeg            elevation at the first station [default: 0]
%       eend            elevation at the last station  [default: 0]
%
% OUTPUT:
% stations          a struct containing following fields
%     network           station network names
%     name              station names (network name + 3-digit index)
%     x                 x-coordinates or longitudes
%     y                 y-coordinates or latitudes
%     elev              elevations
%     z                 z-coordinates
%
% SEE ALSO:
% LOADSTATIONS3D, WRITESTATIONS3D, READSTATIONS3D
%
% Last modified by spipatprathanporn@ucsd.edu, 07/29/2026

% Example: fluid-solid simulation with an OBS and a hydrophone float
if nargin == 0
    stations.network = {'AA'; 'MH'};
    stations.name    = {'OBS01'; 'P0009'};
    stations.x       = [0; 0];
    stations.y       = [0; 0];
    stations.elev    = [0; 0];
    stations.z       = [-5000; -1500];
    return
elseif nargin == 1
    % convert to a struct of array in case the input is an array of struct
    if numel(s) > 1
        s = array2struct(s);
    end
    if ~isfield(s, 'ebeg')
        s.ebeg = zeros(size(s.n));
    end
    if ~isfield(s, 'eend')
        s.eend = zeros(size(s.n));
    end
    % total number of stations
    N = sum(s.n);
    % cumulative number of stations have been added
    cum_n = 0;

    % allocate space for the output
    stations.network = cell(N,1);
    stations.name    = cell(N,1);
    stations.x       = nan(N,1);
    stations.y       = nan(N,1);
    stations.elev    = nan(N,1);
    stations.z       = nan(N,1);
else
    error('Invalid number of inputs')
end

% appending station entries
for ii = 1:numel(s.n)
    ii_begin = cum_n + 1;
    ii_end   = cum_n + s.n(ii);
    
    stations.network(ii_begin:ii_end) = repmat({s.network{ii}}, [s.n(ii) 1]);
    stations.name(ii_begin:ii_end) = strcat(stations.network(ii_begin:ii_end), ...
        cellstr(num2str((1:s.n(ii))', '%03d')));
    stations.x(ii_begin:ii_end) = linspace(s.xbeg(ii), s.xend(ii), s.n(ii))';
    stations.y(ii_begin:ii_end) = linspace(s.ybeg(ii), s.yend(ii), s.n(ii))';
    stations.elev(ii_begin:ii_end) = linspace(s.ebeg(ii), s.eend(ii), s.n(ii))';
    stations.z(ii_begin:ii_end) = linspace(s.zbeg(ii), s.zend(ii), s.n(ii))';

    cum_n = cum_n + s.n(ii);
end
end