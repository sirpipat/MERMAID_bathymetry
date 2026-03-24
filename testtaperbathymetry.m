function varargout = testtaperbathymetry(s2, nu, rho, dx, big, small, taper, plt)
% [bath_big, bath_small, bath_small_old] = ...
%     TESTTAPERBATHYMETRY(s2, nu, rho, dx, big, small, taper)
%
% Generates the tapered bathymetry grids, both the cosine tapered and the
% modified version, as well as the big bathymetry grids serving as the
% ground truth to test the effects of tapered bahtymetry on waveform
% modling.
%
% INPUT:
% s2            The first Matern parameter, aka sigma^2, in field units^2 
% nu            The second Matern parameter (smoothness)
% rho           The third Matern parameter, in the units of the grid
%               (length scale)
% dx            element size
% big           the number of sample points of the big bathymetry grid
% small         the number of sample points of the small bathymetry grid
% taper         struct with a following member variables
%       trim        width of the 4 edges that will be set to be flat
%       cosine      width of the cosine taper
%
% OUTPUT:
% bath_big          big bathymetry grid
% bath_small        small bathymetry grid (new method, modified 2D taper)
% bath_small_old    small bathymetry grid (old method, 2D cosine taper)
%
% SEE ALSO:
% TAPERBATHYMETRY, SIMULOSL
% 
% Last modified by sirawich-at-princeton.edu, 03/23/2026

defval('s2', 10000)
defval('nu', 0.5)
defval('rho', 2000)
defval('dx', 250)
defval('big', 257)
defval('small', 96)
defstruct('taper', 'trim', 8)
defstruct('taper', 'cosine', 16)

% generate a random field
p.NyNx = [big big];
p.dydx = [dx dx];
p.blurs = Inf;
p.taper = 1;
th0 = [s2 nu rho];
topo_big = v2s(simulosl(th0, p));

% taper the big random field
topo_big_trim = topo_big(1+taper.trim:big-taper.trim, ...
    1+taper.trim:big-taper.trim);
topo_big_trim_taper = taperbathymetry(topo_big_trim, ...
    taper.cosine / (big - 2*taper.trim));
topo_big_taper = ones(size(topo_big)) * topo_big_trim_taper(1,1);
topo_big_taper(1+taper.trim:big-taper.trim, ...
    1+taper.trim:big-taper.trim) = topo_big_trim_taper;

% taper the small random field
halfdiff = (big - small) / 2;
topo_small = topo_big(1+halfdiff:big-halfdiff, 1+halfdiff:big-halfdiff);
topo_small_trim = topo_small(1+taper.trim:small-taper.trim, ...
    1+taper.trim:small-taper.trim);
topo_small_trim_taper = taperbathymetry(topo_small_trim, ...
    taper.cosine / (small - 2*taper.trim));
topo_small_taper = ones(size(topo_small)) * topo_small_trim_taper(1,1);
topo_small_taper(1+taper.trim:small-taper.trim, ...
    1+taper.trim:small-taper.trim) = topo_small_trim_taper;

% taper the small random field in an old fashion
w = shanning2d(size(topo_small_trim), ...
    taper.cosine / (small - 2*taper.trim), 0, 0);
topo_small_trim_taper_old = topo_small_trim .* w;
topo_small_taper_old = ones(size(topo_small)) * ...
    topo_small_trim_taper_old(1,1);
topo_small_taper_old(1+taper.trim:small-taper.trim, ...
    1+taper.trim:small-taper.trim) = topo_small_trim_taper_old;

if plt
    figure
    set(gcf, 'Units', 'inches', 'Position', [0 1 10 5])
    

% collect the outputs
varns = {topo_big_taper, topo_small_taper, topo_small_taper_old};
varargout = varns(1:nargout);
end