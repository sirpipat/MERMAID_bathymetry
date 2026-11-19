function varargout = fkpropagatrix(layers, p, x0, x1, stf, tipe, fmax)
% [tout, y, layers] = FKPROPAGATRIX(layers, p, [x0 z0], [x1 z1], [t u], ...
%     tipe, fmax)
%
% Analytically solve wave equation in homogenous flat layered media using
% propagation matrix.
%
% INPUT:
% layers        layered media from top to bottom: struct with these fields
% - rho             density
% - vp              P-wave / acoustic wave speed
% - vs              S-wave speed (0 for acoustic wave)
% - ztop            elevation at the top of the layer
% p             ray parameter (rad s)
% [x0 z0]       wavefront origin
% [x1 z1]       point of interest
% [t u]         source-time function of the input displacement: u(t)
% tipe          'P' or 'SV'
% fmax          maximum frequency used for calculation
%
% OUTPUT:
% tout          output time
% y             wavefield solution
% - ux              horizontal displacement
% - uz              vertical displacement
% - Txz             shear traction
% - Tzz             vertical, normal traction
% - p               hydrostatic stress
% layers        layers struct with additional information
% 
% Last modified by spipatprathanporn@ucsd.edu, 11/18/2026

%% Part I: compute incident angle in each layer
NUMLAYERS = length(layers);
for ii = 1:NUMLAYERS
    layers(ii).theta_p = asin(p * layers(ii).vp);
    layers(ii).theta_s = asin(p * layers(ii).vs);
    layers(ii).p = p;
    layers(ii).eta_p = sqrt(1/layers(ii).vp^2 - p^2);
    layers(ii).issolid = (layers(ii).vs > 0);
    if layers(ii).issolid
        layers(ii).eta_s = sqrt(1/layers(ii).vs^2 - p^2);
    else
        layers(ii).eta_s = nan;
    end
end

%% Part II: decompose the source-time function
t = stf(:,1);
u = stf(:,2);
N = length(t);
fs = (N-1) / (t(end)-t(1));
U = fftshift(fft(u));
F = fs/N * ((0:N-1) - ceil((N-1)/2));       % frequency centered at 0
wh = (abs(F) <= fmax);
U = U(wh);
F = F(wh);
omega = (2*pi) * F;
NUMOMEGA = length(omega);

% displacement and traction at the target point
Y = nan(NUMOMEGA, 5);

%% Part III: compute propagation matrices and propagate
% determine which layer the target point is in
assert(x1(2) < layers(1).ztop, 'Error: z1 is above the free surface')
targetlayer = 1;
while targetlayer < NUMLAYERS && x1(2) < layers(targetlayer+1).ztop
    targetlayer = targetlayer + 1;
end

% loop over angular frequency
for ii = 1:NUMOMEGA
    %% loop over layers to construct G and Q matrices
    for jj = 1:NUMLAYERS
        if jj == NUMLAYERS
            layers(jj).G = layerpropagator(layers(jj), omega(ii), -x0(1), ...
                layers(jj).ztop - x0(2));
        else
            layers(jj).G = layerpropagator(layers(jj), omega(ii), 0, ...
                layers(jj).ztop - layers(jj+1).ztop);
        end
        if layers(jj).issolid
            layers(jj).G_up = layers(jj).G(3:4,3:4);
        else
            layers(jj).G_up = layers(jj).G(2,2);
        end

        % Q for boundary at the top of each layer
        if jj == 1
            [Q,Q_left,Q_right] = potentialtowavefield(layers(jj), ...
                omega(ii), ~layers(jj).issolid, true);
        else
            [Q,Q_left,Q_right] = potentialtowavefield(layers(jj), ...
                omega(ii), ~layers(jj).issolid || ~layers(jj-1).issolid, ...
                false);
        end
        layers(jj).Q_top = Q;
        layers(jj).Q_top_left = Q_left;
        layers(jj).Q_top_right = Q_right;

        % Q for boundary at the bottom of each layer
        if jj == NUMLAYERS
            [Q,Q_left,Q_right] = potentialtowavefield(layers(jj), ...
                omega(ii), ~layers(jj).issolid, false);
        else
            [Q,Q_left,Q_right] = potentialtowavefield(layers(jj), ...
                omega(ii), ~layers(jj).issolid || ~layers(jj+1).issolid, ...
                false);
        end
        layers(jj).Q_bot = Q;
        layers(jj).Q_bot_left = Q_left;
        layers(jj).Q_bot_right = Q_right;
    end
    
    %% iteratively solve for directional transmission/reflection matrices
    layers(1).P_below = -layers(1).Q_top_left \ layers(1).Q_top_right;
    layers(1).P_above = zeros(size(layers(1).P_below));
    for jj = 2:NUMLAYERS
        P = [layers(jj-1).Q_bot_right + layers(jj-1).Q_bot_left * ...
            (layers(jj-1).G_up * layers(jj-1).P_below * ...
            layers(jj-1).G_up), -layers(jj).Q_top_left] \ ...
            layers(jj).Q_top_right;
        if layers(jj).issolid
            layers(jj).P_below = P(end-1:end,:);
        else
            layers(jj).P_below = P(end,:);
        end
        if layers(jj-1).issolid
            layers(jj).P_above = P(1:2,:);
        else
            layers(jj).P_above = P(1,:);
        end
    end

    %% propagate input signal to the location of interest
    % initialize input signal
    if strcmpi(tipe, 'P')
        varphi_up = [sin(layers(NUMLAYERS).theta_p); ...
            cos(layers(NUMLAYERS).theta_p)] * U(ii);
    else
        varphi_up = [cos(layers(NUMLAYERS).theta_p); ...
            sin(layers(NUMLAYERS).theta_p)] * U(ii);
    end
    Q = potentialtowavefield(layers(NUMLAYERS), omega(ii), ...
        false, false);
    varphi_up = Q(1:2,3:4) \ varphi_up;
    % varphi_up = S(ii);

    % propagate to the top of the target layer
    for jj = NUMLAYERS:-1:targetlayer+1
        varphi_up = layers(jj).P_above * layers(jj).G_up * varphi_up;
    end
    varphi_up = layers(targetlayer).G_up * varphi_up;
    % construct downgoing waves at the top of the target layer
    varphi_down = layers(targetlayer).P_below * varphi_up;
    varphi = [varphi_down; varphi_up];
    % propagate to the point of interest
    varphi = layerpropagator(layers(targetlayer), omega(ii), x1(1), ...
        x1(2)-layers(targetlayer).ztop) * varphi;

    % compute the displacement and traction
    Y(ii,1:4) = (potentialtowavefield(layers(targetlayer), omega(ii), ...
        false, false) * varphi).'; % .' is transpose, not adjoint.
    % compute the pressure
    Y(ii,5) = potentialtopressure(layers(targetlayer), omega(ii)) * varphi;
end

%% Part IV: outputs
% output time samples
tout = linspace(0, t(end) - t(1), NUMOMEGA)';

% convert the output to time domain
y.ux  = ifftshift(ifft(ifftshift(Y(:,1))));
y.uz  = ifftshift(ifft(ifftshift(Y(:,2))));
y.Txz = ifftshift(ifft(ifftshift(Y(:,3))));
y.Tzz = ifftshift(ifft(ifftshift(Y(:,4))));
y.p   = ifftshift(ifft(ifftshift(Y(:,5))));

% collect the output
outputs = {tout, y, layers};
varargout = outputs(1:nargout);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G = LAYERPROPAGATOR(layer, omega, dx, dz)
%
% Computes the phase shift for downgoing-P, downgoing-SV, upgoing-P, and
% upgoing-SV given the position shift within a layer and the wave's angular
% frequency.
%
% INPUT:
% layer         struct with these fields
% - rho             density
% - vp              P-wave / acoustic wave speed
% - vs              S-wave speed (0 for acoustic wave)
% - ztop            elevation at the top of the layer
% omega         angular frequency
% dx            horizontal displacement
% dz            vertical displacement
%
% OUTPUT:
% G             diagonal matrix of the phase shifts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = layerpropagator(layer, omega, dx, dz)
if layer.issolid
    G = zeros(4,4);
    G(1,1) = exp(-1i*omega*(layer.p*dx - layer.eta_p*dz));
    G(2,2) = exp(-1i*omega*(layer.p*dx - layer.eta_s*dz));
    G(3,3) = exp(-1i*omega*(layer.p*dx + layer.eta_p*dz));
    G(4,4) = exp(-1i*omega*(layer.p*dx + layer.eta_s*dz));
else
    G = zeros(2,2);
    G(1,1) = exp(-1i*omega*(layer.p*dx - layer.eta_p*dz));
    G(2,2) = exp(-1i*omega*(layer.p*dx + layer.eta_p*dz));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Q,Q_left,Q_right] = POTENTIALTOWAVEFIELD(layer, omega, ...
%     isfluidboundary, isfree)
%
% Computes the matrix that convert wave potential to wavefield
% (displacements and tractions). 
%
% INPUT:
% layer             struct with these fields
% - rho                 density
% - vp                  P-wave / acoustic wave speed
% - vs                  S-wave speed                (solid medium only)
% - p                   ray parameter
% - eta_p               P-wave vertical slowness
% - eta_s               S-wave vertical slowness    (solid medium only)
% - issolid             whether the layer is solid
% omega             angular frequency
% isfluidboundary   whether it is for a fluid boundary
% isfree            whether it is for a free-surface boundary
%
% OUTPUT:
% Q                 the matrix that convert wave potential to wavefield
% Q_left            left half of Q matrix
% Q_right           right half of Q matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q, Q_left, Q_right] = potentialtowavefield(layer, omega, ...
        isfluidboundary, isfree)
defval('isfluidboundary', false)
defval('isfree', false)
if omega == 0
    omega = 1;
end
if layer.issolid
    Q = zeros(4,4);

    % for calculating horiztonal displacement
    Q(1,1) = -1i * omega * layer.p;
    Q(1,2) = -1i * omega * layer.eta_s;
    Q(1,3) = -1i * omega * layer.p;
    Q(1,4) = +1i * omega * layer.eta_s;

    % for calculating vertical displacement
    Q(2,1) = +1i * omega * layer.eta_p;
    Q(2,2) = -1i * omega * layer.p;
    Q(2,3) = -1i * omega * layer.eta_p;
    Q(2,4) = -1i * omega * layer.p;
    
    % for calculating shear traction to a horizontal plane
    Q(3,1) =  2 * omega^2 * layer.rho * layer.vs^2 * layer.p * layer.eta_p;
    Q(3,2) =      omega^2 * layer.rho * (1 - 2 * layer.vs^2 * layer.p^2);
    Q(3,3) = -2 * omega^2 * layer.rho * layer.vs^2 * layer.p * layer.eta_p;
    Q(3,4) =      omega^2 * layer.rho * (1 - 2 * layer.vs^2 * layer.p^2);

    % for calculating normal traction to a horizontal plane
    Q(4,1) = -    omega^2 * layer.rho * (1 - 2 * layer.vs^2 * layer.p^2);
    Q(4,2) =  2 * omega^2 * layer.rho * layer.vs^2 * layer.p * layer.eta_s;
    Q(4,3) = -    omega^2 * layer.rho * (1 - 2 * layer.vs^2 * layer.p^2);
    Q(4,4) = -2 * omega^2 * layer.rho * layer.vs^2 * layer.p * layer.eta_p;
else
    Q = zeros(4,2);

    % for calculating horiztonal displacement
    Q(1,1) = -1i * omega * layer.p;
    Q(1,2) = -1i * omega * layer.p;

    % for calculating vertical displacement
    Q(2,1) = +1i * omega * layer.eta_p;
    Q(2,2) = -1i * omega * layer.eta_p;

    % for calculating shear traction to a horizontal plane
    % third row is always zero

    % for calculating normal traction to a horizontal plane
    Q(4,1) = -omega^2 * layer.rho;
    Q(4,2) = -omega^2 * layer.rho;
end

% some components of displacement are not considered in the boundary
% conditions for free-surface and fluid-solid boundaries
if isfree
    Q = Q(3:4,:);
elseif isfluidboundary
    Q = Q(2:4,:);
end

if size(Q,2) == 4
    Q_left  = Q(:,1:2);
    Q_right = Q(:,3:4);
else
    Q_left  = Q(:,1);
    Q_right = Q(:,2);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p = POTENTIALTOPRESSURE(layer, omega)
%
% Computes the matrix that convert the wave potential to pressure.
%
% INPUT:
% layer             struct with these fields
% - rho                 density
% - vp                  P-wave / acoustic wave speed
% - vs                  S-wave speed
% - issolid             whether the layer is solid
% omega             angular frequency
%
% OUTPUT:
% p                 the matrix that convert the wave potential to pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = potentialtopressure(layer, omega)
% (\lambda + 2/3 \mu) * \nabla^2 \phi
% rho*(alpha^2-2*beta^2) + 2/3*rho*beta^2
% rho*(alpha^2-4/3*beta^2)
% Bulk modulus
kappa = layer.rho * (layer.vp^2 - 4/3 * layer.vs^2);
if layer.issolid
    p = zeros(1,4);
    p(1) = -omega^2 * kappa / layer.vp^2;
    p(2) = -omega^2 * kappa / layer.vs^2;
    p(3) = p(1);
    p(4) = p(2);
else
    p = -omega^2 * kappa / layer.vp^2 * [1 1];
end
end