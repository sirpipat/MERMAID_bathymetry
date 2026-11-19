function varargout = layerarrivals(layers, p, x0, x1, tmax, ampmin)
% [t, x, z, ph, phtex] = LAYERARRIVALS(layers, p, x0, x1, tmax)
%
% Compute arrival times and displacement amplitude of the waves arriving at
% a stations inside a layer with a wave entering the system from the 
% bottom. The first P-wave arrival is zero and its displacement amplitudes
% are one.
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
% tmax          cutoff time
% ampmin        cutoff amplitude relative to the largest arrival
%
% OUTPUT:
% t             arrival time
% x             displacement amplitude in the horizontal direction
% z             displacement amplitude in the vertical direction
% ph            phase names
% phtex         phase names in tex syntax
%
% SEE ALSO:
% FLUIDSOLIDARRIVALS, FKPROPAGATRIX
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

%% Part II: compute the amplitude ratio
for ii = 1:NUMLAYERS
   if ii == 1
       D = amplituderatio([], layers(ii));
       layers(ii).D_top = D;
   else
       D = amplituderatio(layers(ii-1), layers(ii));
       layers(ii).D_top = D;
       layers(ii-1).D_bot = D;
   end
end

%% Part III: compute travel time from one point to another
% determine which layer the target point is in
assert(x1(2) < layers(1).ztop, 'Error: z1 is above the free surface')
targetlayer = 1;
while targetlayer < NUMLAYERS && x1(2) < layers(targetlayer+1).ztop
    targetlayer = targetlayer + 1;
end

% propagation time from top to bottom and vice-versa
for ii = 1:NUMLAYERS-1
    layers(ii).tP = (layers(ii).ztop - layers(ii+1).ztop) * layers(ii).eta_p;
    layers(ii).tS = (layers(ii).ztop - layers(ii+1).ztop) * layers(ii).eta_s;
end
% propagation from the origin point to the top of half space
layers(NUMLAYERS).tP = -x0(1) * p + (layers(NUMLAYERS).ztop - x0(2)) * ...
    layers(NUMLAYERS).eta_p;
layers(NUMLAYERS).tS = -x0(1) * p + (layers(NUMLAYERS).ztop - x0(2)) * ...
    layers(NUMLAYERS).eta_s;

% propagation from the interfaces right above and/or below the target
if targetlayer < NUMLAYERS
    tP_up = x1(1) * p + (x1(2) - layers(targetlayer+1).ztop) * ...
        layers(targetlayer).eta_p;
    tS_up = x1(1) * p + (x1(2) - layers(targetlayer+1).ztop) * ...
        layers(targetlayer).eta_s;
else
    tP_up = nan;
    tS_up = nan;
end
tP_down = x1(1) * p + (layers(targetlayer).ztop - x1(2)) * ...
    layers(targetlayer).eta_p;
tS_down = x1(1) * p + (layers(targetlayer).ztop - x1(2)) * ...
    layers(targetlayer).eta_s;

%% Part IV: Propagate
paths = [];
h = MinHeap(2);
rp = RayPath(sprintf('%dupP', NUMLAYERS), layers(NUMLAYERS).tP, 0, 1);
rp = rp.setexptime(layers, targetlayer, tP_up, tP_down);
h.push([rp.exptime abs(rp.amp)], rp);

aP_direct = 1;
for ii = NUMLAYERS:-1:targetlayer+1
    phasein = sprintf('%dupP', ii);
    phaseout = sprintf('%dupP', ii-1);
    aP_direct = ampcoeff(layers(ii), phasein, phaseout) * aP_direct;
end
ampmincheck = 0.1 * abs(ampmin * aP_direct);

% loop
MAXNUMPATHS = 10000;
while length(paths) < MAXNUMPATHS && ~h.isEmpty()
    [key, rp] = h.pop();
    n = rp.currlayer();
    dire = rp.currdirection();
    if n == NUMLAYERS && strcmp(dire, 'down')
        continue
    end
    % add reflection phase into phaselist
    if layers(n).issolid
        if strcmp(dire, 'up')
            phaseouts = {sprintf('%ddownP', n), sprintf('%ddownS', n)};
            if n == targetlayer
                phaseouts(end+1:end+2) = {sprintf('%ddownP-end', n), ...
                    sprintf('%ddownS-end', n)};
            end
        else
            phaseouts = {sprintf('%dupP', n), sprintf('%dupS', n)};
            if n == targetlayer
                phaseouts(end+1:end+2) = {sprintf('%dupP-end', n), ...
                    sprintf('%dupS-end', n)};
            end
        end
    else
        if strcmp(dire, 'up')
            phaseouts = {sprintf('%ddownP', n)};
            if n == targetlayer
                phaseouts(end+1) = {sprintf('%ddownP-end', n)};
            end
        else
            phaseouts = {sprintf('%dupP', n)};
            if n == targetlayer
                phaseouts(end+1) = {sprintf('%dupP-end', n)};
            end
        end
        
    end
    % add transmission phase into phaselist
    if strcmp(dire, 'up')
        if n > 1
            if layers(n-1).issolid
                phaseouts(end+1:end+2) = {sprintf('%dupP', n-1), ...
                    sprintf('%dupS', n-1)};
            else
                phaseouts(end+1) = {sprintf('%dupP', n-1)};
            end
            if n-1 == targetlayer
                phaseouts(end+1) = {sprintf('%dupP-end', n-1)};
                if layers(n-1).issolid
                    phaseouts(end+1) = {sprintf('%dupS-end', n-1)};
                end
            end
        end
    else
        if n <= NUMLAYERS-1
            if layers(n+1).issolid
                phaseouts(end+1:end+2) = {sprintf('%ddownP', n+1), ...
                    sprintf('%ddownS', n+1)};
            else
                phaseouts(end+1) = {sprintf('%ddownP', n+1)};
            end
            if n+1 == targetlayer
                phaseouts(end+1) = {sprintf('%ddownP-end', n+1)};
                if layers(n+1).issolid
                    phaseouts(end+1) = {sprintf('%ddownS-end', n+1)};
                end
            end
        end
    end
    
    % propagate!
    for ii = 1:length(phaseouts)
        ampratio = ampcoeff(layers(n), rp.currphase(), phaseouts{ii});
        n_out = sscanf(phaseouts{ii}, '%d', 1);
        if contains(phaseouts{ii}, 'end')
            if contains(phaseouts{ii}, 'P')
                if contains(phaseouts{ii}, 'up')
                    t = tP_up;
                else
                    t = tP_down;
                end
            else
                if contains(phaseouts{ii}, 'up')
                    t = tS_up;
                else
                    t = tS_down;
                end
            end
        else
            if contains(phaseouts{ii}, 'P')
                t = layers(n_out).tP;
            else
                t = layers(n_out).tS;
            end
        end
        % make a new RayPath object
        rp_new = rp.addphase(phaseouts{ii}, t, ampratio);
        % add to the (done) path lists if it is done
        if rp_new.isdone()
            paths = [paths; rp_new];
        elseif rp_new.isokay(tmax, ampmincheck)
            % otherwise, calculate the expected arrival time and then
            % push into the heap
            rp_new = rp_new.setexptime(layers, targetlayer, tP_up, ...
                tP_down);
            h.push(rp_new.key(), rp_new);
        end
    end
end

%% Part V: convert paths to outputs
[~,isortt] = sort([paths.time]);
paths = paths(isortt);
t = paths(1).time;
x = paths(1).getx(layers);
z = paths(1).getz(layers);
ph = {paths(1).phasename()};
phtex = {paths(1).phasename(true)};
for ii = 2:length(paths)
    % same path: sum together, concatenate phase names
    if abs(paths(ii).time - t(end)) < 1e-5
        x(end) = x(end) + paths(ii).getx(layers);
        z(end) = z(end) + paths(ii).getz(layers);
        ph{end} = [ph{end} ' + ' paths(ii).phasename()];
        phtex{end} = [phtex{end} ' + ' paths(ii).phasename(true)];
    % new path: add new entry
    else
        t(end+1) = paths(ii).time;
        x(end+1) = paths(ii).getx(layers);
        z(end+1) = paths(ii).getz(layers);
        ph{end+1} = paths(ii).phasename();
        phtex{end+1} = paths(ii).phasename(true);
    end
end
% normalize
x = x / max(abs(x));
z = z / max(abs(z));
wh = sqrt(x.^2 + z.^2) >= ampmin;
t = t(wh);
x = x(wh);
z = z(wh);
ph = ph(wh);
phtex = phtex(wh);

outputs = {t',x',z',ph',phtex'};
varargout = outputs(1:nargout);
end

%% Helper functions

function [Q,Q_left,Q_right] = potentialtowavefield(layer, omega, ...
        isfluidboundary, isfree)
defval('isfluidboundary', false)
defval('isfree', false)
if omega == 0
    omega = 1;
end
if layer.issolid
    Q = zeros(4,4);

    Q(1,1) = -1i * omega * layer.p;
    Q(1,2) = -1i * omega * layer.eta_s;
    Q(1,3) = -1i * omega * layer.p;
    Q(1,4) = +1i * omega * layer.eta_s;

    Q(2,1) = +1i * omega * layer.eta_p;
    Q(2,2) = -1i * omega * layer.p;
    Q(2,3) = -1i * omega * layer.eta_p;
    Q(2,4) = -1i * omega * layer.p;
    
    Q(3,1) =  2 * omega^2 * layer.rho * layer.vs^2 * layer.p * layer.eta_p;
    Q(3,2) =      omega^2 * layer.rho * (1 - 2 * layer.vs^2 * layer.p^2);
    Q(3,3) = -2 * omega^2 * layer.rho * layer.vs^2 * layer.p * layer.eta_p;
    Q(3,4) =      omega^2 * layer.rho * (1 - 2 * layer.vs^2 * layer.p^2);

    Q(4,1) = -    omega^2 * layer.rho * (1 - 2 * layer.vs^2 * layer.p^2);
    Q(4,2) =  2 * omega^2 * layer.rho * layer.vs^2 * layer.p * layer.eta_s;
    Q(4,3) = -    omega^2 * layer.rho * (1 - 2 * layer.vs^2 * layer.p^2);
    Q(4,4) = -2 * omega^2 * layer.rho * layer.vs^2 * layer.p * layer.eta_s;
else
    Q = zeros(4,2);

    Q(1,1) = -1i * omega * layer.p;
    Q(1,2) = -1i * omega * layer.p;

    Q(2,1) = +1i * omega * layer.eta_p;
    Q(2,2) = -1i * omega * layer.eta_p;

    % third row is always zero

    Q(4,1) = -omega^2 * layer.rho;
    Q(4,2) = -omega^2 * layer.rho;
end

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

function P = potentialratio(layer_above, layer_below)
defval('layer_below', []);

if isempty(layer_above)
    Q_left_above = [];
    Q_right_above = [];
else
    [~,Q_left_above, Q_right_above] = potentialtowavefield(layer_above, ...
        1, ~layer_above.issolid && (~isempty(layer_below) && ...
        ~layer_below.issolid), isempty(layer_below));
end
if isempty(layer_below)
    Q_left_below = [];
    Q_right_below = [];
else
    [~, Q_left_below, Q_right_below] = potentialtowavefield(layer_below, ...
        1, ~layer_below.issolid && (~isempty(layer_above) && ...
        ~layer_above.issolid), isempty(layer_above));
end
% calculate P
P = [Q_right_above, -Q_left_below] \ [-Q_left_above, Q_right_below];
P = real(P);
end

function D = amplituderatio(layer_above, layer_below)
defval('layer_below', []);
P = potentialratio(layer_above, layer_below);

C = zeros(4,4);
if ~isempty(layer_above)
    C(1,1) = layer_above.vp;
    C(2,2) = layer_above.vs;
end
if ~isempty(layer_below)
    C(3,3) = layer_below.vp;
    C(4,4) = layer_below.vs;
end
wh = any(C);
C = C(wh, wh);
D = real(C \ P * C);
end

% phasein/phaseout options
% 'upP'
% 'downP'
% 'upSV'
% 'downSV'
function val = ampcoeff(layer, phasein, phaseout)
% parse the inputs
if layer.issolid
    isolid = 1;

else
    isolid = 0;
end

% Incoming wave to the top interface
if contains(phasein, 'up')
    D = layer.D_top;
    if contains(phasein, 'P')
        jj = size(D,2) - isolid;
    else
        assert(isolid==1, 'No S-wave in fluid layer.')
        jj = size(D,2);
    end

    % Outgoing wave is going up -> Transmission
    if contains(phaseout, 'up')
        if contains(phaseout, 'P')
            assert(size(D,1)-isolid >= 2, 'No P-wave in vacuum.')
            ii = 1;
        else
            assert(size(D,1)-isolid == 3, 'No S-wave in vacuum/fluid.')
            ii = 2;
        end
    % Outgoing wave is going down -> Reflection
    else
        if contains(phaseout, 'P')
            ii = size(D,1) - isolid;
        else
            assert(isolid==1, 'No S-wave in fluid layer.')
            ii = size(D,1);
        end
    end
% Incoming wave to the bottom interface
else
    D = layer.D_bot;
    if contains(phasein, 'P')
        jj = 1;
    else
        assert(isolid==1, 'No S-wave in fluid layer.')
        jj = 2;
    end

    % Outgoing wave is going up -> Reflection
    if contains(phaseout, 'up')
        if contains(phaseout, 'P')
            ii = 1;
        else
            assert(isolid==1, 'No S-wave in fluid layer.')
            ii = 2;
        end
    % Outgoing wave is going down -> Transmission
    else
        if contains(phaseout, 'P')
            assert(size(D,1)-isolid >= 2, 'No P-wave in vacuum.')
            ii = 2 + isolid;
        else
            assert(size(D,1)-isolid == 3, 'No S-wave in vacuum/fluid.')
            ii = 3 + isolid;
        end
    end
end



val = D(ii,jj);
end