function [t, x, z, ph] = fluidsolidarrivals(z_interface, z_station, rho_f, vp_f, rho_s, vp_s, vs_s, theta, t_max)
% [t, x, z, ph] = FLUIDSOLIDARRIVALS(z_interface, z_station, rho_f, ...
%     vp_f, rho_s, vp_s, vs_s, theta, t_max)
%
% Compute arrival times and displacement amplitude of the waves arriving at
% a stations inside a fluid layer on top of a solid half-space with the
% P-wave entering the system from the bottom. The first P-wave arrival is
% zero and its displacement amplitudes are one.
%
% INPUT:
% z_interface   depth of the interface (given the top of fluid layer is 0)
% z_station     depth of the station
% rho_f         density of the fluid
% vp_f          P-wave velocity in the fluid
% rho_s         density of the solid
% vp_s          P-wave velocity in the solid
% vs_s          S-wave velocity in the solid
% theta         incidence angle of the P-wave entering from the bottom
%               (in radians)
% t_max         maximum time for consideration
%
% OUTPUT:
% t             arrival time
% x             displacement amplitude in the horizontal direction
% z             displacement amplitude in the vertical direction
% ph            phase (1 = P, 2 = SV, + = up, - = down)
%
% EXAMPLE:
% % run a demo
% fluidsolidarrivals('demo');
%
% Last modified by sirawich-at-princeton.edu 11/18/2026

%% demos
if ischar(z_interface) && strcmp(z_interface, 'demo')
    z_interface = 4000;
    z_station = 6000;
    rho_f = 1020;
    vp_f = 1500;
    rho_s = 2500;
    vp_s = 3400;
    vs_s = 1963;
    theta = 20 * pi/180;
    t_max = 100;
    
    [t, x, z, ph] = fluidsolidarrivals(z_interface, z_station, rho_f, ...
        vp_f, rho_s, vp_s, vs_s, theta, t_max);
    
    figure(1)
    set(gcf, 'Units', 'inches', 'Position', [0 1 8 4])
    clf
    ax1 = subplot('Position', [0.09 0.55 0.88 0.33]);
    stem(t(abs(ph)==1), x(abs(ph)==1) * sin(theta), 'LineWidth', 1)
    hold on
    stem(t(abs(ph)==2), x(abs(ph)==2) * sin(theta), 'LineWidth', 1)
    xlim([-2 52])
    ylim([-1.1 1.1])
    grid on
    nolabels(ax1, 1)
    ylabel('x-disp')
    legend('P-wave', 'SV-wave')
    title_string = sprintf(['$$\\rho_s = %.f\\ \\textnormal{kg/m}^3, ' ...
        '\\alpha_s = %.f\\ \\textnormal{m/s}, \\beta_s = %.f\\ ' ...
        '\\textnormal{m/s}, \\rho_f = %.f\\ \\textnormal{kg/m}^3, ' ...
        '\\alpha_f = %.f\\ \\textnormal{m/s} $$ \n $$ \\theta = %.2f ' ...
        '^{\\circ}, z_\\textnormal{interface} = %.f \\textnormal{m}, ' ...
        'z_\\textnormal{station} = %.f \\textnormal{m}$$'], ...
        rho_s, vp_s, vs_s, rho_f, vp_f, theta * 180/pi, ...
        z_interface, z_station);
    title(title_string, 'Interpreter', 'latex')
    set(ax1, 'TickDir', 'out', 'FontSize', 12)
    
    ax2 = subplot('Position', [0.09 0.14 0.88 0.33]);
    stem(t(abs(ph)==1), z(abs(ph)==1) * cos(theta), 'LineWidth', 1)
    hold on
    stem(t(abs(ph)==2), z(abs(ph)==2) * cos(theta), 'LineWidth', 1)
    xlim([-2 52])
    ylim([-1.1 1.1])
    grid on
    xlabel('time (s)')
    ylabel('z-disp')
    set(ax2, 'TickDir', 'out', 'FontSize', 12)
    
    set(gcf, 'Renderer', 'painters')
    figdisp(sprintf('%s_demo', mfilename), [], [], 2, [], 'epstopdf');
    return
end

%% angles of incident for each wave
theta_s = asin(vs_s / vp_s * sin(theta));
theta_f = asin(vp_f / vp_s * sin(theta));

%% displacement amplitude coefficients
% solid to fluid
[~, D_S2F] = fluidsolidcoefficients2(rho_f, vp_f, rho_s, vp_s, vs_s, theta);
D_S2F_RP   = D_S2F(1);   % reflected P
D_S2F_RSV  = D_S2F(2);   % reflected SV
D_S2F_T    = D_S2F(3);   % transimitted acoustic

% free-surface reflection in an acoustic media
D_FREE_R   = -1;

% fluid to solid
[~, D_F2S] = fluidsolidcoefficients3(rho_f, vp_f, rho_s, vp_s, vs_s, theta_f);
D_F2S_R    = D_F2S(1);   % reflected acoustic
D_F2S_TP   = D_F2S(2);   % transmitted P
D_F2S_TSV  = D_F2S(3);   % transmitted SV

%% sign to convert down-going wave to displacement
SGN_P_X      =  1;
SGN_P_Z      = -1;
SGN_SV_X     =  1;
SGN_SV_Z     =  1;

%% determine which layer the station is in
if (z_interface > z_station)
    layer = "fluid";
else
    layer = "solid";
end

t = 0;
x = 1;
z = 1;
ph = 1;
if strcmp(layer, "solid")
    % travel time of first down-going P-wave
    tPP = 2 * (z_station - z_interface) * cos(theta) / vp_s;
    
    % travel time of first down-going SV-wave
    tPSV = (z_station - z_interface) * (cos(theta) / vp_s + ...
        cos(theta_s) / vs_s);
    
    % round-trip travel time of acoustic wave in the fluid layer
    tFF = 2 * z_interface * cos(theta_f) / vp_f;
    
    % number of rounds
    nP = floor((t_max - tPP) / tFF);
    nSV = floor((t_max - tPSV) / tFF);
    
    % down-going P-wave
    % [P; P^bathP; P(A^A)P]     % A denote acoustic phase
    t = [t; tPP + (0:nP)' * tFF];
    x = [x; [D_S2F_RP; D_S2F_T * (D_FREE_R .^ (1:nP)') .* ...
        (D_F2S_R .^ (0:nP-1)') * D_F2S_TP] * SGN_P_X];
    z = [z; [D_S2F_RP; D_S2F_T * (D_FREE_R .^ (1:nP)') .* ...
        (D_F2S_R .^ (0:nP-1)') * D_F2S_TP] * SGN_P_Z];
    ph = [ph; -ones(nP+1, 1)];
    
    % down-going SV-wave
    % [P; P^bathS; P(A^A)S]     % A denote acoustic phase
    t = [t; tPSV + (0:nSV)' * tFF];
    x = [x; [D_S2F_RSV; D_S2F_T * (D_FREE_R .^ (1:nSV)') .* ...
        (D_F2S_R .^ (0:nSV-1)') * D_F2S_TSV] * SGN_SV_X * ...
        cos(theta_s) / sin(theta)];
    z = [z; [D_S2F_RSV; D_S2F_T * (D_FREE_R .^ (1:nSV)') .* ...
        (D_F2S_R .^ (0:nSV-1)') * D_F2S_TSV] * SGN_SV_Z * ...
        sin(theta_s) / cos(theta)];
    ph = [ph; -repmat(2, [nSV+1, 1])];
else
    % travel time of the first down-going acoustic wave
    tAA = 2 * z_station * cos(theta_f) / vp_f;
    
    % round-trip travel time of acoustic wave in the fluid layer
    tFF = 2 * z_interface * cos(theta_f) / vp_f;
    
    % number of rounds
    n_up = floor(t_max / tFF);
    n_down = floor((t_max - tAA) / tFF);
    
    % up-going acoustic wave arrivals
    t = [t; (1:n_up)' * tFF];
    x = [x; (D_FREE_R * D_F2S_R) .^ (1:n_up)'];
    z = [z; (D_FREE_R * D_F2S_R) .^ (1:n_up)'];
    ph = [ph; ones(n_up, 1)];
    
    % down-going acoustic wave arrivals
    t = [t; tAA + (0:n_down)' * tFF];
    x = [x; D_FREE_R * (D_F2S_R * D_FREE_R) .^ (0:n_down)' * SGN_P_X];
    z = [z; D_FREE_R * (D_F2S_R * D_FREE_R) .^ (0:n_down)' * SGN_P_Z];
    ph = [ph; -ones(n_down+1, 1)];
end

% sort arrival by time
[t, it] = sort(t);
x = x(it);
z = z(it);
ph = ph(it);

% vertical incident does not have any horizontal displacement
if abs(theta) == 0
    x = 0 * t;
end
end