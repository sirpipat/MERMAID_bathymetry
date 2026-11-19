function varargout = layerarrivals_demo(model)
% [layers, p, t, x, z, ph, phtex] = LAYERARRIVALS_DEMO(model)
%
% Demos of LAYERARRIVALS to demonstrate its ability to identify phases
% arrivals in layered media as well as their arrival times and relative
% amplitudes. It labels arrival times on top of the seismograms from
% FKPROPAGATRIX.
%
% INPUT:
% model         which model used in the demonstration
%     1             Ocean + Crust from Pipatprathanporn & Simons (2024)
%     2             Ocean + Crust + Mantle
%     3             Ocean + Sediment + Crust + Mantle
%
% OUTPUT:
% layers        layered model used in the demonstration
% p             ray parameter
% t             arrival time
% x             displacement amplitude in the horizontal direction
% z             displacement amplitude in the vertical direction
% ph            phase names
% phtex         phase names in tex syntax
%
% SEE ALSO:
% LAYERARRIVALS, FKPROPAGATRIX
%
% Last modified by spipatprathanporn@ucsd.edu, 11/18/2026

defval('model', 1)
if model == 1
    layers = struct('rho', {1020 2500}', ...
        'vp', {1500 3400}', ...
        'vs', {0 1963}', ...
        'ztop', {0 -4129.4}');
    p = sin(10 * pi/180) / 1500;
elseif model == 2
    layers = struct('rho', {1020 2850 3300}', ...
        'vp', {1500 6500 8080}', ...
        'vs', {0 3700 4490}', ...
        'ztop', {0 -4129.4 -12940.2}');
    p = sin(25.45 * pi/180) / 8080;
elseif model == 3
    layers = struct('rho', {1020 1820 2850 3300}', ...
        'vp', {1500 1750 6500 8080}', ...
        'vs', {0 540 3700 4490}', ...
        'ztop', {0 -4129.4 -4329.4 -12940.2}');
    p = sin(25.45 * pi/180) / 8080;
end

% Propagation matrix method
% 100-Hz Gaussian: almost Dirac Comb
% No phase labeled
t = (-100:0.005:100)'; x = exp(-(100*(t)).^2 / 2);
[tout, yout] = fkpropagatrix(layers, p, [-12000 -20000], ...
    [0 -24250], [t x], 'P', 100);

% Reflectivity method
% with phases
[t,x,z,ph,phtex] = layerarrivals(layers,p,[-12000 -20000],[0 -24250], 30, 0.15);

figure
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 8 4])
subplot('Position', [0.08 0.65 0.84 0.2563])
plot(tout, real(yout.ux) / max(abs(yout.ux)), 'LineWidth', 1, 'Color', 'k')
grid on
xlim([0 30])
ylim([-1 1] * 1.3)
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'on')
nolabels(gca, 1);
title('u_{x}')
hold on
scatter(t, x, 'LineWidth', 1, 'MarkerEdgeColor', [0.3 0.7 1])

subplot('Position', [0.08 0.06 0.84 0.50])
plot(tout, real(yout.uz) / max(abs(yout.uz)), 'LineWidth', 1, 'Color', 'k')
grid on
xlim([0 30])
ylim([-1 1] * 1.3)
set(gca, 'TickDir', 'out', 'FontSize', 11, 'Box', 'on')
xlabel('time (s)')
title('u_{z}')
hold on
scatter(t, z, 'LineWidth', 1, 'MarkerEdgeColor', [0.3 0.7 1])
legend('propagator matrix (100 Hz Gaussian)', ...
    'reflectivity method (click a circle for phase name)', ...
    'Location', 'southoutside')
set(gcf, 'Renderer', 'painters')

% add callback function to display phases
dcm = datacursormode(gcf);
set(dcm, 'Enable', 'on');
set(dcm, 'UpdateFcn', @(obj, event_obj) myDataTipCallback(obj, ...
    event_obj, pptex));

% outputs
outputs = {layers, p, t, x, z, ph, phtex};
varargout = outputs(1:nargout);
end

% --- Callback Function ---
function output_txt = myDataTipCallback(~, event_obj, extraInfo)
    % Get the current X and Y coordinates
    pos = get(event_obj, 'Position');
    xVal = pos(1);
    yVal = pos(2);
    
    % Find the index of the clicked point (optional, to pull related data)
    target = get(event_obj, 'Target');
    xData = get(target, 'XData');
    [~, idx] = min(abs(xData - xVal));
    
    % Build the custom text cell array to display in the data tip
    output_txt = {
        sprintf('X: %.2f', xVal), ...
        sprintf('Y: %.2f', yVal), ...
        extraInfo{idx} % or look up extra info mapped to 'idx'
    };
end