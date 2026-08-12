function [x, y, xxm, yym, xxn, yyn, zzm, zzn] = MESHvsND(M, N, f)
% [x, y, xxm, yym, xxn, yyn, zzm, zzn] = NDVMESH(M, N, f)
%
% Compares/contrasts NDGRID and MESHGRID, because they are so confusing,
% neither provides good explanation of how we should use them.
%
% INPUT:
% M         number of x-coordinates (x = 1:M)
% N         number of y-coordinates (y = 1:N)
% f         function of the elevation z = f(xx,yy)
%           For example, f = @(x,y) (sin(x/10).*cos(y/20))
%
% OUTPUT:
% x         x-cooridnates of points
% y         y-coordinates of points
% xxm       x-coordinates over grid from MESHGRID
% yym       y-coordinates over grid from MESHGRID
% xxn       x-coordinates over grid from NDGRID
% yyn       y-coordinates over grid from NDGRID
% zzm       elevation from f(xxm, zzm)
% zzn       elevation from f(xxn, zzn)
%
% SEE ALSO:
% MESHGRID, NDGRID
% 
% Last modified by spipatprathanporn@ucsd.edu, 08/12/2026

defval('M', 61)
defval('N', 73)
%defval('f', @(x,y) (sin(x/20 * 2*pi) .* exp(-(y-30).^2 / 100)))
defval('f', @(x,y) (sin(x/20 * 2*pi) .* cos(y/40 * 2*pi)))

x = 1:M;
y = 1:N;
[xxm, yym] = meshgrid(x, y);
[xxn, yyn] = ndgrid(x, y);
zzm = f(xxm, yym);
zzn = f(xxn, yyn);

figure
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 6 10])

subplot('Position', subplotposition(4,2,1))
imagesc(xxm)
grid on
cb = colorbar;
set(get(cb, 'Label'), 'String', 'xxm_{ij}')
axis equal
axis tight
xlabel(sprintf('j = 1:%d', M))
ylabel(sprintf('i = 1:%d', N))
title('xxm')
set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

subplot('Position', subplotposition(4,2,3))
imagesc(yym)
grid on
cb = colorbar;
set(get(cb, 'Label'), 'String', 'yym_{ij}')
axis equal
axis tight
xlabel(sprintf('j = 1:%d', M))
ylabel(sprintf('i = 1:%d', N))
title('yym')
set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

subplot('Position', subplotposition(4,2,5))
imagesc(zzm)
grid on
cb = colorbar;
set(get(cb, 'Label'), 'String', 'zzm_{ij}')
axis equal
axis tight
xlabel(sprintf('j = 1:%d', M))
ylabel(sprintf('i = 1:%d', N))
title('zzm')
set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

subplot('Position', subplotposition(4,2,7))
surface(xxm, yym, zzm)
cb = colorbar;
set(get(cb, 'Label'), 'String', 'zzm')
axis equal
axis tight
xlabel(sprintf('x = 1:%d', M))
ylabel(sprintf('y = 1:%d', N))
title('surface(xxm, yym, zzm)')
set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

%
subplot('Position', subplotposition(4,2,2))
imagesc(xxn)
grid on
cb = colorbar;
set(get(cb, 'Label'), 'String', 'xxn_{ij}')
axis equal
axis tight
xlabel('i')
ylabel('j')
xlabel(sprintf('j = 1:%d', N))
ylabel(sprintf('i = 1:%d', M))
title('xxn')
set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

subplot('Position', subplotposition(4,2,4))
imagesc(yyn)
grid on
cb = colorbar;
set(get(cb, 'Label'), 'String', 'yyn_{ij}')
axis equal
axis tight
xlabel('i')
ylabel('j')
xlabel(sprintf('j = 1:%d', N))
ylabel(sprintf('i = 1:%d', M))
title('yyn')
set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

subplot('Position', subplotposition(4,2,6))
imagesc(zzn)
grid on
cb = colorbar;
set(get(cb, 'Label'), 'String', 'zzn_{ij}')
axis equal
axis tight
xlabel('i')
ylabel('j')
xlabel(sprintf('j = 1:%d', N))
ylabel(sprintf('i = 1:%d', M))
title('zzn')
set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

subplot('Position', subplotposition(4,2,8))
surface(xxn, yyn, zzn)
cb = colorbar;
set(get(cb, 'Label'), 'String', 'zzn')
axis equal
axis tight
xlabel(sprintf('x = 1:%d', M))
ylabel(sprintf('y = 1:%d', N))
title("surface(xxn, yyn, zzn)")
set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

ax0 = subplot('Position', [0.08 0.95 0.84 0.01]);
set(ax0, 'Color', 'none')
set(get(ax0, 'XAxis'), 'Visible', 'off')
set(get(ax0, 'YAxis'), 'Visible', 'off')
func_name = func2str(f);
func_name = replace(func_name, '@(x,y)', '');
func_name = func_name(2:end-1);
title(sprintf('x=1:%d, y=1:%d, z = %s', M, N, func_name))
set(gca, 'FontSize', 12)

set(gcf, 'Renderer', 'painters')
figdisp(mfilename, [], [], 2, [], 'epstopdf')
end