function w = shanning2d(n, r, shape, sac)
% w = SHANNING2D([ny nx], r, shape, sac)
%
% Calculates 2D Hanning windows of a certain length and a certain width
% fraction the exact way SAC says it is doing it...
%
% INPUT:
%
% [ny nx]       The required length of the window
% r             The fraction of the window that is tapered [default 0.5]
% shape         shape of the window
%               0 Use rectangular window [default]
%               1 Use circlular window, diameter is the shorter dimension
% sac           0 Use MATLAB [default]
%               1 Use actual SAC
%
% OUTPUT:
%
% w       The Hanning window 
%
% SEE ALSO:
% SHANNING
%
% Last modified by sirawich-at-princeton.edu, 10/06/2026

defval('r',0.5)
defval('shape', 0)
defval('sac',0)

if length(n) == 1
    nx = n;
    ny = n;
else
    nx = n(2);
    ny = n(1);
end


wx = shanning(nx, r, sac);
wy = shanning(ny, r, sac);

if shape == 0
    w = wx' .* wy;
else
    % construct the radius on a grid
    x = linspace(-nx+1, nx-1, nx)';
    y = linspace(-ny+1, ny-1, ny)';
    [xx, yy] = meshgrid(x, y);
    rr = sqrt(xx.^2 + yy.^2);
    
    % construct the radius with known Hanning window values
    if nx < ny
        wr = wx;
        wr = [wr; 0];
        r0 = [x; 2*ny];
    else
        wr = wy;
        wr = [wr; 0];
        r0 = [y; 2*nx];
    end
    
    % interpolate
    w = interp1(r0, wr, rr, 'Linear');
end
    