function arro = addarrow(x, y, u, v, l, h, a, z, varargin)
% arro = ADDARROW(x, y, u, v, l, h, a, z, Name, Value)
%
% Adds an arrow to a plot with customized arrow length and arrowhead size.
%
% INPUTS:
% x,y           coordinates of the arrow tail
% u,v           direction of the arrow
% l             arrow length
% h             arrowhead size
% a             angle (degrees) of the arrowhead to the arrow
% z             elevation (only for 3-D plot. Leave it blank for 2-D plot)
% Name,Value    Name-Value arguments as in PLOT
%
% OUTPUT:
% arro          arrow object handle
%
% SEE ALSO:
% QUIVER, ARROW
%
% Last modified by spipatprathanporn@ucsd.edu, 08/12/2026

% unit vector of the arrow
V = [u; v] / norm([u; v]);

% unit vector of the arrowhead
a = deg2rad(a);
W1 = [-cos(a) -sin(a); sin(a) -cos(a)] * V;
W2 = [-cos(a) sin(a); -sin(a) -cos(a)] * V;

% points of the arrow
X1 = [x; y];
X2 = X1 + l * V;
X3 = X2 + h * W1;
X4 = X2 + h * W2;
XN = [NaN; NaN];
X = [X1 X2 X3 XN X2 X4];

% plot the arrow
if isempty(z)
    arro = plot(X(1,:), X(2,:), varargin{:});
else
    arro = plot3(X(1,:), X(2,:), z*ones(1,size(X,2)), varargin{:});
end
end