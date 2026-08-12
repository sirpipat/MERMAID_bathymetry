function varargout = figuretitle(fig, titlestr, substr, varargin)
% [ax, t, s] = FIGURETITLE(fig, titlestr, substr, Name, Value)
%
% Adds the tile of the figure as a separate axes. Please make sure you left
% some space near the top of the figure. The position of the axes is always
% [0.08 0.93 0.84 0.01].
%
% INPUTS:
% fig           target figure       [Default: gcf]
% titlestr      title string        [Default: 'title']
% substr        subtitle string     [Default: '']
% Name, Value   Name-Value arguments for TITLE
%
% OUTPUTS:
% ax            axes handle
% t             title object handle
% s             subtitle object handle
%
% SEE ALSO:
% TITLE
%
% Last modified by spipatprathanporn@ucsd.edu, 08/12/2026

defval('fig', gcf)
defval('titlestr', 'title')
defval('substr', '')

% sets the target figure
figure(fig);

% adds the title axes
ax = subplot('Position', [0.08 0.93 0.84 0.01]);
[t,s] = title(titlestr, substr, varargin{:});
nolabels(ax, 3);
set(get(ax, 'XAxis'), 'Visible', 'off')
set(get(ax, 'YAxis'), 'Visible', 'off')
set(ax, 'Color', 'none')

% collects outputs
outputs = {ax, t, s};
varargout = outputs(1:nargout);
end