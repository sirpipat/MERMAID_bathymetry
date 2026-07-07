function condition = checkslope(itfs, maxslope)
% condition = CHECKSLOPE(itfs)
% condition = CHECKSLOPE(fname)
% condition = CHECKSLOPE(___, maxslope)
%
% Checks whether the maximum of absolute slope of the elevation grid along 
% x and y elevation is always below or equal to a certain value.
%
% INPUT:
% itfs          interfaces, a cell array of struct with following fields
%       SPACING_XI      spacing in x-direction
%       SPACING_ETA     spacing in y-direction
%       Z               elevation grid at (X,Y) or (LON, LAT)
% fname         name of the interface file
% maxslope      maximum absolute slope allowed [default: 1+sqrt(2)]
%
% OUTPUT:
% condition     whether the elevation grid satisfy the criteria
%
% SEE ALSO:
% LOADINTERFACEFILES3D
%
% Last modified by spipatprathanporn@ucsd.edu, 07/06/2026

% maximum slope: equivalent to 75 degrees
MAX_SLOPE = 1 + sqrt(2);

defval('maxslope', MAX_SLOPE)

% if the input is a directory, read the interfacefile
if isstring(itfs) || ischar(itfs)
    itfs = loadinterfacefiles3d(fname);
end

condition = false(size(itfs));
% compute the absolute slope in X and Y direction and check the condition
for ii = 1:length(itfs)
    slope_XI = abs(diff(itfs{ii}.Z, 1, 1)) / itfs{ii}.SPACING_XI;
    slope_ETA = abs(diff(itfs{ii}.Z, 1, 2)) / itfs{ii}.SPACING_ETA;
    condition(ii) = all(slope_XI <= maxslope, 'all') && ...
        all(slope_ETA <= maxslope, 'all');
end
end