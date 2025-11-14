%PCOLORMX  Extended version of mapping toolbox function pcolorm
%
%   The function pcolorm in the mapping toolbox has a couple of bugs and
%   undesireable features:
%      1) It treats lats and lons as if they denote the bottom-left corner
%         of the gridbox instead of the center
%      2) It deletes the last row and column
%   This function fixes those problems.

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%
% Todo:
% * Verify this works for equal-area grids
%==============================================================================%
function hp = pcolormx(lat, lon, vals, alphas)

dlat = lat(2) - lat(1);
dlon = lon(2) - lon(1);
xlat = [lat - dlat/2; lat(end) + dlat/2];
xlon = [lon - dlon/2; lon(end) + dlon/2];

% Adjust for pole-centered coordinates
xlat = max(xlat, -90);
xlat = min(xlat,  90);

xvals = vals;
xvals = [xvals; xvals(1,:)];
xvals = [xvals, xvals(:,1)];

if (nargin == 3)
  hp = pcolorm(xlat, xlon, xvals);
elseif (nargin == 4)

  avals = alphas;
  if (size(alphas) == size(vals))
    avals = [avals; avals(1,:)];
    avals = [avals, avals(:,1)];
  end

  hp = pcolorm(xlat, xlon, xvals, 'alphadata', avals, ...
              'facealpha', 'flat');
end
