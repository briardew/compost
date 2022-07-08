%GETDP  Build pressure levels from surface pressure
%
%   DP = GETDP(PS,NLEV) returns pressure increments DP.
%
%   [DP,PE] = GETDP(PS,NLEV) returns the pressure edges PE as well.
%
%   [DP,PE,PL] = GETDP(PS,NLEV) returns the pressure level mid-points PL as
%   well.

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%==============================================================================%
function [dp, pe, pl] = getdp(ps, nlev)

geoslevs = load('geoslevs');
ak = geoslevs.(['ak',num2str(nlev)]);
bk = geoslevs.(['bk',num2str(nlev)]);

nlat = size(ps, 2);
nlon = size(ps, 1);

pe = zeros(nlon, nlat, nlev+1);
for kk = 1:nlev+1
% Assumes ps is TOTAL AIR and in Pa
  pe(:,:,kk) = ak(kk) + bk(kk)*ps;
end
% GEOS ordering is top to bottom
pe = flip(pe, 3);
pl = 0.5*(pe(:,:,1:end-1) + pe(:,:,2:end));

dp = diff(pe, [], 3);
