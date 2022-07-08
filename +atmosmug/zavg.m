%ZAVG  Conservatively transform from one vertical grid to another
%
%   VVq = ZAVG(PE, VV, PEq) returns a vector VLq defined on the layers of the
%   vertical grid with edge pressures PEq that has the same vertical averages
%   as the vector VV defined on the layers of the vertical grid with edge
%   pressures PPq.
%
%   [VVq, WW] = ZAVG(PE, VV, PEq) returns the weights WW of the
%   interpolation.

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%
% Changelog:
% 2019/12/30	Enabled grids with different surface pressures, etc.
%==============================================================================%
function [avg, wgt] = zavg(pemod, f, peobs)

LDEBUG = 0;

nsig  = numel(pemod) - 1;
nlays = numel(peobs) - 1;

grdpe = interp1(pemod, [1:nsig+1]', peobs, 'linear', 'extrap');
% Restrict obs layers to model range
grdpe = max(min(grdpe,nsig+1),1);

top2bot = (grdpe(1) < grdpe(end));

avg = zeros(nlays,1);
wgt = zeros(nlays,nsig);

for j = 1:nlays
   if (top2bot)
     dz1 = grdpe(j+1);
     dz2 = grdpe(j);
   else
     dz1 = grdpe(j);
     dz2 = grdpe(j+1);
   end

   iz1 = floor(dz1);
   iz2 = floor(dz2);
   if (nsig < iz1), iz1 = nsig; end

%  Obs layer falls within single model layer
   if (iz1 == iz2)
      avg(j) = f(iz1);
      wgt(j,iz1) = 1;
      continue;
   end

%  Obs layer covers multiple model layers
   dlay = 0;
   for l = iz1:-1:iz2
      delz = 1;
      if (l == iz1), delz =         dz1 - iz1;  end
      if (l == iz2), delz = delz - (dz2 - iz2); end

      delp = pemod(l) - pemod(l+1);

      dlay     = dlay   +      delp*delz;
      avg(j)   = avg(j) + f(l)*delp*delz;
      wgt(j,l) =               delp*delz;
   end

   avg(j)   = avg(j)   / dlay;
   wgt(j,:) = wgt(j,:) / dlay;

   if (LDEBUG)
      sprintf('iz1  = %e\t iz2 = %e\n', iz1,  iz2);
      sprintf('dlay = %e\t avg = %e\n', dlay, avg(j));
   end
end
