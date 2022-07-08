%INTERPDPOS  Conservative interpolation of daily fluxes from monthly fluxes
%
%   XD = INTERPDPOS(XM, H) computes the smoothest daily interpolant XD from
%   the monthly values XM such that the monthly means of XD are XM.
%
%   The matrix H is the system of linear equations for the unconstrained
%   problem.  It is required as an input to avoid the cost of reconstructing
%   it many times.

% Author: Brad Weir
%==============================================================================%
function dayvals = interpdpos(movals, QQ, AAeq)

MAXITS  = 32;
TOTMONS = numel(movals);
TOTDAYS = size(QQ,1);

HH = [QQ'*QQ, AAeq'; AAeq, sparse(TOTMONS, TOTMONS)];
bb = [zeros(TOTDAYS,1); movals];
xx = HH \ bb;
dayvals = xx(1:TOTDAYS);

% b. Iterate until constrained solution is found
ioff = [];
noff = numel(ioff);
for it = 1:MAXITS
  nold = noff;
  ioff = union(ioff, find(dayvals < 0));
  noff = numel(ioff);

% Finished when no new negatives are added
  if (noff == nold), break; end

% Compute solution while fixing all negative values at 0
  AAp = sparse(noff, TOTDAYS, noff);
  for np = 1:noff
    AAp(np,ioff(np)) = 1;
  end

  HHp = [QQ'*QQ, [AAeq; AAp]'; [AAeq; AAp], zeros(TOTMONS+noff)];
  bbp = [zeros(TOTDAYS,1); movals; zeros(noff,1)];
  xx  = HHp \ bbp;

  dayvals = xx(1:TOTDAYS);
  dayvals = max(dayvals, 0);
end
