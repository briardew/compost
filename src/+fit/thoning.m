%FIT.THONING  Compute Thoning fit to obs - model residuals
%
%   I don't know

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%
% Changelog:
% 2019/04/16	New version based off of something I gave to Abhishek
%
% TODO:
%==============================================================================%
function bias = thoning(dnobs, omfin, isok)

JDAYS  = 365.25;
NPOLY  = 2;
NHARM2 = 4;

% 1. SETUP THE MATRIX OF PREDICTORS AA
% ====================================

dvec0 = datevec(dnobs(  1));
dvecF = datevec(dnobs(end));
dnum0 = datenum(dvec0(1), 01, 01);
tt = (dnobs - dnum0)/JDAYS;
tt = reshape(tt, numel(tt), 1);

%% Chose npoly based on number of years
%totyrs = dvecF(1) - dvec0(1) + 1;
%if (6 < totyrs), NPOLY = 3; end

AA = [];
for nn = 1:NPOLY
  AA = [AA, tt.^(nn-1)];
end

for nn = 1:NHARM2
  AA = [AA, sin(2*nn*pi*tt), cos(2*nn*pi*tt)];
end

size(AA)
[tt(1:2); tt(end)]

% 2. SOLVE FOR THE COEFFICIENTS
% =============================
zz   = AA(isok,:) \ omfin(isok);
zz
bias = AA*zz;
