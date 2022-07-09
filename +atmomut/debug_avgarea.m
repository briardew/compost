%DEBUG_AVGAREA  Performs a series of tests to make sure avgarea is working
%   correctly
%
%   1. HI -> LO conserves area averages
%   2. LO -> HI has no ringing
%   3. HI -> LO has no ringing
%   4. A step function stays between 0 and 1
%      (Equivalent to positivity of smoothing matrices)
%
%   Note: There is no reason the expect non-nested grids to pass the
%   LO -> HI -> LO test since the area averaging has to decide how to
%   distribute values across cells.  This can be accomplished by making one
%   transformation the Moore-Penrose pseudoinverse of the other.  However,
%   this transformation has negative values, which will then not preserve
%   positivity.  AFAIK, there is no pseudoinverse with all positive entries,
%   but cannot prove it.

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%==============================================================================%

addpath('/discover/nobackup/bweir/matlab/nanstats');

%DHI  = 6;
%DLO  = 12;
DHI  = 3;
DLO  = 6;
%DHI  = 0.25;
%DLO  = 0.5;

lon  = [-180:8:172]';
NLON = numel(lon);

% Output some junk
sbar = repmat('=', 1, 80);
disp(sbar);
disp('Area averaging debugging script');
disp(sbar);
disp(' ');

% 1. POLE NODE TO POLE NODE UPSCALING
% -----------------------------------
disp('1. Pole   NODE to pole   NODE');
disp('-----------------------------');

lathi = [-90+DHI/2:DHI:90-DHI/2]';
latlo = [-90+DLO/2:DLO:90-DLO/2]';

lonhi = lon;
lonlo = lon;

diag_avgarea;

% 2. POLE NODE TO POLE CENTER UPSCALING
% -------------------------------------
disp('2. Pole   NODE to pole CENTER');
disp('-----------------------------');

lathi = [-90+DHI/2:DHI:90-DHI/2]';
latlo = [-90:DLO:90]';

lonhi = lon;
lonlo = lon;

diag_avgarea;

% 3. POLE CENTER TO POLE NODE UPSCALING
% -------------------------------------
disp('3. Pole CENTER to pole   NODE');
disp('-----------------------------');

lathi = [-90:DHI:90]';
latlo = [-90+DLO/2:DLO:90-DLO/2]';

lonhi = lon;
lonlo = lon;

diag_avgarea;

% 4. POLE CENTER TO POLE NODE UPSCALING
% -------------------------------------
disp('4. Pole CENTER to pole CENTER');
disp('-----------------------------');

lathi = [-90:DHI:90]';
latlo = [-90:DLO:90]';

lonhi = lon;
lonlo = lon;

diag_avgarea

% 5. DIFFERENT LON SPACINGS
% -------------------------
disp('5. Different lon spacings');
disp('-------------------------');

lathi = [-90+DLO/2:DLO:90-DLO/2]';
latlo = [-90+DLO/2:DLO:90-DLO/2]';

lonhi = [-178:4:178]';
lonlo = [-176:8:176]';

diag_avgarea;

% 6. DIFFERENT LAT & LON SPACINGS
% -------------------------------
disp('6. Different lat & lon spacings');
disp('-------------------------------');

lathi = [-90+DHI/2:DHI:90-DHI/2]';
latlo = [-90+DLO/2:DLO:90-DLO/2]';

lonhi = [-178:4:178]';
lonlo = [-176:8:176]';

diag_avgarea;
