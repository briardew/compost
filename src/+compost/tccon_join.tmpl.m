%TCCON_JOIN  Join model comparisons to TCCON retrievals for different years

% Author(s):	Brad Weir <brad.weir@nasa.gov>
%
% Changelog:
% 2019-03-18	New version
%
% TODO:
%===============================================================================

HEADID = '>>>HEADID<<<';
YEAR0  = >>>YEAR0<<<;
YEARF  = >>>YEARF<<<;

FHEAD  = [HEADID, '.'];

% Fill the local workspace with the variables we need
load([FHEAD, num2str(YEAR0), '.mat']);

% Fill in subsequent years
for nyear = YEAR0+1:YEARF
    vin1 = load([FHEAD, num2str(nyear), '.mat']);

    iyrm = find(dnmod == datenum(nyear,01,01));

    % Add in this year's data
    for ic = 1:NSITES
        cell_psmod{ic}(iyrm:end)   = vin1.cell_psmod{ic}(iyrm:end);
        cell_xgasmod{ic}(iyrm:end) = vin1.cell_xgasmod{ic}(iyrm:end);
        cell_xgasalt{ic}(iyrm:end) = vin1.cell_xgasalt{ic}(iyrm:end);
    end
end

clear nyear vnew vin1 iyrm;

% Interpolate to obs time step
for ic = 1:NSITES
    dnobs = cell_dnobs{ic};

    cell_psmod{ic}   = interp1(dnmod, cell_psmod{ic},   dnobs);
    cell_xgasmod{ic} = interp1(dnmod, cell_xgasmod{ic}, dnobs);
    cell_xgasalt{ic} = interp1(dnmod, cell_xgasalt{ic}, dnobs);
end
