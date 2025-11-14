%ACOS_JOIN  Join model comparisons to lite file retrievals for different years

% Author(s):	Brad Weir <brad.weir@nasa.gov>
%
% Changelog:
% 2019-03-18	New version
%
% TODO:
% * Get rid of this; just use codas_compare.tmpl.m
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
    psmod(iyrm:end)    = vin1.psmod(iyrm:end);
    xco2mod(iyrm:end)  = vin1.xco2mod(iyrm:end);
    priormod(iyrm:end) = vin1.priormod(iyrm:end);
end

clear nyear vnew vin1 iyrm;
