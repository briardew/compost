%OBSPACK_JOIN  Join model comparisons to NOAA ObsPack data for different
%   years and different types

% Author(s):	Brad Weir <brad.weir@nasa.gov>
%
% Changelog:
% 2019-03-18	New version
% 2025-11-13	Added support for no mobile
%
% TODO:
%===============================================================================

HEADID = '>>>HEADID<<<';
YEAR0  = >>>YEAR0<<<;
YEARF  = >>>YEARF<<<;

FIXHEAD = [HEADID, '_station.'];
MOBHEAD = [HEADID, '_mobile.'];

% Fill the local workspace with the variables we need
load([FIXHEAD, num2str(YEAR0), '.mat']);

vfix = load([FIXHEAD, num2str(YEAR0), '.mat']);
vmob = [];
try, vmob = load([MOBHEAD, num2str(YEAR0), '.mat']); end

% Fill in subsequent years
for nyear = YEAR0+1:YEARF
    vfix1 = load([FIXHEAD, num2str(nyear), '.mat']);
    vmob1 = [];
    try, vmob1 = load([MOBHEAD, num2str(nyear), '.mat']); end

    iyrm = find(dnmod == datenum(nyear,01,01));

    for ic = 1:NSITES
        % Station data
        vfix.cell_prsmod{ic}(iyrm:end) = vfix1.cell_prsmod{ic}(iyrm:end);
        vfix.cell_qqmod{ic}(iyrm:end)  = vfix1.cell_qqmod{ic}(iyrm:end);
        vfix.cell_gasmod{ic}(iyrm:end) = vfix1.cell_gasmod{ic}(iyrm:end);

        % Mobile data
        if ~isempty(vmob) && ~isempty(vmob1)
            iyro = find(vmob.cell_dnobs{ic} < datenum(nyear,01,01), 1, 'last') + 1;
            if isempty(iyro), iyro = 1; end

            vmob.cell_prsmod{ic}(iyro:end) = vmob1.cell_prsmod{ic}(iyro:end);
            vmob.cell_qqmod{ic}(iyro:end)  = vmob1.cell_qqmod{ic}(iyro:end);
            vmob.cell_gasmod{ic}(iyro:end) = vmob1.cell_gasmod{ic}(iyro:end);
        end
    end
end

clear nyear vnew vfix1 vmob1 iyrm iyro;

% Interpolate and merge
for ic = 1:NSITES
    dnobs = cell_dnobs{ic};

    % Start with station data
    cell_prsmod{ic} = interp1(dnmod, vfix.cell_prsmod{ic}, dnobs);
    cell_qqmod{ic}  = interp1(dnmod, vfix.cell_qqmod{ic},  dnobs);
    cell_gasmod{ic} = interp1(dnmod, vfix.cell_gasmod{ic}, dnobs);

    % Use mobile data where available
    if ~isempty(vmob) && 0 < nnz(~isnan(vmob.cell_gasmod{ic}))
        cell_prsmod{ic} = vmob.cell_prsmod{ic};
        cell_qqmod{ic}  = vmob.cell_qqmod{ic};
        cell_gasmod{ic} = vmob.cell_gasmod{ic};
    end
end

clear vfix vmob;
