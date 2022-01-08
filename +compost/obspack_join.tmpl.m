%OBSPACK_JOIN  Join model comparisons to NOAA ObsPack data for different
%   years and different types

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%
% Changelog:
% 2019/03/18	New version
%
% TODO:
%==============================================================================%

HEADID = '>>>HEADID<<<';
YEAR0  = >>>YEAR0<<<;
YEARF  = >>>YEARF<<<;

LOHEAD = [HEADID, '_station__'];
HIHEAD = [HEADID, '_mobile__'];

% Fill the local workspace with the variables we need
load([HIHEAD, num2str(YEAR0), '.mat']);

vlo = load([LOHEAD, num2str(YEAR0), '.mat']);
vhi = load([HIHEAD, num2str(YEAR0), '.mat']);

% Fill in subsequent years
for nyear = YEAR0+1:YEARF
  vlo1 = load([LOHEAD, num2str(nyear), '.mat']);
  vhi1 = load([HIHEAD, num2str(nyear), '.mat']);

  iyrm = find(dnmod == datenum(nyear,01,01));

  for ic = 1:NSITES
%   a. Add in station data
    vlo.cell_prsmod{ic}(iyrm:end) = vlo1.cell_prsmod{ic}(iyrm:end);
    vlo.cell_qqmod{ic}(iyrm:end)  = vlo1.cell_qqmod{ic}(iyrm:end);
    vlo.cell_gasmod{ic}(iyrm:end) = vlo1.cell_gasmod{ic}(iyrm:end);

%   b. Add in mobile data
    iyro = find(vhi.cell_dnobs{ic} < datenum(nyear,01,01), 1, 'last') + 1;
    if (isempty(iyro)), iyro = 1; end

    vhi.cell_prsmod{ic}(iyro:end) = vhi1.cell_prsmod{ic}(iyro:end);
    vhi.cell_qqmod{ic}(iyro:end)  = vhi1.cell_qqmod{ic}(iyro:end);
    vhi.cell_gasmod{ic}(iyro:end) = vhi1.cell_gasmod{ic}(iyro:end);
  end
end

clear nyear vnew vlo1 vhi1 iyrm iyro;

% Start with lo-res data, use hi-res data for aircraft
for ic = 1:NSITES
  dnobs = cell_dnobs{ic};

  cell_prsmod{ic} = interp1(dnmod, vlo.cell_prsmod{ic}, dnobs);
  cell_qqmod{ic}  = interp1(dnmod, vlo.cell_qqmod{ic},  dnobs);
  cell_gasmod{ic} = interp1(dnmod, vlo.cell_gasmod{ic}, dnobs);

  if (nnz(~isnan(vhi.cell_gasmod{ic})) > 0)
    cell_prsmod{ic} = vhi.cell_prsmod{ic};
    cell_qqmod{ic}  = vhi.cell_qqmod{ic};
    cell_gasmod{ic} = vhi.cell_gasmod{ic};
  end
end

clear vlo vhi;
