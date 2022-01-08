%OBSPACK_STATION_COMPARE  Compare model results to NOAA ObsPack station data
%
%   This comparison is fast, but only updates the point in space the model is
%   evaluated at every 3 hours.  This is acceptable for stationary sites, but
%   is UNACCEPTABLE for MOBILE (viz. AIRCRAFT, SHIPBOARD, and AIRCORE data).
%   For those data types, use the obspack_mobile_compare utility.

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%
% Changelog:
% 2019/03/18	New version
% 2019/04/18	Changed qcflags to all x's if variable doesn't exist
% 2019/11/13	Many changes, mostly zero-diff
%
% TODO:
% * Replace constants with MAPL definitions
% * Pack variables into fsites structure instead of cell arrays
% * Improve qcflag handling (non-trivial, varies across datasets)
% * Keep a record somehow of obs lons that get changed or tossed
% * Maybe save model elevations?


% "Virtual-air":
% --------------
% GEOS holds tracers as total-air mass mixing ratios multiplied by a constant,
% which is usually the ratio of molar masses of the gas and dry-air.  The
% result is something that looks like a total-air molar mixing ratio, but
% isn't.  We call such values "virtual-air" molar mixing ratios.
%
% Basic naming conventions:
% -------------------------
% lat, lon	Latitude, longitude at obs location
% alt, prs	Altitude, pressure     ""    ""
% pl, zl	Pressure, height of model at mid-layer points
% pe, ze	Pressure, height of model on layer edges
% qq		Sum of all water components (qv + ql + qi)
%
% gasnat	Model values on model levels
% gasmod	Model values on obs   levels
% gasobs	Obs   values on obs   levels
% gasapr	Prior values on obs   levels
%
% xgasmod       Model value with averaging kernel applied
% xgasobs	Obs   value for  averaging kernel applications
% xgasapr	Prior value for  averaging kernel applications
%
% NB: The pl, pe, zl, and ze variables only really make sense on the model
% grids.  If we are comparing to observed values, we change the names to prs
% and alt since the e and l extensions are to denote layer edges or centers.
%==============================================================================%

% User supplied definitions
% -------------------------
ISDRY   = >>>ISDRY<<<;					% Are model mixing ratios dry air or total?
SCLOBS  = >>>SCLOBS<<<;					% Scaling from   obs units to output units
SCLMOD  = >>>SCLMOD<<<;					% Scaling from model units to output units
VAROBS  = >>>VAROBS<<<;					% Obs name for trace gas
VARMOD  = >>>VARMOD<<<;					% Model name for trace gas

VARPHIS = >>>VARPHIS<<<;				% Model name for surface geopotential
VARPS   = >>>VARPS<<<;					% Model name for surface pressure
VARZL   = >>>VARZL<<<;					% Model name for geopotential height
VARQW   = >>>VARQW<<<;					% Model name for total water mixing ratio

DIROBS  = >>>DIROBS<<<;					% ObsPack directory
DIRMOD  = >>>DIRMOD<<<;					% Directory holding model files
HDMET   = >>>HDMET<<<;					% Head of model meteo files
HDGAS   = >>>HDGAS<<<;					% Head of model trace gas files
TTMET   = >>>TTMET<<<;					% Time format of model meteo files
TTGAS   = >>>TTGAS<<<;					% Time format of model trace gas files

NLEV  = 72;						% Number of model levels
DNUM0 = datenum(2000, 01, 01);				% Starting date of comparison (arbitrary)
DNUMF = datenum(2025, 12, 31);				% Ending   date of comparison (arbitrary)
DSKIP = 1/8;						% Fraction of day between files


% PRINT HEADER & SET UP ENVIRONMENT
%==============================================================================%
sline = ['#', repmat('=', 1, 79)];
disp(sline);
disp('# NOAA ObsPack comparison (STATION version)');
disp('#');
if (ISDRY)
  disp('# *** Treating model values as     DRY-air MOLE fractions ...');
  disp('#                                  ---     ----              ');
else
  disp('# *** Treating model values as VIRTUAL-air MOLE fractions ...');
  disp('                               -------     ----              ');
end
disp(sline);
disp(' ');

% Some error checking
% -------------------
if (isempty(VARQW) & ~ISDRY)
  error('No water vapor variable present to dry trace gas.');
end

if (isempty(VARZL) & isempty(VARPHIS))
  error('Need either level heights (VARZL) or surface geopotential (VARPHIS).');
end

if (DIROBS(end) ~= '/'), DIROBS = [DIROBS, '/']; end
if (DIRMOD(end) ~= '/'), DIRMOD = [DIRMOD, '/']; end

% Set up environment
% ------------------
addpath('/discover/nobackup/bweir/matlab/globutils');

% Replace with MAPL value load
CP = 1.0046e+3;
RD = 2.8705e+2;
RV = 4.6150e+2;
RDOVERCP = RD/CP;
EPS  = RD/RV;
KAP1 = RDOVERCP + 1;
KAPR = 1/RDOVERCP;
GRAV = 9.80665;

BPOW = 8.31447*0.0065/(GRAV*0.0289644);
BCON = 288.15/0.0065;


% 1. READ OBS DATA
%==============================================================================%
disp('--- Observational data ---');
disp(['Using data in ', DIROBS, ' ...']);
disp(' ');

contents  = dir([DIROBS, VAROBS, '_*.nc']);
cell_fobs = {contents(:).name};
NSITES    = numel(cell_fobs);

cell_dnobs   = cell(NSITES, 1);
cell_obflags = cell(NSITES, 1);
cell_qcflags = cell(NSITES, 1);
cell_opnums  = cell(NSITES, 1);

cell_lat     = cell(NSITES, 1);
cell_lon     = cell(NSITES, 1);
cell_alt     = cell(NSITES, 1);

cell_zsobs   = cell(NSITES, 1);
cell_gasobs  = cell(NSITES, 1);

subset = @(v,i) v(i);
for ic = 1:NSITES
  fobs = [DIROBS, cell_fobs{ic}];

% Print obs filename
  fprintf(['Reading ', cell_fobs{ic}, ' ...\n']);

% Pull out time values and convert to date fractions
  dvec  = ncread(fobs, 'time_components');
  year  = double(dvec(1,:));
  month = double(dvec(2,:));
  day   = double(dvec(3,:));
  hour  = double(dvec(4,:));
  minut = double(dvec(5,:));
  sec   = double(dvec(6,:));

  dnobs = datenum(year, month, day, hour, minut, sec)';
  lats  = double(ncread(fobs, 'latitude'));
  lons  = double(ncread(fobs, 'longitude'));

% Make sure all lons are in [-180,180)
% (*** This should throw a warning or keep a record or something ***)
  i180 = find(180 <= lons);
  lons(i180) = lons(i180) - 360;

% Toss out obs with bad lat/lon data
% (*** This should throw a warning or keep a record or something ***)
  iok = find(-1000 < lats & -1000 < lons);

  cell_dnobs{ic} = dnobs(iok);

% Pull out trace gas values
  vals = double(ncread(fobs, 'value'));
  vals(vals < 0) = NaN;
  cell_gasobs{ic} = SCLOBS * vals(iok);

% Pull out vertical info and flags
  cell_lat{ic} = lats(iok);
  cell_lon{ic} = lons(iok);

% Pull out altitude data
  try
    cell_alt{ic} = double(subset(ncread(fobs, 'altitude'),     iok));
  catch
    cell_alt{ic} = double(subset(ncread(fobs, 'gps_altitude'), iok));
  end

  cell_zsobs{ic} = NaN*ones(size(iok));
  try
    cell_zsobs{ic} = double(subset(ncread(fobs, 'elevation'), iok));
  end

  qcflags = repmat('x', [10, numel(vals)]);
  try
    qcf_in = ncread(fobs, 'qcflag');
    if (size(qcf_in) == size(qcflags)), qcflags = qcf_in; end
  end
  cell_qcflags{ic} = qcflags(:,iok);

  cell_obflags{ic} = NaN*ones(size(iok));
  try
    cell_obflags{ic} = subset(ncread(fobs, 'obs_flag'), iok);
  end

  cell_opnums{ic} = NaN*ones(size(iok));
  try
    cell_opnums{ic} = subset(ncread(fobs, 'obspack_num'), iok);
  end
end
fprintf('\n');

clear subset;


% 2. COMPUTE MODEL COMPARISON
%==============================================================================%
disp('--- Model comparison ---');

% Set up model time and grid variables
% ------------------------------------
dnmod = [DNUM0:DSKIP:DNUMF]';

contents = dir([DIRMOD, HDMET, '*.nc4']);
if (isempty(contents)), exit; end
fmet = [contents(1).folder, '/', contents(1).name];

grdlat = ncread(fmet, 'lat');
grdlon = ncread(fmet, 'lon');
% Extend longitudinal dim for periodic interp
grdlon = [grdlon; 180];

if (~isempty(VARPHIS))
  grdzs = 1/GRAV * ncread(fmet, VARPHIS);
% Extend longitudinal dim for periodic interp
  grdzs = [grdzs; grdzs(1,:)];
end

% Allocate output arrays and fill w/ nans
% ---------------------------------------
cell_prsmod = cell(NSITES, 1);
cell_qqmod  = cell(NSITES, 1);
cell_gasmod = cell(NSITES, 1);

for ic = 1:NSITES
  cell_prsmod{ic} = NaN*ones(size(dnmod));
  cell_qqmod{ic}  = NaN*ones(size(dnmod));
  cell_gasmod{ic} = NaN*ones(size(dnmod));
end

for it = 1:numel(dnmod)
  dnnow = dnmod(it);
  date  = datestr(dnnow, 'yyyymmdd');
  time  = datestr(dnnow, 'HH');

  fgas = [DIRMOD, HDGAS, date, '_', datestr(dnnow,TTGAS), 'z.nc4'];
  fmet = [DIRMOD, HDMET, date, '_', datestr(dnnow,TTMET), 'z.nc4'];

% A. Read model output, skipping missing/broken files
% ---------------------------------------------------
  try
    gas = SCLMOD*ncread(fgas, VARMOD);
    dp  = 1e-2*getdp(ncread(fmet, VARPS), NLEV);
    qq  = zeros(size(gas));
    if (~isempty(VARQW)), qq = ncread(fmet, VARQW); end
    if (~isempty(VARZL)), zl = ncread(fmet, VARZL); end

%   Print message if read succeeds
    fprintf(['Computing values on ', date, ' at ', time, 'Z', ' ...\n']);
  catch ME
    continue;
  end

% Extend longitudinal dim for periodic interp
  gas = [gas; gas(1,:,:)];
  dp  = [ dp;  dp(1,:,:)];
  qq  = [ qq;  qq(1,:,:)];
  if (~isempty(VARZL)), zl = [zl; zl(1,:,:)]; end

% B. Interpolate to time and space of observations
% ------------------------------------------------
  for ic = 1:NSITES
    dnobs = cell_dnobs{ic};

%   Skip sites that have no data for this time
    if (numel(dnobs) == 0), continue; end
    if (dnnow+DSKIP < dnobs(1) | dnobs(end) < dnnow-DSKIP), continue; end

%   Determine spatial location for interpolation
%   (This could use some improvement, but note that this program is only
%    appropriate when the obs location doesn't change much)
    [junk,iu] = unique(dnobs, 'stable');

    if (1 < numel(iu))
      alt = interp1(dnobs(iu), cell_alt{ic}(iu), dnnow, 'linear', 'extrap');
      lat = interp1(dnobs(iu), cell_lat{ic}(iu), dnnow, 'linear', 'extrap');
      lon = interp1(dnobs(iu), cell_lon{ic}(iu), dnnow, 'linear', 'extrap');
    else
      alt = cell_alt{ic}(iu);
      lat = cell_lat{ic}(iu);
      lon = cell_lon{ic}(iu);
    end

%   Determine spatial interpolation points and weights
    iLA  = find(lat < grdlat, 1);
    iLO  = find(lon < grdlon, 1);

    wgtB = (grdlat(iLA) - lat)/(grdlat(iLA) - grdlat(iLA-1));
    wgtT = 1 - wgtB;
    wgtL = (grdlon(iLO) - lon)/(grdlon(iLO) - grdlon(iLO-1));
    wgtR = 1 - wgtL;

%   Compute pressures at model interfaces (pemod) and layer centers
%   (plmod)
    pemod = 0.01 + wgtT*wgtR*cumsum(squeeze(dp(iLO,iLA,:)))   + ...
                   wgtB*wgtR*cumsum(squeeze(dp(iLO,iLA-1,:))) + ...
                   wgtT*wgtL*cumsum(squeeze(dp(iLO-1,iLA,:))) + ...
                   wgtB*wgtL*cumsum(squeeze(dp(iLO-1,iLA-1,:)));
    pemod = [0.01; pemod];
    plmod = (   (pemod(2:end).^KAP1 - pemod(1:end-1).^KAP1) ...
             ./ (KAP1*(pemod(2:end) - pemod(1:end-1)))      ).^KAPR;

%   Determine mid-layer heights
    if (~isempty(VARZL))
      zlmod = wgtT*wgtR*squeeze(zl(iLO,iLA,:))   + ...
              wgtB*wgtR*squeeze(zl(iLO,iLA-1,:)) + ...
              wgtT*wgtL*squeeze(zl(iLO-1,iLA,:)) + ...
              wgtB*wgtL*squeeze(zl(iLO-1,iLA-1,:));
    else
      zsmod = wgtT*wgtR*grdzs(iLO,iLA)   + ...
              wgtB*wgtR*grdzs(iLO,iLA-1) + ...
              wgtT*wgtL*grdzs(iLO-1,iLA) + ...
              wgtB*wgtL*grdzs(iLO-1,iLA-1);
      zlmod = zsmod + BCON*(1 - (plmod/plmod(end)).^BPOW);
    end

%   Compute model specific humidity and trace gas on native levels
    qqnat   = wgtT*wgtR*squeeze(qq(iLO,iLA,:))   + ...
              wgtB*wgtR*squeeze(qq(iLO,iLA-1,:)) + ...
              wgtT*wgtL*squeeze(qq(iLO-1,iLA,:)) + ...
              wgtB*wgtL*squeeze(qq(iLO-1,iLA-1,:));
    gasnat  = wgtT*wgtR*squeeze(gas(iLO,iLA,:))   + ...
              wgtB*wgtR*squeeze(gas(iLO,iLA-1,:)) + ...
              wgtT*wgtL*squeeze(gas(iLO-1,iLA,:)) + ...
              wgtB*wgtL*squeeze(gas(iLO-1,iLA-1,:));

%   Convert to dry-air mole fractions
%   ---------------------------------
%   ***           All GEOS-5 trace gas molar mixing ratios are            ***
%   ***           actually "virtual molar mixing ratios"                  ***
%   Thus, they are treated as rescalings of the mass mixing ratio by the
%   DRY-AIR molecular masses, i.e.:
%       1. XX_dry = XX_total / (1 - qq),         regardless of mass/molar
%       2. XX_mol = XX_mass * MAIR_DRY/MOLMASS   regardless of dry/total
    if (~ISDRY), gasnat = gasnat./(1 - qqnat); end

%   Interpolate to vertical location of obs
    prsmod = interp1(zlmod, plmod,  alt, 'linear', plmod(end));
    qqmod  = interp1(zlmod, qqnat,  alt, 'linear', qqnat(end));
    gasmod = interp1(zlmod, gasnat, alt, 'linear', gasnat(end));

    cell_prsmod{ic}(it) = prsmod;
    cell_qqmod{ic}(it)  = qqmod;
    cell_gasmod{ic}(it) = gasmod;
  end
end
fprintf('\n');

% Keep file sizes small
clearvars -except NSITES dnmod cell_*
