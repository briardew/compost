%TCCON_COMPARE  Compare model results to TCCON retrievals (fast version)
%
%   This comparison is fast, but only updates the point in space the model is
%   evaluated at every 3 hours.  This is acceptable for stationary sites, but
%   is UNACCEPTABLE for AIRCRAFT & SATELLITE observations.

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%
% Changelog:
% 2019/11/18	New version
%
% TODO:
% * Replace constants with MAPL definitions
% * Pack variables into fsites structure instead of cell arrays


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

VARPS   = >>>VARPS<<<;					% Model name for surface pressure
VARQW   = >>>VARQW<<<;					% Model name for total water mixing ratio

DIROBS  = >>>DIROBS<<<;					% Directory holding obs files
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
disp('# TCCON comparison');
disp('#');
if (ISDRY)
  disp('# Treating model values as     DRY-air MOLE fractions ...');
  disp('#                              ---     ----              ');
else
  disp('# Treating model values as VIRTUAL-air MOLE fractions ...');
  disp('#                          -------     ----              ');
end
disp(sline);
disp(' ');

% Some error checking
% -------------------
if (isempty(VARQW) & ~ISDRY)
  error('No water vapor variable present to dry trace gas.');
end

if (DIROBS(end) ~= '/'), DIROBS = [DIROBS, '/']; end
if (DIRMOD(end) ~= '/'), DIRMOD = [DIRMOD, '/']; end

% Set up environment
% ------------------
atmosmug.constants;

% 1. READ OBS DATA
%==============================================================================%
disp('--- Observational data ---');
disp(['Using data in ', DIROBS, ' ...']);
disp(' ');

fsites = dir([DIROBS, VAROBS, '_*.nc']);
NSITES = numel(fsites);

cell_fobs  = cell(NSITES, 1);
cell_dnobs = cell(NSITES, 1);
cell_lat   = cell(NSITES, 1);
cell_lon   = cell(NSITES, 1);

cell_peavg   = cell(NSITES, 1);
cell_avgker  = cell(NSITES, 1);
cell_h2oapr  = cell(NSITES, 1);
cell_gasapr  = cell(NSITES, 1);
cell_xgasapr = cell(NSITES, 1);
cell_xgasobs = cell(NSITES, 1);
cell_xgaserr = cell(NSITES, 1);

for ic = 1:NSITES
  fobs = [fsites(ic).folder, '/', fsites(ic).name];

  fprintf(['Reading ', fsites(ic).name, ' ...\n']);

  cell_fobs{ic}  = fsites(ic).name;
  cell_dnobs{ic} = ncread(fobs, 'dnum');
  cell_lat{ic}   = double(ncread(fobs, 'lat'));
  cell_lon{ic}   = double(ncread(fobs, 'lon'));

  cell_peavg{ic}   =          ncread(fobs, 'peavg');
  cell_avgker{ic}  =  squeeze(ncread(fobs, 'avgker'));
  cell_h2oapr{ic}  =          ncread(fobs, 'priorh2o');
  cell_gasapr{ic}  = SCLOBS * ncread(fobs, 'priorpro');
  cell_xgasapr{ic} = SCLOBS * ncread(fobs, 'priorobs')';
  cell_xgasobs{ic} = SCLOBS * ncread(fobs, 'obs')';
  cell_xgaserr{ic} = SCLOBS * ncread(fobs, 'uncert')';
end
fprintf('\n');

NLAVG = size(cell_avgker{ic}, 1);


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

NLAT = numel(grdlat);
NLON = numel(grdlon);
[LA,  LO]       = meshgrid(grdlat, grdlon);
[LAz, LOz, KKz] = meshgrid(grdlat, grdlon, [1:NLEV]');

% Allocate output arrays and fill w/ nans
% ---------------------------------------
cell_psmod   = cell(NSITES, 1);
cell_xgasmod = cell(NSITES, 1);
cell_xgasalt = cell(NSITES, 1);
cell_tpwapr  = cell(NSITES, 1);
cell_tpwmod  = cell(NSITES, 1);

for ic = 1:NSITES
  cell_psmod{ic}   = NaN*ones(size(dnmod));
  cell_xgasmod{ic} = NaN*ones(size(dnmod));
  cell_xgasalt{ic} = NaN*ones(size(dnmod));
  cell_tpwapr{ic}  = NaN*ones(size(dnmod));
  cell_tpwmod{ic}  = NaN*ones(size(dnmod));
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
    dp  = 1e-2*atmosmug.getdp(ncread(fmet, VARPS), NLEV);
    qq  = zeros(size(gas));
    if (~isempty(VARQW)), qq = ncread(fmet, VARQW); end

%   Print message if read succeeds
    fprintf(['Computing values on ', date, ' at ', time, 'Z', ' ...\n']);
  catch ME
    continue;
  end

% Extend longitudinal dim for periodic interp
  gas = [gas; gas(1,:,:)];
  dp  = [ dp;  dp(1,:,:)];
  qq  = [ qq;  qq(1,:,:)];

% B. Interpolate to time and space of observations
% ------------------------------------------------
  for ic = 1:NSITES
    dnobs = cell_dnobs{ic};

%   Determine spatial location for interpolation
%   *** Only appropriate when obs location/pressure is nearly constant ***
    iob0 = find(dnnow-DSKIP <= dnobs, 1, 'first');
    iobF = find(dnobs <  dnnow+DSKIP, 1, 'last');
    if (isempty(iobF)), iobF = numel(dnobs); end
    iobs = [iob0:iobF]';

%   Skip sites that have no data for this time
    if (isempty(iobs)), continue; end

    lat = mean(cell_lat{ic}(iobs));
    lon = mean(cell_lon{ic}(iobs));

    peavg   = mean(cell_peavg{ic}(:,iobs),  2);
    avgker  = mean(cell_avgker{ic}(:,iobs), 2);
%   What's effect of mean?  Lots of scatter at tropical sites
    h2oapr  = mean(cell_h2oapr{ic}(:,iobs), 2);
    gasapr  = mean(cell_gasapr{ic}(:,iobs), 2);
    xgasapr = mean(cell_xgasapr{ic}(iobs));

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

%   Compute model specific humidity and trace gas on native levels
    qqnat   = wgtT*wgtR*squeeze(qq(iLO,iLA,:))   + ...
              wgtB*wgtR*squeeze(qq(iLO,iLA-1,:)) + ...
              wgtT*wgtL*squeeze(qq(iLO-1,iLA,:)) + ...
              wgtB*wgtL*squeeze(qq(iLO-1,iLA-1,:));
    gasnat  = wgtT*wgtR*squeeze(gas(iLO,iLA,:))   + ...
              wgtB*wgtR*squeeze(gas(iLO,iLA-1,:)) + ...
              wgtT*wgtL*squeeze(gas(iLO-1,iLA,:)) + ...
              wgtB*wgtL*squeeze(gas(iLO-1,iLA-1,:));

%
%   *** Once you settle this drying issue, you can move all this outside  ***
%   *** the site loop ***
%

%   Convert to dry-air mole fractions
%   ---------------------------------
%   ***           All GEOS-5 trace gas molar mixing ratios are            ***
%   ***           actually "virtual molar mixing ratios"                  ***
%   Thus, they are treated as rescalings of the mass mixing ratio by the
%   DRY-AIR molecular masses, i.e.:
%       1. XX_dry = XX_total / (1 - qq),         regardless of mass/molar
%       2. XX_mol = XX_mass * MAIR_DRY/MOLMASS   regardless of dry/total
    if (~ISDRY), gasnat = gasnat./(1 - qqnat); end

%   Interpolate to obs levels and apply averaging kernel
    gasmod  = interp1(plmod, gasnat, peavg, 'linear', gasnat(end));
    xgasmod = xgasapr + sum(avgker.*(gasmod - gasapr), 1)';

%   Wet model profile with model qv, then dry with prior qv
    qqmod   = interp1(plmod,  qqnat, peavg, 'linear',  qqnat(end));
    gasalt  = gasmod .* (1 - qqmod) ./ (1 - h2oapr);
    xgasalt = xgasapr + sum(avgker.*(gasalt - gasapr), 1)';

    cell_psmod{ic}(it)   = pemod(end);
    cell_xgasmod{ic}(it) = xgasmod;
    cell_xgasalt{ic}(it) = xgasalt;

    cell_tpwapr{ic}(it) = sum(0.5*(h2oapr(1:end-1) + h2oapr(2:end)).*diff(peavg));
    cell_tpwmod{ic}(it) = sum(0.5*( qqmod(1:end-1) +  qqmod(2:end)).*diff(peavg));
  end
end
fprintf('\n');

% Keep file sizes small
clearvars -except NSITES dnmod cell_*
