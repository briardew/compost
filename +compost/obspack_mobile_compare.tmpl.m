%OBSPACK_MOBILE_COMPARE  Compare model results to NOAA ObsPack mobile data
%
%   This comparison is very slow because it updates the point in space the model
%   is evaluated at with every new observation.  This is necessary for mobile
%   (viz. aircraft and shipboard) observations, but is unncessary for stationary
%   sites.  For stationary sites, use the obspack_station_compare utility.

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%
% Changelog:
% 2019/03/18	New version
% 2019/04/18	Changed qcflags to all x's if variable doesn't exist
% 2019/04/25	Tweaks to handle extrapolation for measurements lower than
%		model altitude
% 2019/11/13	Many changes, mostly zero-diff
% 2023/02/15	Adding capability to read packed, daily files (e.g., MERRA-2 GMI)
%		Includes file format change, update to shell script
%
% TODO:
% * Replace constants with MAPL definitions
% * Pack variables into structure instead of cell arrays
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
DSKIP   = >>>DSKIP<<<;					% Fraction of day between files

NLEV  = 72;						% Number of model levels
DNUM0 = datenum(2000, 01, 01);				% Starting date of comparison (arbitrary)
DNUMF = datenum(2025, 12, 31);				% Ending   date of comparison (arbitrary)


% PRINT HEADER & SET UP ENVIRONMENT
%==============================================================================%
sline = ['#', repmat('=', 1, 79)];
disp(sline);
disp('# NOAA ObsPack comparison (MOBILE version)');
disp('#');
if (ISDRY)
  disp('# *** Treating model values as     DRY-air MOLE fractions ...');
  disp('#                                  ---     ----              ');
else
  disp('# *** Treating model values as VIRTUAL-air MOLE fractions ...');
  disp('#                              -------     ----              ');
end
disp('#');
disp('# *** Only comparing to aircraft, aircore, and shipboard data');
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
atmomut.constants;

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


% 2. FLATTEN
%==============================================================================%
array_dnobs  = [];
array_lat    = [];
array_lon    = [];
array_alt    = [];
array_inds   = [];
array_sitens = [];
array_gasobs = [];	% Not necessary, here for debugging

for ic = 1:NSITES
  if (isempty(strfind(cell_fobs{ic}, 'aircraft'))  & ...
      isempty(strfind(cell_fobs{ic}, 'aircore'))   & ...
      isempty(strfind(cell_fobs{ic}, 'shipboard')))
    continue;
  end

  array_dnobs  = [array_dnobs;  cell_dnobs{ic}];
  array_lat    = [array_lat;    cell_lat{ic}];
  array_lon    = [array_lon;    cell_lon{ic}];
  array_alt    = [array_alt;    cell_alt{ic}];

  nobs = numel(cell_dnobs{ic});
  array_inds   = [array_inds;   [1:nobs]'];
  array_sitens = [array_sitens; ic*ones(nobs, 1)];

  array_gasobs = [array_gasobs; cell_gasobs{ic}];
end

% Sort by obs time
[yy,inds] = sort(array_dnobs);

array_dnobs  = array_dnobs(inds);
array_lat    = array_lat(inds);
array_lon    = array_lon(inds);
array_alt    = array_alt(inds);
array_inds   = array_inds(inds);
array_sitens = array_sitens(inds);
array_gasobs = array_gasobs(inds);


% 3. COMPUTE MODEL COMPARISON
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

[LA,  LO]       = meshgrid(grdlat, grdlon);
[LAz, LOz, KKz] = meshgrid(grdlat, grdlon, [1:NLEV]');

if (~isempty(VARPHIS))
  grdzsin = 1/GRAV * ncread(fmet, VARPHIS);
% Hack to deal with packed daily files
  if (size(grdzsin,3) == 1)
    grdzs = grdzsin;
  else
    grdzs = grdzsin(:,:,1);
  end
% Extend longitudinal dim for periodic interp
  grdzs = [grdzs; grdzs(1,:)];
end

% Allocate output arrays and fill w/ nans
% ---------------------------------------
array_prsmod = NaN*ones(size(array_dnobs));
array_qqmod  = NaN*ones(size(array_dnobs));
array_gasmod = NaN*ones(size(array_dnobs));

lread = 1;				% Switch to minimize reads
dnprv = Inf;				% Hack to skip first iteration
for it = 1:numel(dnmod)
  dnnow = dnmod(it);
  date  = datestr(dnnow, 'yyyymmdd');

  fgas = [DIRMOD, HDGAS, date, datestr(dnnow,TTGAS), '.nc4'];
  fmet = [DIRMOD, HDMET, date, datestr(dnnow,TTMET), '.nc4'];

% A. Read model output, skipping missing/broken files
% ---------------------------------------------------
  if (lread == 1), try
    gasin = SCLMOD*ncread(fgas, VARMOD);
    psin  = ncread(fmet, VARPS);
    qqin  = zeros(size(gasin));
    zlin  = zeros(size(gasin));
    if (~isempty(VARQW)), qqin = ncread(fmet, VARQW); end
    if (~isempty(VARZL)), zlin = ncread(fmet, VARZL); end

%   Print daily message if read succeeds
    if (floor(dnnow) == dnnow)
      fprintf(['Computing values for ', date, ' ...\n']);
    end
  catch ME
    continue;
  end, end

% Hack to deal with packed daily files
  if (size(gasin,4) == 1)
    gasnow = gasin;
    psnow  = psin;
    qqnow  = qqin;
    zlnow  = zlin;
  else
    nn  = round((dnnow - floor(dnnow))/DSKIP) + 1;
    gasnow = gasin(:,:,:,nn);
    psnow  =  psin(:,:,  nn);
    qqnow  =  qqin(:,:,:,nn);
    zlnow  =  zlin(:,:,:,nn);

%   Try to minimize reads
    lread = 0;
    if (floor(dnnow + DSKIP) == dnnow + DSKIP), lread = 1; end
  end
  dpnow = 1e-2*atmomut.getdp(psnow, NLEV);

% Extend longitudinal dim for periodic interp
  gasnow = [gasnow; gasnow(1,:,:)];
  dpnow  = [ dpnow;  dpnow(1,:,:)];
  qqnow  = [ qqnow;  qqnow(1,:,:)];
  zlnow  = [ zlnow;  zlnow(1,:,:)];

% Compute pressures at model interfaces (penow) and layer centers (plnow)
  penow = 0.01 + cumsum(dpnow, 3);
  penow = cat(3, 0.01*ones(numel(grdlon), numel(grdlat)), penow);
  plnow = (   (penow(:,:,2:end).^KAP1 - penow(:,:,1:end-1).^KAP1) ...
           ./ (KAP1*(penow(:,:,2:end) - penow(:,:,1:end-1)))      ).^KAPR;

% Use barometric formula for mid-layer heights if not provided
  if (isempty(VARZL))
     zlnow = zeros(size(plnow));
     for kk = 1:size(zlnow,3)
       zlnow(:,:,kk) = grdzs + BCON*(1 - (plnow(:,:,kk)./plnow(:,:,end)).^BPOW);
     end
  end

% Convert to dry-air mole fractions
% ---------------------------------
% ***           All GEOS-5 trace gas molar mixing ratios are            ***
% ***                   "virtual molar mixing ratios"                   ***
% They are treated as rescalings of the mass mixing ratio by the
% DRY-AIR molecular masses, i.e.:
%     1. XX_dry = XX_total / (1 - qq),         regardless of mass/molar
%     2. XX_mol = XX_mass * MAIR_DRY/MOLMASS   regardless of dry/total
  if (~ISDRY), gasnow = gasnow./(1 - qqnow); end

% B. Interpolate to time and space of obs
% ---------------------------------------
% Determine first (iob0) and last (iobF) obs indices bounded by model date
% numbers dnprv and dnnow (NB: to troubleshoot, run datestr on elements of
% array_dnobs)
  iob0 = find(dnprv <= array_dnobs, 1, 'first');
  iobF = find(array_dnobs <  dnnow, 1, 'last');
  iobs = [iob0:iobF]';

  if (~isempty(iobs))
    alphas = (dnnow - array_dnobs(iobs))/(dnnow - dnprv);

    nobs = numel(iobs);
    lats = array_lat(iobs);
    lons = array_lon(iobs);
    alts = array_alt(iobs);

%   Interpolate model altitudes (zlprv, zlnow) to obs lats and lons
    zlout0 = zeros(nobs, NLEV);
    zlout1 = zeros(nobs, NLEV);
    for kk = 1:NLEV
      zlout0(:,kk) = interp2(LA, LO, zlprv(:,:,kk), lats, lons);
      zlout1(:,kk) = interp2(LA, LO, zlnow(:,:,kk), lats, lons);
    end
    zlout0 = zlout0';
    zlout1 = zlout1';

%   Translate altitudes (alts) to fractions indicating level index (kkobs)
    kkobs0 = zeros(nobs, 1);
    kkobs1 = zeros(nobs, 1);
    for ii = 1:nobs
      zlext0 = [1e9; zlout0(:,ii); -1e9];
      zlext1 = [1e9; zlout1(:,ii); -1e9];
      kkext  = [  1,       1:NLEV, NLEV]';

      kkobs0(ii) = interp1(zlext0, kkext, alts(ii));
      kkobs1(ii) = interp1(zlext1, kkext, alts(ii));
    end

%   Do a 3d interp on now and prv variables
    prsmod0 = interp3(LAz, LOz, KKz, plprv, lats, lons, kkobs0);
    prsmod1 = interp3(LAz, LOz, KKz, plnow, lats, lons, kkobs1);
    prsmods = alphas.*prsmod0 + (1 - alphas).*prsmod1;

    qqmod0  = interp3(LAz, LOz, KKz, qqprv, lats, lons, kkobs0);
    qqmod1  = interp3(LAz, LOz, KKz, qqnow, lats, lons, kkobs1);
    qqmods  = alphas.*qqmod0 + (1 - alphas).*qqmod1;

    gasmod0 = interp3(LAz, LOz, KKz, gasprv, lats, lons, kkobs0);
    gasmod1 = interp3(LAz, LOz, KKz, gasnow, lats, lons, kkobs1);
    gasmods = alphas.*gasmod0 + (1 - alphas).*gasmod1;

    array_prsmod(iobs) = prsmods;
    array_qqmod(iobs)  = qqmods;
    array_gasmod(iobs) = gasmods;
  end

% Save fields for next iteration
  dnprv  = dnnow;
  plprv  = plnow;
  zlprv  = zlnow;
  qqprv  = qqnow;
  gasprv = gasnow;
end
fprintf('\n');


% 4. UN-FLATTEN
%==============================================================================%
cell_prsmod = cell(NSITES, 1);
cell_qqmod  = cell(NSITES, 1);
cell_gasmod = cell(NSITES, 1);

for ic = 1:NSITES
  cell_prsmod{ic} = NaN*ones(size(cell_dnobs{ic}));
  cell_qqmod{ic}  = NaN*ones(size(cell_dnobs{ic}));
  cell_gasmod{ic} = NaN*ones(size(cell_dnobs{ic}));
end

for ii = 1:numel(array_dnobs)
  ic = array_sitens(ii);
  jj = array_inds(ii);

  cell_prsmod{ic}(jj) = array_prsmod(ii);
  cell_qqmod{ic}(jj)  = array_qqmod(ii);
  cell_gasmod{ic}(jj) = array_gasmod(ii);
end

% Keep file sizes small
clearvars -except NSITES dnmod cell_*
