%LITE_COMPARE  Compare model results to retrievals of XCO2 in the lite-file
%   format, viz., ACOS-GOSAT, OCO, and BESD-SCIAMACHY

% Author(s):	Brad Weir (brad.weir@nasa.gov)
%
% Changelog:
% 2019/11/18	New version
%
% TODO:
% * Replace constants with MAPL definitions


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
SCLOBS  = >>>SCLOBS<<<;					% Scaling from mol/mol to output units
SCLMOD  = >>>SCLMOD<<<;					% Scaling from mol/mol to output units
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
DSKIP   = >>>DSKIP<<<;					% Fraction of day between files

DNUM0 = datenum(2000, 01, 01);				% Starting date of comparison (arbitrary)
DNUMF = datenum(2025, 12, 31);				% Ending   date of comparison (arbitrary)


% PRINT HEADER & SET UP ENVIRONMENT
%==============================================================================%
sline = ['#', repmat('=', 1, 79)];
disp(sline);
disp('# LITE FILE (ACOS-GOSAT, OCO, and BESD-SCIAMACHY) comparison');
disp('#');
if (ISDRY)
  disp('# Treating model values as     DRY-air MOLE fractions ...');
  disp('#                              ---     ----              ');
else
  disp('# Treating model values as VIRTUAL-air MOLE fractions ...');
  disp('                           -------     ----              ');
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
atmomut.constants;

% Set up model time and space grid information
% --------------------------------------------
dnmod = [DNUM0:DSKIP:DNUMF]';

contents = dir([DIRMOD, HDMET, '*.nc4']);
if (isempty(contents)), exit; end
fmet = [contents(1).folder, '/', contents(1).name];

grdlat = ncread(fmet, 'lat');
grdlon = ncread(fmet, 'lon');
% Extend lat & lon for interp
grdlat = [grdlat(1) - (grdlat(2) - grdlat(1)); grdlat; ...
    grdlat(end) + (grdlat(end) - grdlat(end-1))];
grdlon = [grdlon(1) - (grdlon(2) - grdlon(1)); grdlon; ...
    grdlon(end) + (grdlon(end) - grdlon(end-1))];

NLAT = numel(grdlat);
NLON = numel(grdlon);
NLEV = 72;

% Determine read range for obs
% ----------------------------
fmet  = [contents(  1).folder, '/', contents(  1).name];
date0 = ncreadatt(fmet, 'time', 'begin_date');
dnin0 = datenum(num2str(date0), 'yyyymmdd') - 1;
fmet  = [contents(end).folder, '/', contents(end).name];
dateF = ncreadatt(fmet, 'time', 'begin_date');
dninF = datenum(num2str(dateF), 'yyyymmdd') + 1;

% 1. READ OBS DATA
%==============================================================================%
disp('--- Observational data ---');
disp(['Using data in ', DIROBS, ' ...']);
disp(' ');

% Set up arrays
% -------------
dnobs    = [];
obslat   = [];
obslon   = [];
xco2obs  = [];
flags    = [];
uncert   = [];
priorobs = [];
priorlev = [];
avgker   = [];
peavgs   = [];
psobs    = [];

foots  = [];
modes  = [];
surfts = [];
psapr  = [];
szas   = [];

aoddu  = [];
aodss  = [];
aodoc  = [];
aodbc  = [];
aodsu  = [];
aodql  = [];
aodqi  = [];
aodtot = [];

% Read 6-hourly data
% ------------------
for dnin = dnin0:6/24:dninF
  dvec  = datevec(dnin);
  flist = dir([DIROBS, '/Y', num2str(dvec(1)), '/*_', ...
               datestr(dvec,'yyyymmdd_hh'), 'z.nc*']);
  if (isempty(flist)), continue; end
% If there are multiple files, we only pick the newest
  fin = [flist(1).folder, '/', flist(1).name];

% Threshold ("necessary") variables
  dsecs_in = ncread(fin, 'time');
  dnobs_in = datenum(1970,01,01) + double(dsecs_in)/60./60./24.;

  obslat_in   = ncread(fin, 'latitude');
  obslon_in   = ncread(fin, 'longitude');
  xco2obs_in  = ncread(fin, 'xco2_final');
  flags_in    = ncread(fin, 'qcflag');
  uncert_in   = ncread(fin, 'xco2_uncert');
  priorobs_in = ncread(fin, 'xco2_apriori');
  priorlev_in = ncread(fin, 'co2_profile_apriori')';
  avgker_in   = ncread(fin, 'xco2_avgker')';
  peavgs_in   = ncread(fin, 'pressure_levels')';
  psobs_in    = max(peavgs_in(:,1), peavgs_in(:,end));

  dnobs    = [dnobs;    dnobs_in];
  obslat   = [obslat;   obslat_in];
  obslon   = [obslon;   obslon_in];
  xco2obs  = [xco2obs;  xco2obs_in];
  flags    = [flags;    flags_in];
  uncert   = [uncert;   uncert_in];
  priorobs = [priorobs; priorobs_in];
  priorlev = [priorlev; priorlev_in];
  avgker   = [avgker;   avgker_in];
  peavgs   = [peavgs;   peavgs_in];
  psobs    = [psobs;    psobs_in];

% Baseline ("unnecessary") variables -- retrieval
  foots_in  = NaN*ones(size(dnobs_in));
  modes_in  = NaN*ones(size(dnobs_in));
  surfts_in = NaN*ones(size(dnobs_in));
  psapr_in  = NaN*ones(size(dnobs_in));
  szas_in   = NaN*ones(size(dnobs_in));

  try, foots_in  = ncread(fin, 'footprint');          end
  try, modes_in  = ncread(fin, 'operation_mode');     end
  try, surfts_in = ncread(fin, 'surface_type');       end
  try, psapr_in  = ncread(fin, 'psurf_apriori');      end
  try, szas_in   = ncread(fin, 'solar_zenith_angle'); end

  foots  = [foots;  foots_in];
  modes  = [modes;  modes_in];
  surfts = [surfts; surfts_in];
  psapr  = [psapr;  psapr_in];
  szas   = [szas;   szas_in];

% Baseline ("unnecessary") variables -- aerosols
  aoddu_in  = NaN*ones(size(dnobs_in));
  aodss_in  = NaN*ones(size(dnobs_in));
  aodoc_in  = NaN*ones(size(dnobs_in));
  aodbc_in  = NaN*ones(size(dnobs_in));
  aodsu_in  = NaN*ones(size(dnobs_in));
  aodql_in  = NaN*ones(size(dnobs_in));
  aodqi_in  = NaN*ones(size(dnobs_in));
  aodtot_in = NaN*ones(size(dnobs_in));

  try, aoddu_in  = ncread(fin, 'aod_dust');    end
  try, aodss_in  = ncread(fin, 'aod_seasalt'); end
  try, aodoc_in  = ncread(fin, 'aod_oc');      end
  try, aodbc_in  = ncread(fin, 'aod_bc');      end
  try, aodsu_in  = ncread(fin, 'aod_sulfate'); end
  try, aodql_in  = ncread(fin, 'aod_water');   end
  try, aodqi_in  = ncread(fin, 'aod_ice');     end
  try, aodtot_in = ncread(fin, 'aod_total');   end

  aoddu  = [aoddu;  aoddu_in];
  aodss  = [aodss;  aodss_in];
  aodoc  = [aodoc;  aodoc_in];
  aodbc  = [aodbc;  aodbc_in];
  aodsu  = [aodsu;  aodsu_in];
  aodql  = [aodql;  aodql_in];
  aodqi  = [aodqi;  aodqi_in];
  aodtot = [aodtot; aodtot_in];

  fprintf(['Reading ', flist(1).name, ' ...\n']);
end
fprintf('\n');

if (numel(dnobs) == 0), return; end

% Paranoia
[~,inds] = sort(dnobs);

dnobs    = dnobs(inds);
obslat   = obslat(inds);
obslon   = obslon(inds);
xco2obs  = xco2obs(inds);

flags    = flags(inds);
uncert   = uncert(inds);
priorobs = priorobs(inds);
priorlev = priorlev(inds,:);
avgker   = avgker(inds,:);
peavgs   = peavgs(inds,:);
psobs    = psobs(inds);

foots  = foots(inds);
modes  = modes(inds);
surfts = surfts(inds);
psapr  = psapr(inds);
szas   = szas(inds);

aoddu  = aoddu(inds);
aodss  = aodss(inds);
aodoc  = aodoc(inds);
aodbc  = aodbc(inds);
aodsu  = aodsu(inds);
aodql  = aodql(inds);
aodqi  = aodqi(inds);
aodtot = aodtot(inds);


% 2. COMPUTE MODEL COMPARISON
%==============================================================================%
disp('--- Model comparison ---');

psmod    = NaN*ones(size(dnobs));
tco2mod  = NaN*ones(size(dnobs));
tco2use  = NaN*ones(size(dnobs));
xco2mod  = NaN*ones(size(dnobs));
priormod = NaN*ones(size(dnobs));

% Hack for first iteration
dnprv  = Inf;
dpprv  = NaN*ones(NLON, NLAT, NLEV);
qqprv  = NaN*ones(NLON, NLAT, NLEV);
co2prv = NaN*ones(NLON, NLAT, NLEV);

for it = 1:numel(dnmod)
  dnnow = dnmod(it);
  date  = datestr(dnnow, 'yyyymmdd');
  time  = datestr(dnnow, 'HH');

  fgas = [DIRMOD, HDGAS, date, datestr(dnnow,TTGAS), '.nc4'];
  fmet = [DIRMOD, HDMET, date, datestr(dnnow,TTMET), '.nc4'];

% A. Read model output, skipping missing/broken files
% ---------------------------------------------------
  try
    co2now = SCLMOD*ncread(fgas, VARMOD);
    dpnow  = 1e-2*atmomut.getdp(ncread(fmet, VARPS), NLEV);
    qqnow  = zeros(size(co2now));
    if (~isempty(VARQW)), qqnow = ncread(fmet, VARQW); end

%   Print message if read succeeds
    fprintf(['Computing values on ', date, ' at ', time, 'Z', ' ...\n']);
  catch ME
    continue;
  end

% Extend lat for constant interp
  co2now = cat(2, co2now(:,1,:), co2now, co2now(:,end,:));
  dpnow  = cat(2,  dpnow(:,1,:),  dpnow,  dpnow(:,end,:));
  qqnow  = cat(2,  qqnow(:,1,:),  qqnow,  qqnow(:,end,:));
% Extend lon for periodic interp
  co2now = cat(1, co2now(end,:,:), co2now, co2now(1,:,:));
  dpnow  = cat(1,  dpnow(end,:,:),  dpnow,  dpnow(1,:,:));
  qqnow  = cat(1,  qqnow(end,:,:),  qqnow,  qqnow(1,:,:));


% Convert to dry-air mole fractions
% ---------------------------------
% ***           All GEOS-5 trace gas molar mixing ratios are            ***
% ***                   "virtual molar mixing ratios"                   ***
% They are treated as rescalings of the mass mixing ratio by the
% DRY-AIR molecular masses, i.e.:
%     1. XX_dry = XX_total / (1 - qq),         regardless of mass/molar
%     2. XX_mol = XX_mass * MAIR_DRY/MOLMASS   regardless of dry/total
  if (~ISDRY), co2now = co2now./(1 - qqnow); end

% Apply the observation operator
  iob0 = find(dnprv <= dnobs, 1, 'first');
  iobF = find(dnobs <  dnnow, 1, 'last');
  nobs = iobF - iob0 + 1;

  iLOs = zeros(nobs, 1);
  iLAs = zeros(nobs, 1);

  wgtBs = zeros(nobs, NLEV);
  wgtTs = zeros(nobs, NLEV);
  wgtLs = zeros(nobs, NLEV);
  wgtRs = zeros(nobs, NLEV);

  wgtprvs = zeros(nobs, NLEV);
  wgtnows = zeros(nobs, NLEV);

% Determine interpolation points and weights
  for iob = iob0:iobF
    iLA  = find(obslat(iob) < grdlat, 1);
    iLO  = find(obslon(iob) < grdlon, 1);

    if (isempty(iLA)), iLA = NLAT; end
    if (isempty(iLO)), iLO = NLON; end

    wgtB = (grdlat(iLA) - obslat(iob))/(grdlat(iLA) - grdlat(iLA-1));
    wgtT = 1 - wgtB;
    wgtL = (grdlon(iLO) - obslon(iob))/(grdlon(iLO) - grdlon(iLO-1));
    wgtR = 1 - wgtL;

    wgtprv = 8*(dnobs(iob) - dnprv);
    wgtnow = 1 - wgtprv;

    iLOs(iob-iob0+1) = iLO;
    iLAs(iob-iob0+1) = iLA;

    wgtBs(iob-iob0+1,:)  = wgtB;
    wgtTs(iob-iob0+1,:)  = wgtT;
    wgtLs(iob-iob0+1,:)  = wgtL;
    wgtRs(iob-iob0+1,:)  = wgtR;

    wgtnows(iob-iob0+1,:) = wgtnow;
    wgtprvs(iob-iob0+1,:) = wgtprv;
  end

  iptTRs = iLOs   + (iLAs-1)*NLON;
  iptBRs = iLOs   + (iLAs-2)*NLON;
  iptTLs = iLOs-1 + (iLAs-1)*NLON;
  iptBLs = iLOs-1 + (iLAs-2)*NLON;

% Reshape arrays to make indexing easier
  rprv_dp = reshape(dpprv, NLON*NLAT, NLEV);
  rnow_dp = reshape(dpnow, NLON*NLAT, NLEV);

  rprv_qq = reshape(qqprv, NLON*NLAT, NLEV);
  rnow_qq = reshape(qqnow, NLON*NLAT, NLEV);

  rprv_co2 = reshape(co2prv, NLON*NLAT, NLEV);
  rnow_co2 = reshape(co2now, NLON*NLAT, NLEV);

% Interpolate
  oprv_dp = wgtTs.*wgtRs.*rprv_dp(iptTRs,:) + ...
            wgtBs.*wgtRs.*rprv_dp(iptBRs,:) + ...
            wgtTs.*wgtLs.*rprv_dp(iptTLs,:) + ...
            wgtBs.*wgtLs.*rprv_dp(iptBLs,:);
  onow_dp = wgtTs.*wgtRs.*rnow_dp(iptTRs,:) + ...
            wgtBs.*wgtRs.*rnow_dp(iptBRs,:) + ...
            wgtTs.*wgtLs.*rnow_dp(iptTLs,:) + ...
            wgtBs.*wgtLs.*rnow_dp(iptBLs,:);

  oprv_qq = wgtTs.*wgtRs.*rprv_qq(iptTRs,:)   + ...
            wgtBs.*wgtRs.*rprv_qq(iptBRs,:)   + ...
            wgtTs.*wgtLs.*rprv_qq(iptTLs,:)   + ...
            wgtBs.*wgtLs.*rprv_qq(iptBLs,:);
  onow_qq = wgtTs.*wgtRs.*rnow_qq(iptTRs,:)   + ...
            wgtBs.*wgtRs.*rnow_qq(iptBRs,:)   + ...
            wgtTs.*wgtLs.*rnow_qq(iptTLs,:)   + ...
            wgtBs.*wgtLs.*rnow_qq(iptBLs,:);

  oprv_co2 = wgtTs.*wgtRs.*rprv_co2(iptTRs,:) + ...
             wgtBs.*wgtRs.*rprv_co2(iptBRs,:) + ...
             wgtTs.*wgtLs.*rprv_co2(iptTLs,:) + ...
             wgtBs.*wgtLs.*rprv_co2(iptBLs,:);
  onow_co2 = wgtTs.*wgtRs.*rnow_co2(iptTRs,:) + ...
             wgtBs.*wgtRs.*rnow_co2(iptBRs,:) + ...
             wgtTs.*wgtLs.*rnow_co2(iptTLs,:) + ...
             wgtBs.*wgtLs.*rnow_co2(iptBLs,:);

  dpmod   = wgtprvs.*oprv_dp  + wgtnows.*onow_dp;
  qqmod   = wgtprvs.*oprv_qq  + wgtnows.*onow_qq;
  co2lmod = wgtprvs.*oprv_co2 + wgtnows.*onow_co2;

% Compute pressures at model interfaces (pemod) and layer centers (plmod)
  pemod = 0.01 + cumsum(dpmod,2);
  pemod = [0.01*ones(nobs,1), pemod];
  plmod = (   (pemod(:,2:end).^KAP1 - pemod(:,1:end-1).^KAP1) ...
           ./ (KAP1*(pemod(:,2:end) - pemod(:,1:end-1)))      ).^KAPR;

% Save some stuff for evaluation
  psmod(iob0:iobF)   = pemod(:,end);
  tco2mod(iob0:iobF) = sum(dpmod.*co2lmod,2);

% Compute pressures at obs edges (peobs) and layer centers (plobs)
  peobs = peavgs(iob0:iobF,:);
  plobs = 0.5*(peobs(:,2:end) + peobs(:,1:end-1));

% Compute tracer averages over obs levels
  modelavg = zeros(size(plobs));

  pegrd = zeros(size(peobs));
  for iob = 1:nobs
    pegrd(iob,:) = interp1([0, pemod(iob,:), 1e8], [1, [1:NLEV+1], NLEV+1], ...
                           peobs(iob,:), 'linear');
  end

% Compute pressure weight and tracer average at each obs layer center
  for kl = 1:size(plobs,2)
%   I think the protections should go here, not in top and ibot
%   because petop and pebot are used unprotected below;  nb. this is a
%   philosophical issue because they are protected above
    petop = min(pegrd(:,kl), pegrd(:,kl+1));
    pebot = max(pegrd(:,kl), pegrd(:,kl+1));

    itop = max(floor(petop),    1);
    ibot = min(floor(pebot), NLEV);

    avg = 0;
    dpr = 0;
    for ind = itop:ibot
      wgtz = 1;
      if (ind == itop), wgtz = wgtz - (petop - itop); end 
      if (ind == ibot), wgtz =         pebot - ibot;  end
      avg = avg + wgtz.*dpmod(:,ind).*co2lmod(:,ind);
      dpr = dpr + wgtz.*dpmod(:,ind);
    end
    avg = avg ./ dpr;

    modelavg(:,kl) = avg;
  end

  for iob = 1:nobs
    for kl = size(plobs,2):-1:2
      if (~isnan(modelavg(iob,kl))), break; end
    end
    modelavg(iob,kl+1:end) = modelavg(iob,kl);
  end

  tco2use(iob0:iobF) = sum(abs(diff(peobs,[],2)).*modelavg,2);

% Some gross error checks
  if (nnz(isnan(modelavg))  > 0), error('Computed value is NaN ...');      end
  if (nnz(modelavg <     0) > 0), error('Computed value is too low ...');  end
  if (nnz(modelavg > 10000) > 0), error('Computed value is too high ...'); end

% Compute model values on pressure edges
  modellev = zeros(size(peobs));
  modellev(:,1      ) =      modelavg(:,1);
  modellev(:,2:end-1) = 0.5*(modelavg(:,1:end-1) + modelavg(:,2:end));
  modellev(:,  end  ) =      modelavg(:,end);

% Choose between averaging kernels defined on edge or layer grids
  if (size(avgker,2) == size(peavgs,2)), modeluse = modellev; end
  if (size(avgker,2) ~= size(peavgs,2)), modeluse = modelavg; end

  xco2mod(iob0:iobF) = priorobs(iob0:iobF) ...
                     + sum(avgker(iob0:iobF,:).*(modeluse - priorlev(iob0:iobF,:)), 2);

% More gross error checks
  if (nnz(isnan(xco2mod(iob0:iobF))) > 0), error('Computed value is NaN ...'); end

% Save fields for next iteration
  dnprv  = dnnow;
  dpprv  = dpnow;
  qqprv  = qqnow;
  co2prv = co2now;
end
fprintf('\n');
