%OBSPACK_STATION_COMPARE  Compare model results to NOAA ObsPack station data
%
%   This comparison is fast, but only updates the point in space the model is
%   evaluated at every 3 hours.  This is acceptable for stationary sites, but
%   is UNACCEPTABLE for MOBILE (viz. AIRCRAFT, SHIPBOARD, and AIRCORE data).
%   For those data types, use the obspack_mobile_compare utility.

% Author(s):	Brad Weir <brad.weir@nasa.gov>
%
% Changelog:
% 2019-03-18	New version
% 2019-04-18	Changed qcflags to all x's if variable doesn't exist
% 2019-11-13	Many changes, mostly zero-diff
% 2023-02-15	Adding capability to read packed, daily files (e.g., MERRA-2 GMI)
%		Includes file format change, update to shell script
% 2025-12-03	Read NLEV from collection
% 2025-12-03	Cosmetic changes to follow new "standards" (4-space indents, no
%		parens around ifs, etc.); should be zero-diff
% 2025-12-03	Add support for edge collections
% 2025-12-04	Change TTGAS and TTMET to be entire date-time tag
%
% TODO:
% * Pack variables into structure instead of cell arrays
% * Improve qcflag handling (non-trivial, varies across datasets)
% * Keep a record somehow of obs lons that get changed or tossed
% * Maybe save model elevations?

% "Virtual-air":
% ---
% GEOS holds tracers as total-air mass mixing ratios multiplied by a constant,
% which is usually the ratio of molar masses of the gas and dry-air.  The
% result is something that looks like a total-air molar mixing ratio, but
% isn't.  We call such values "virtual-air" molar mixing ratios.
%
% Basic naming conventions:
% ---
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
%===============================================================================

% User supplied definitions
% ---
ISDRY   = >>>ISDRY<<<;					% Are model mixing ratios dry air or total?
SCLOBS  = >>>SCLOBS<<<;					% Scaling from   obs units to output units
SCLMOD  = >>>SCLMOD<<<;					% Scaling from model units to output units
VAROBS  = >>>VAROBS<<<;					% Obs name for trace gas
VARMOD  = >>>VARMOD<<<;					% Model name for trace gas

VARQW   = >>>VARQW<<<;					% Model name for total water mixing ratio
VARPS   = >>>VARPS<<<;					% Model name for surface pressure
% Disabled for now
% (Is there a way to have these default to empty so that I don't have to make
% every RC file define them?)
VARPL   = ''; %VARPL   = >>>VARPL<<<;					% Model name for pressure levels
VARPHIS = >>>VARPHIS<<<;				% Model name for surface geopotential
VARZL   = >>>VARZL<<<;					% Model name for geopotential height

DIROBS  = >>>DIROBS<<<;					% ObsPack directory
DIRMOD  = >>>DIRMOD<<<;					% Directory holding model files
HDMET   = >>>HDMET<<<;					% Head of model meteo files
HDGAS   = >>>HDGAS<<<;					% Head of model trace gas files
TTMET   = >>>TTMET<<<;					% Time format of model meteo files
TTGAS   = >>>TTGAS<<<;					% Time format of model trace gas files
DSKIP   = >>>DSKIP<<<;					% Fraction of day between files
YEAR0   = >>>YEAR0<<<;					% First year
YEARF   = >>>YEARF<<<;					% Final year

DNUM0   = datenum(YEAR0, 01, 01);			% Starting date of comparison (arbitrary)
DNUMF   = datenum(YEARF, 12, 31);			% Ending   date of comparison (arbitrary)

% PRINT HEADER & SET UP ENVIRONMENT
%===============================================================================
sline = ['#', repmat('=', 1, 79)];
disp(sline);
disp('# NOAA ObsPack comparison (STATION version)');
disp('#');
if ISDRY
    disp('# *** Treating model values as     DRY-air MOLE fractions ...');
    disp('#                                  ---     ----              ');
else
    disp('# *** Treating model values as VIRTUAL-air MOLE fractions ...');
    disp('#                              -------     ----              ');
end
disp(sline);
disp(' ');

% Some error checking
% ---
if isempty(VARQW) & ~ISDRY
    error('No water vapor variable present to dry trace gas.');
end

if isempty(VARZL) & isempty(VARPHIS)
    error('Need either level heights (VARZL) or surface geopotential (VARPHIS).');
end

% Set up environment
% ---
atmomut.constants;

% 1. READ OBS DATA
%===============================================================================
disp('--- Observational data ---');
disp(['Using data in ', DIROBS, ' ...']);
disp(' ');

contents  = dir([DIROBS, '/', VAROBS, '_*.nc']);
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
    fobs = [DIROBS, '/', cell_fobs{ic}];
    fprintf(['Reading ', cell_fobs{ic}, ' ...\n']);

    dvec  = ncread(fobs, 'time_components');
    % Hack for empty files (private data)
    if isempty(dvec), continue; end

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
        if size(qcf_in) == size(qcflags), qcflags = qcf_in; end
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
%===============================================================================
disp('--- Model comparison ---');

% Get time and grid variables from gas file
% ---
dnmod = [DNUM0:DSKIP:DNUMF]';

contents = dir([DIRMOD, '/', HDGAS, '*.nc4']);
if isempty(contents), exit; end
fgas = [contents(1).folder, '/', contents(1).name];

grdlat = ncread(fgas, 'lat');
grdlon = ncread(fgas, 'lon');
% Extend lat & lon for interp (constant in lat, periodic in lon)
grdlat = [grdlat(1) - (grdlat(2) - grdlat(1)); grdlat; ...
    grdlat(end) + (grdlat(end) - grdlat(end-1))];
grdlon = [grdlon(1) - (grdlon(2) - grdlon(1)); grdlon; ...
    grdlon(end) + (grdlon(end) - grdlon(end-1))];

NLEV = numel(ncread(fgas, 'lev'));

% Get geopotential height from met file (if requested)
% ---
if ~isempty(VARPHIS)
    contents = dir([DIRMOD, '/', HDMET, '*.nc4']);
    if isempty(contents), exit; end
    fmet = [contents(1).folder, '/', contents(1).name];

    grdzsin = 1/GRAV * ncread(fmet, VARPHIS);
    % Hack to deal with packed daily files
    if size(grdzsin,3) == 1
        grdzs = grdzsin;
    else
        grdzs = grdzsin(:,:,1);
    end
    % Extend lat & lon for interp
    grdzs = cat(2, grdzs(  :,1,:), grdzs, grdzs(:,end,:));
    grdzs = cat(1, grdzs(end,:,:), grdzs, grdzs(1,  :,:));
end

% Allocate output arrays and fill w/ nans
% ---
cell_gasmod = cell(NSITES, 1);
cell_qqmod  = cell(NSITES, 1);
cell_prsmod = cell(NSITES, 1);

for ic = 1:NSITES
    cell_gasmod{ic} = NaN*ones(size(dnmod));
    cell_qqmod{ic}  = NaN*ones(size(dnmod));
    cell_prsmod{ic} = NaN*ones(size(dnmod));
end

lread = 1;				% Switch to minimize reads
for it = 1:numel(dnmod)
    dnnow = dnmod(it);
    date  = datestr(dnnow, 'yyyymmdd');

    fgas = [DIRMOD, '/', HDGAS, datestr(dnnow,TTGAS), '.nc4'];
    fmet = [DIRMOD, '/', HDMET, datestr(dnnow,TTMET), '.nc4'];

    % A. Read model output, skipping missing/broken files
    % ---
    if lread == 1, try
        gasin = SCLMOD*ncread(fgas, VARMOD);

        qqin  = NaN*ones(size(gasin));
        plin  = NaN*ones(size(gasin));
        zlin  = NaN*ones(size(gasin));

        if ~isempty(VARQW), qqin = ncread(fmet, VARQW); end
        if ~isempty(VARPL)
            plin = ncread(fmet, VARPL);
        elseif ~isempty(VARPS)
            psin = ncread(fmet, VARPS);
            % Hack to support edge collections
            if NLEV == 73
                [dpin, pein, plin] = atmomut.getdp(psin, NLEV-1);
                plin = pein;
            else
                [dpin, pein, plin] = atmomut.getdp(psin, NLEV);
            end
        end
        if ~isempty(VARZL)
            zlin = ncread(fmet, VARZL);
        end

        % Print daily message if read succeeds
        if floor(dnnow) == dnnow
            fprintf(['Computing values for ', date, ' ...\n']);
        end
    catch ME
        continue;
    end, end

    % Hack to deal with packed daily files
    if size(gasin,4) == 1
        gas = gasin;
        qq  = qqin;
        pl  = plin;
        zl  = zlin;
    else
        nn  = round((dnnow - floor(dnnow))/DSKIP) + 1;
        gas = gasin(:,:,:,nn);
        qq  =  qqin(:,:,:,nn);
        pl  =  plin(:,:,:,nn);
        zl  =  zlin(:,:,:,nn);

        % Try to minimize reads
        lread = 0;
        if floor(dnnow + DSKIP) == dnnow + DSKIP, lread = 1; end
    end
    % Convert from Pa to hPa
    pl = 1e-2*pl;

    % Extend lat & lon for interp
    gas = cat(2, gas(  :,1,:), gas, gas(:,end,:));
    qq  = cat(2,  qq(  :,1,:),  qq,  qq(:,end,:));
    pl  = cat(2,  pl(  :,1,:),  pl,  pl(:,end,:));
    zl  = cat(2,  zl(  :,1,:),  zl,  zl(:,end,:));
    gas = cat(1, gas(end,:,:), gas, gas(1,  :,:));
    qq  = cat(1,  qq(end,:,:),  qq,  qq(1,  :,:));
    pl  = cat(1,  pl(end,:,:),  pl,  pl(1,  :,:));
    zl  = cat(1,  zl(end,:,:),  zl,  zl(1,  :,:));

    % B. Interpolate to time and space of observations
    % ---
    for ic = 1:NSITES
        dnobs = cell_dnobs{ic};

        % Skip sites that have no data for this time
        if numel(dnobs) == 0 | dnnow+DSKIP < dnobs(1) | ...
            dnobs(end) < dnnow-DSKIP
            continue;
        end

        % Determine spatial location for interpolation
        % (This could use some improvement, but this program is only
        % appropriate when the obs location doesn't change much)
        [junk, iu] = unique(dnobs, 'stable');

        if 1 < numel(iu)
            alt = interp1(dnobs(iu), cell_alt{ic}(iu), dnnow, 'linear', 'extrap');
            lat = interp1(dnobs(iu), cell_lat{ic}(iu), dnnow, 'linear', 'extrap');
            lon = interp1(dnobs(iu), cell_lon{ic}(iu), dnnow, 'linear', 'extrap');
        else
            alt = cell_alt{ic}(iu);
            lat = cell_lat{ic}(iu);
            lon = cell_lon{ic}(iu);
        end

        % Determine spatial interpolation points and weights
        iLA  = find(lat < grdlat, 1);
        iLO  = find(lon < grdlon, 1);

        wgtB = (grdlat(iLA) - lat)/(grdlat(iLA) - grdlat(iLA-1));
        wgtT = 1 - wgtB;
        wgtL = (grdlon(iLO) - lon)/(grdlon(iLO) - grdlon(iLO-1));
        wgtR = 1 - wgtL;

        % Compute model values on native levels
        gasnat = wgtT*wgtR*squeeze(gas(iLO,iLA,:))   + ...
                 wgtB*wgtR*squeeze(gas(iLO,iLA-1,:)) + ...
                 wgtT*wgtL*squeeze(gas(iLO-1,iLA,:)) + ...
                 wgtB*wgtL*squeeze(gas(iLO-1,iLA-1,:));
        qqnat  = wgtT*wgtR*squeeze(qq(iLO,iLA,:))   + ...
                 wgtB*wgtR*squeeze(qq(iLO,iLA-1,:)) + ...
                 wgtT*wgtL*squeeze(qq(iLO-1,iLA,:)) + ...
                 wgtB*wgtL*squeeze(qq(iLO-1,iLA-1,:));
        plnat  = wgtT*wgtR*squeeze(pl(iLO,iLA,:))   + ...
                 wgtB*wgtR*squeeze(pl(iLO,iLA-1,:)) + ...
                 wgtT*wgtL*squeeze(pl(iLO-1,iLA,:)) + ...
                 wgtB*wgtL*squeeze(pl(iLO-1,iLA-1,:));
        if ~isempty(VARZL)
            zlnat = wgtT*wgtR*squeeze(zl(iLO,iLA,:))   + ...
                    wgtB*wgtR*squeeze(zl(iLO,iLA-1,:)) + ...
                    wgtT*wgtL*squeeze(zl(iLO-1,iLA,:)) + ...
                    wgtB*wgtL*squeeze(zl(iLO-1,iLA-1,:));
        else
            zsnat = wgtT*wgtR*grdzs(iLO,iLA)   + ...
                    wgtB*wgtR*grdzs(iLO,iLA-1) + ...
                    wgtT*wgtL*grdzs(iLO-1,iLA) + ...
                    wgtB*wgtL*grdzs(iLO-1,iLA-1);
            zlnat = zsnat + BCON*(1 - (plnat/plnat(end)).^BPOW);
        end

        % Convert to dry-air mole fractions
        % ---
        % ***           All GEOS trace gas molar mixing ratios are            ***
        % ***           actually "virtual molar mixing ratios"                ***
        % Thus, they are treated as rescalings of the mass mixing ratio by the
        % DRY-AIR molecular masses, i.e.:
        %     1. XX_dry = XX_total / (1 - qq),         regardless of mass/molar
        %     2. XX_mol = XX_mass * MAIR_DRY/MOLMASS   regardless of dry/total
        if ~ISDRY, gasnat = gasnat./(1 - qqnat); end

        % Interpolate to vertical location of obs
        gasmod = interp1(zlnat, gasnat, alt, 'linear', gasnat(end));
        qqmod  = interp1(zlnat, qqnat,  alt, 'linear', qqnat(end));
        prsmod = interp1(zlnat, plnat,  alt, 'linear', plnat(end));

        cell_prsmod{ic}(it) = prsmod;
        cell_qqmod{ic}(it)  = qqmod;
        cell_gasmod{ic}(it) = gasmod;
    end
end
fprintf('\n');

% Keep file sizes small
clearvars -except NSITES dnmod cell_*
