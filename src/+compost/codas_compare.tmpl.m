%CODAS_COMPARE  Compare model results to retrievals in CoDAS format

% Author(s):	Brad Weir <brad.weir@nasa.gov>
%
% Changelog:
% 2025-09-10	First version, based on lite_compare
%
% TODO:
%===============================================================================

% User supplied definitions
% ---
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
HEADID = '>>>HEADID<<<';
YEAR0  = >>>YEAR0<<<;
YEARF  = >>>YEARF<<<;

DNUM0 = datenum(YEAR0, 01, 01);				% Starting date of comparison (arbitrary)
DNUMF = datenum(YEARF, 12, 31);				% Ending   date of comparison (arbitrary)

% PRINT HEADER & SET UP ENVIRONMENT
%===============================================================================
sline = ['#', repmat('=', 1, 79)];
disp(sline);
disp('# CoDAS file comparison');
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
if isempty(VARQW) & ~ISDRY
    error('No water vapor variable present to dry trace gas.');
end

if 0 < numel(DIROBS) && DIROBS(end) == '/', DIROBS = DIROBS(1:end-1); end
if 0 < numel(DIRMOD) && DIRMOD(end) == '/', DIRMOD = DIRMOD(1:end-1); end

% Set up environment
% ---
atmomut.constants;

% Set up model time and space grid information
% ---
dnmod = [DNUM0:DSKIP:DNUMF]';

contents = dir([DIRMOD, '/', HDGAS, '*.nc4']);
if isempty(contents), exit; end
fgas = [contents(1).folder, '/', contents(1).name];

grdlat = ncread(fgas, 'lat');
grdlon = ncread(fgas, 'lon');
% Extend lat & lon for interp
grdlat = [grdlat(1) - (grdlat(2) - grdlat(1)); grdlat; ...
    grdlat(end) + (grdlat(end) - grdlat(end-1))];
grdlon = [grdlon(1) - (grdlon(2) - grdlon(1)); grdlon; ...
    grdlon(end) + (grdlon(end) - grdlon(end-1))];

NLAT = numel(grdlat);
NLON = numel(grdlon);
NLEV = 72;

disp(['Using model data in ', DIRMOD, ' ...']);
disp(['Using obs   data in ', DIROBS, ' ...']);
disp(' ');

% Hack for first iteration
dnprv  = Inf;
dpprv  = NaN*ones(NLON, NLAT, NLEV);
qqprv  = NaN*ones(NLON, NLAT, NLEV);
gasprv = NaN*ones(NLON, NLAT, NLEV);

for dnin = datenum(YEAR0,01,01):6/24:datenum(YEARF,12,31,18,00,00)
%   1. READ OBS DATA
%===============================================================================
    dvec = datevec(dnin);
    % Matches all netCDF files with or without minutes in time string
    flist1 = dir([DIROBS, '/Y', num2str(dvec(1)), '/*.', ...
        datestr(dnin, 'yyyymmdd_HH'), 'z.nc*']);
    flist2 = dir([DIROBS, '/Y', num2str(dvec(1)), '/*.', ...
        datestr(dnin, 'yyyymmdd_HHMM'), 'z.nc*']);
    flist = [flist1; flist2];

    if numel(flist) == 0, continue; end

    % Pick newest file if there are multiple matches
    fin = [flist(1).folder, '/', flist(1).name];
    fprintf(['Reading Y', num2str(dvec(1)), '/', flist(1).name, ' ...\n']);

    % Read necessary variables
    date  = ncread(fin, 'date');
    time  = ncread(fin, 'time');
    dnobs = datenum([num2str(date,'%06u'), num2str(time,'%06u')], ...
         'yyyymmddHHMMSS');

    lat = ncread(fin, 'lat');
    lon = ncread(fin, 'lon');
    obs = ncread(fin, 'obs');
    avgker = ncread(fin, 'avgker')';
    peavg  = ncread(fin, 'peavg')';

    isbad    = zeros(size(obs));
    uncert   = zeros(size(obs));
    priorobs = zeros(size(obs));
    priorpro = zeros(size(obs));

    % Read optional variables
    try, isbad    = ncread(fin, 'isbad');     end
    try, uncert   = ncread(fin, 'uncert');    end
    try, priorobs = ncread(fin, 'priorobs');  end
    try, priorpro = ncread(fin, 'priorpro')'; end

    gesobs = NaN*ones(size(obs));
    plavg = 0.5*(peavg(:,2:end) + peavg(:,1:end-1));

%   2. READ MODEL DATA
%===============================================================================
    for nn = 1:2
        dnnow = dnin + (nn - 1)*3/24;
        date = datestr(dnnow, 'yyyymmdd');
        time = datestr(dnnow, 'HH');

        % A bit of complexity to deal with HH and HHMM formats
        fgas = [DIRMOD, '/', HDGAS, date, datestr(dnnow,TTGAS), '.nc4'];
        fmet = [DIRMOD, '/', HDMET, date, datestr(dnnow,TTMET), '.nc4'];

        % A. Read model output, skipping missing/broken files
        % ---
        try
            gasnow = SCLMOD*ncread(fgas, VARMOD);
            dpnow  = 1e-2*atmomut.getdp(ncread(fmet, VARPS), NLEV);
            qqnow  = zeros(size(gasnow));
            if ~isempty(VARQW), qqnow = ncread(fmet, VARQW); end

            fprintf(['Reading model on ', date, ' at ', time, 'Z', ' ...\n']);
        catch ME
            fprintf([ME.message, '\n']);
            continue;
        end

        % Extend lat for constant interp
        gasnow = cat(2, gasnow(:,1,:), gasnow, gasnow(:,end,:));
        dpnow  = cat(2,  dpnow(:,1,:),  dpnow,  dpnow(:,end,:));
        qqnow  = cat(2,  qqnow(:,1,:),  qqnow,  qqnow(:,end,:));
        % Extend lon for periodic interp
        gasnow = cat(1, gasnow(end,:,:), gasnow, gasnow(1,:,:));
        dpnow  = cat(1,  dpnow(end,:,:),  dpnow,  dpnow(1,:,:));
        qqnow  = cat(1,  qqnow(end,:,:),  qqnow,  qqnow(1,:,:));

        % Convert to dry-air mole fractions
        % ---
        % ***           All GEOS-5 trace gas molar mixing ratios are            ***
        % ***                   "virtual molar mixing ratios"                   ***
        % They are treated as rescalings of the mass mixing ratio by the
        % DRY-AIR molecular masses, i.e.:
        %     1. X_dry = X_total / (1 - qq),         regardless of mass/molar
        %     2. X_mol = X_mass * MAIR_DRY/MOLMASS   regardless of dry/total
        if ~ISDRY, gasnow = gasnow./(1 - qqnow); end

%       3. COMPUTE MODEL vs OBS COMPARISON
%===============================================================================
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
        for iob = 1:nobs
            iLA  = find(lat(iob0+iob-1) < grdlat, 1);
            iLO  = find(lon(iob0+iob-1) < grdlon, 1);

            if isempty(iLA), iLA = NLAT; end
            if isempty(iLO), iLO = NLON; end

            wgtB = (grdlat(iLA) - lat(iob0+iob-1))/(grdlat(iLA) - grdlat(iLA-1));
            wgtT = 1 - wgtB;
            wgtL = (grdlon(iLO) - lon(iob0+iob-1))/(grdlon(iLO) - grdlon(iLO-1));
            wgtR = 1 - wgtL;

            wgtprv = 8*(dnobs(iob0+iob-1) - dnprv);
            wgtnow = 1 - wgtprv;

            iLOs(iob) = iLO;
            iLAs(iob) = iLA;

            wgtBs(iob,:)  = wgtB;
            wgtTs(iob,:)  = wgtT;
            wgtLs(iob,:)  = wgtL;
            wgtRs(iob,:)  = wgtR;

            wgtnows(iob,:) = wgtnow;
            wgtprvs(iob,:) = wgtprv;
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

        rprv_gas = reshape(gasprv, NLON*NLAT, NLEV);
        rnow_gas = reshape(gasnow, NLON*NLAT, NLEV);

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

        oprv_gas = wgtTs.*wgtRs.*rprv_gas(iptTRs,:) + ...
                   wgtBs.*wgtRs.*rprv_gas(iptBRs,:) + ...
                   wgtTs.*wgtLs.*rprv_gas(iptTLs,:) + ...
                   wgtBs.*wgtLs.*rprv_gas(iptBLs,:);
        onow_gas = wgtTs.*wgtRs.*rnow_gas(iptTRs,:) + ...
                   wgtBs.*wgtRs.*rnow_gas(iptBRs,:) + ...
                   wgtTs.*wgtLs.*rnow_gas(iptTLs,:) + ...
                   wgtBs.*wgtLs.*rnow_gas(iptBLs,:);

        dpmod   = wgtprvs.*oprv_dp  + wgtnows.*onow_dp;
        qqmod   = wgtprvs.*oprv_qq  + wgtnows.*onow_qq;
        geslmod = wgtprvs.*oprv_gas + wgtnows.*onow_gas;

        % Compute pressures at model interfaces (pemod) and decimal grid
        % index (pegrd)
        pemod = 0.01 + cumsum(dpmod,2);
        pemod = [0.01*ones(nobs,1), pemod];

        % Compute pressure weight and tracer average at each obs layer center
        geslavg = zeros(size(plavg(iob0:iobF,:)));
        for iob = 1:nobs
            geslavg(iob,:) = atmomut.zavg(pemod(iob,:), geslmod(iob,:), ...
                peavg(iob0+iob-1,:));
        end

        % Compute model values on pressure edges
        geseavg = zeros(size(peavg(iob0:iobF,:)));
        geseavg(:,1      ) =      geslavg(:,1);
        geseavg(:,2:end-1) = 0.5*(geslavg(:,1:end-1) + geslavg(:,2:end));
        geseavg(:,  end  ) =      geslavg(:,end);

        % Choose between averaging kernels defined on edge or layer grids
        if size(peavg,2) == size(avgker,2) + 1
            gesuse = geslavg;
        else
            gesuse = geseavg;
        end
        gesobs(iob0:iobF) = priorobs(iob0:iobF) + ...
            sum(avgker(iob0:iobF,:).*(gesuse - priorpro(iob0:iobF,:)), 2);

        % Save fields for next iteration
        dnprv  = dnnow;
        dpprv  = dpnow;
        qqprv  = qqnow;
        gasprv = gasnow;
    end

    if nnz(~isnan(gesobs)) == 0, continue; end

    % Write
    dirout = [HEADID, '/Y', num2str(dvec(1))];
    if ~isfolder(dirout)
        [status, result] = system(['mkdir -p ', dirout]);
    end

    fout = [dirout, '/', HEADID, '.', datestr(dnin,'yyyymmdd_HH'), 'z.nc'];
    fprintf(['Writing ', fout, ' ...\n']);

    if ~isfile(fout)
        nccreate(fout, 'lat',   'dimensions',{'nsound',numel(obs)}, ...
            'format','netcdf4', 'deflate',9, 'shuffle',true);
        nccreate(fout, 'lon',   'dimensions',{'nsound',numel(obs)}, ...
            'format','netcdf4', 'deflate',9, 'shuffle',true);
        nccreate(fout, 'obs',   'dimensions',{'nsound',numel(obs)}, ...
            'format','netcdf4', 'deflate',9, 'shuffle',true);
        nccreate(fout, 'mod',   'dimensions',{'nsound',numel(obs)}, ...
            'format','netcdf4', 'deflate',9, 'shuffle',true);
        nccreate(fout, 'isbad', 'dimensions',{'nsound',numel(obs)}, ...
            'format','netcdf4', 'deflate',9, 'shuffle',true);
        gesout = gesobs;
    else
        % Retain any non-NaN values that are NaNs (for time boundaries)
        gesout = ncread(fout, 'mod');
        inds = find(~isnan(gesobs));
        gesout(inds) = gesobs(inds);
    end

    ncwrite(fout, 'lat',   lat);
    ncwrite(fout, 'lon',   lon);
    ncwrite(fout, 'obs',   obs);
    ncwrite(fout, 'mod',   gesout);
    ncwrite(fout, 'isbad', isbad);
end
fprintf('\n');
