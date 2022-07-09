% Bug testing
fb = '../covid_free/covid_free__obspack_co2_mip_v3.mat';
fa = 'm2cc_ana__obspack_co2_mip_v3.mat';
expid = 'm2cc';
units = 'ppmv';
%fb = input('Backgound filename: ', 's');
%fa = input('Analysis filename: ', 's');
%expid = input('Experiment id: ', 's');
%units = input('Units: ', 's');

disp(['Comparing ', fb, ' to ']);
disp(['          ', fa, ' ...']);

load(fb);
va = load(fa);

cell_gasana = va.cell_gasmod;
clear va;

% Get gas name
fobs = cell_fobs{1};
ii   = strfind(fobs, '_');
gas  = fobs(1:ii-1);

addpath('/discover/nobackup/bweir/matlab');
cmapd = flipud(brewermap(128,'BrBg'));
style = hgexport('readstyle', 'hires');

COLOR1 = [0   135 255]/255;
COLOR2 = [255 2   51 ]/255;
COLOR3 = [153 102 255]/255;

nn = 0;

% ATom/CONTRAIL/ORCAS
% -------------------
nn = nn + 1;
xx(nn).tag   = [gas, '_atom1'];
xx(nn).title = 'JJA: ATom 1/CON';
xx(nn).names = {[gas, '_tom_aircraft-insitu_1_allvalid.nc'], ...
                [gas, '_con_aircraft-flask_42_allvalid'], ...
                [gas, '_con_aircraft-insitu_42_allvalid']};
xx(nn).DLIMS = [datenum(2016,07,01), datenum(2016,09,01)];
xx(nn).XLIMS = [-90, 90];
xx(nn).YLIMS = [0, 15];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));


nn = nn + 1;
xx(nn).tag   = [gas, '_atom2'];
xx(nn).title = 'DJF: ATom 2/CON/ORCAS';
xx(nn).names = {[gas, '_tom_aircraft-insitu_1_allvalid.nc'], ...
                [gas, '_orc_aircraft-insitu_3_allvalid-merge10.nc'], ...
                [gas, '_con_aircraft-flask_42_allvalid'], ...
                [gas, '_con_aircraft-insitu_42_allvalid']};
xx(nn).DLIMS = [datenum(2016,01,01), datenum(2016,04,01); ...
                datenum(2017,01,01), datenum(2017,03,01)];
xx(nn).XLIMS = [-90, 90];
xx(nn).YLIMS = [0, 15];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));

nn = nn + 1;
xx(nn).tag   = [gas, '_atom3'];
xx(nn).title = 'SON: ATom 3/CON';
xx(nn).names = {[gas, '_tom_aircraft-insitu_1_allvalid.nc'], ...
                [gas, '_con_aircraft-flask_42_allvalid'], ...
                [gas, '_con_aircraft-insitu_42_allvalid']};
xx(nn).DLIMS = [datenum(2017,09,01), datenum(2017,11,01)];
xx(nn).XLIMS = [-90, 90];
xx(nn).YLIMS = [0, 15];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));

nn = nn + 1;
xx(nn).tag   = [gas, '_atom4'];
xx(nn).title = 'MAM: ATom 4/CON';
xx(nn).names = {[gas, '_tom_aircraft-insitu_1_allvalid.nc'], ...
                [gas, '_con_aircraft-flask_42_allvalid'], ...
                [gas, '_con_aircraft-insitu_42_allvalid']};
% Hack to include CONTRAIL (isnt in Obspack for 2018 onward)
xx(nn).DLIMS = [datenum(2017,04,01), datenum(2017,06,01); ...
                datenum(2018,04,01), datenum(2018,06,01)];
%xx(nn).DLIMS = [datenum(2018,04,01), datenum(2018,06,01)];
xx(nn).XLIMS = [-90, 90];
xx(nn).YLIMS = [0, 15];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));

% ACT Clear-sky & OCO-2 underflights
% ----------------------------------
nn = nn + 1;
xx(nn).tag   = [gas, '_act2016summer'];
xx(nn).title = 'July-Aug ''16: ACT-America';
xx(nn).names = {[gas, '_act_aircraft-insitu_428_allvalid-b200.nc'], ...
                [gas, '_act_aircraft-insitu_428_allvalid-c130.nc'], ...;
                [gas, '_act_aircraft-pfp_1_allvalid-b200.nc'], ...
                [gas, '_act_aircraft-pfp_1_allvalid-c130.nc']};
xx(nn).DLIMS = [datenum(2016,07,27), datenum(2016,07,28); ... % OCO-2
                datenum(2016,08,05), datenum(2016,08,06); ... % OCO-2
                datenum(2016,08,13), datenum(2016,08,14); ... % Clear
                datenum(2016,08,14), datenum(2016,08,15); ... % Clear
                datenum(2016,08,22), datenum(2016,08,23); ... % Clear
                datenum(2016,08,28), datenum(2016,08,29); ... % Clear
                datenum(2016,08,27), datenum(2016,08,28)];    % OCO-2
xx(nn).XLIMS = [30, 50];
xx(nn).YLIMS = [0, 10];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));

nn = nn + 1;
xx(nn).tag   = [gas, '_act2017winter'];
xx(nn).title = 'Feb-Mar ''17: ACT-America';
xx(nn).names = {[gas, '_act_aircraft-insitu_428_allvalid-b200.nc'], ...
                [gas, '_act_aircraft-insitu_428_allvalid-c130.nc'], ...;
                [gas, '_act_aircraft-pfp_1_allvalid-b200.nc'], ...
                [gas, '_act_aircraft-pfp_1_allvalid-c130.nc']};
xx(nn).DLIMS = [datenum(2017,02,13), datenum(2017,02,14); ... % OCO-2
                datenum(2017,02,15), datenum(2017,02,16); ... % OCO-2
                datenum(2017,03,08), datenum(2017,03,09)];    % OCO-2
xx(nn).XLIMS = [30, 50];
xx(nn).YLIMS = [0, 10];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));

nn = nn + 1;
xx(nn).tag   = [gas, '_act2017fall'];
xx(nn).title = 'Oct-Nov ''17: ACT-America';
xx(nn).names = {[gas, '_act_aircraft-insitu_428_allvalid-b200.nc'], ...
                [gas, '_act_aircraft-insitu_428_allvalid-c130.nc'], ...;
                [gas, '_act_aircraft-pfp_1_allvalid-b200.nc'], ...
                [gas, '_act_aircraft-pfp_1_allvalid-c130.nc']};
xx(nn).DLIMS = [datenum(2017,10,03), datenum(2017,10,04); ... % Clear
                datenum(2017,10,22), datenum(2017,10,23); ... % OCO-2
                datenum(2017,10,27), datenum(2017,10,28); ... % OCO-2
                datenum(2017,11,06), datenum(2017,11,07); ... % OCO-2
                datenum(2017,11,09), datenum(2017,11,10)];    % OCO-2
xx(nn).XLIMS = [30, 50];
xx(nn).YLIMS = [0, 10];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));

nn = nn + 1;
xx(nn).tag   = [gas, '_act2018spring'];
xx(nn).title = 'Apr-May ''18: ACT-America';
xx(nn).names = {[gas, '_act_aircraft-insitu_428_allvalid-b200.nc'], ...
                [gas, '_act_aircraft-insitu_428_allvalid-c130.nc'], ...;
                [gas, '_act_aircraft-pfp_1_allvalid-b200.nc'], ...
                [gas, '_act_aircraft-pfp_1_allvalid-c130.nc']};
xx(nn).DLIMS = [datenum(2018,04,25), datenum(2018,04,26); ... % OCO-2
                datenum(2018,04,27), datenum(2018,04,28); ... % OCO-2
                datenum(2018,04,29), datenum(2018,04,30); ... % OCO-2
                datenum(2018,05,08), datenum(2018,05,09)];    % Clear
xx(nn).XLIMS = [30, 50];
xx(nn).YLIMS = [0, 10];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));

% ABoVE (Arctic-CAP)
% ------------------
nn = nn + 1;
xx(nn).tag   = [gas, '_above2017spring'];
xx(nn).title = 'Apr-May ''17: Arctic-CAP';
xx(nn).names = {[gas, '_above_aircraft-insitu_1_allvalid.nc'], ...
                [gas, '_above_aircraft-pfp_1_allvalid.nc']};
xx(nn).DLIMS = [datenum(2017,04,01), datenum(2017,06,01)];
xx(nn).XLIMS = [55, 75];
xx(nn).YLIMS = [0, 6];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));

nn = nn + 1;
xx(nn).tag   = [gas, '_above2017summer'];
xx(nn).title = 'Jun-Aug ''17: Arctic-CAP';
xx(nn).names = {[gas, '_above_aircraft-insitu_1_allvalid.nc'], ...
                [gas, '_above_aircraft-pfp_1_allvalid.nc']};
xx(nn).DLIMS = [datenum(2017,06,01), datenum(2017,09,01)];
xx(nn).XLIMS = [55, 75];
xx(nn).YLIMS = [0, 6];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));

nn = nn + 1;
xx(nn).tag   = [gas, '_above2017fall'];
xx(nn).title = 'Sep-Nov ''17: Arctic-CAP';
xx(nn).names = {[gas, '_above_aircraft-insitu_1_allvalid.nc'], ...
                [gas, '_above_aircraft-pfp_1_allvalid.nc']};
xx(nn).DLIMS = [datenum(2017,09,01), datenum(2017,12,01)];
xx(nn).XLIMS = [55, 75];
xx(nn).YLIMS = [0, 6];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));

% ABoVE (ASCENDS)
% ---------------
nn = nn + 1;
xx(nn).tag   = [gas, '_ascends2017summer'];
xx(nn).title = 'Jun-Aug ''17: ASCENDS';
xx(nn).names = {[gas, '_gsfc_aircraft-insitu_430_allvalid.nc']};
xx(nn).DLIMS = [datenum(2017,06,01), datenum(2017,09,01)];
xx(nn).XLIMS = [55, 75];
xx(nn).YLIMS = [0, 12];
xx(nn).CLIMS = [-3, 3];
xx(nn).OKFUN = @(dd,lat,lon) ones(size(dd));

for nn = 1:numel(xx)
  alts = [];
  lats = [];
  omf1s = [];
  omf2s = [];

  for jj = 1:numel(xx(nn).names)
%   Find site number
    sname = xx(nn).names{jj};
    for ic = 1:NSITES
      kk = strfind(cell_fobs{ic}, sname);
      if (~isempty(kk)), break; end
    end

    if (isempty(kk))
      error(['Could not find matching site name: ', sname, ' ...']);
    end

    dnobs   = cell_dnobs{ic};
    lat     = cell_lat{ic};
    lon     = cell_lon{ic};
    alt     = cell_alt{ic};
    qcflags = cell_qcflags{ic};
    gasobs  = cell_gasobs{ic};
    gasmod  = cell_gasmod{ic};
    gasana  = cell_gasana{ic};

    iokd = find(qcflags(1,:) == '.' & qcflags(2,:) == '.');
    iokx = find(qcflags(1,:) == 'x' & qcflags(2,:) == 'x');
    iok1 = union(iokd, iokx);

    iok2 = find(~isnan(gasmod));
    iok3 = find(~isnan(gasana));
    iok4 = find(xx(nn).OKFUN(dnobs, lat, lon));

    iok  = intersect(intersect(intersect(iok1, iok2), iok3), iok4);

    for ii = 1:size(xx(nn).DLIMS,1)
      DNUM0 = xx(nn).DLIMS(ii,1);
      DNUMF = xx(nn).DLIMS(ii,2);

      iuse = find(DNUM0 <= dnobs & dnobs < DNUMF);
      iuse = intersect(iuse, iok);

      dmin = floor((dnobs(iuse) - datenum(1900,01,01))*24.*60.);
      [cc, id, it] = unique(dmin, 'stable');

%     Apply additive mass fix to simulations
%     --------------------------------------
      fix1  = 0.063*(dnobs - datenum(2015,01,01))/365.25;
      fix2  = 0*dnobs;

%     alts  = [alts; alt(iuse)];
%     lats  = [lats; lat(iuse)];
%     omf1s = [omf1s; gasobs(iuse) - (gasmod(iuse) - fix1(iuse))];
%     omf2s = [omf2s; gasobs(iuse) - (gasana(iuse) - fix2(iuse))];

      alta  = accumarray(it, alt(iuse), [], @mean);
      lata  = accumarray(it, lat(iuse), [], @mean);
      omf1a = accumarray(it, gasobs(iuse) - (gasmod(iuse) - fix1(iuse)), [], @mean);
      omf2a = accumarray(it, gasobs(iuse) - (gasana(iuse) - fix2(iuse)), [], @mean);

      alts  = [alts;  alt(iuse)];
      lats  = [lats;  lat(iuse)];
      omf1s = [omf1s; omf1a(it)];
      omf2s = [omf2s; omf2a(it)];
    end
  end

  leg1 = ['(', num2str(     mean(omf1s),    '%.2f'), ', ', ...
               num2str(sqrt(mean(omf1s.^2)),'%.2f'), ')'];
  leg2 = ['(', num2str(     mean(omf2s),    '%.2f'), ', ', ...
               num2str(sqrt(mean(omf2s.^2)),'%.2f'), ')'];

% "Simulation" pane
% -----------------
  figure(1); clf;
  scatter(lats, alts*1e-3, 8, omf1s, 's', 'filled');

  xlim(xx(nn).XLIMS);
  ylim(xx(nn).YLIMS);
  caxis(xx(nn).CLIMS);
  colormap(cmapd);
  grid on;
  set(gca, 'fontsize', 20);

% NB: leg1 and 'SIM'
  legend(leg1, 'location', 'northeast');
  text(0.03, 0.94, 'SIM', 'units', 'normalized', ...
       'fontsize', 30, 'fontweight', 'bold');

% Only for top pane
  text(0.03, 1.08, xx(nn).title, 'units', 'normalized', ...
       'fontsize', 30, 'fontweight', 'bold');
  set(gca, 'xticklabel', ''); ylabel('Altitude (km)');
  colorbar horz;

  set(gca, 'Position', [0.1100, 0.4373-0.28, 0.7750+0.09, 0.4700+0.28]);
  hgexport(gcf, ['figs/m2cc_noaa_curtain_sim_', xx(nn).tag, '.png'], style);

% "Analysis" pane
% ---------------
  figure(2); clf;
  scatter(lats, alts*1e-3, 8, omf2s, 's', 'filled');

  xlim(xx(nn).XLIMS);
  ylim(xx(nn).YLIMS);
  caxis(xx(nn).CLIMS);
  colormap(cmapd);
  grid on;
  set(gca, 'fontsize', 20);

% NB: leg2 and 'ANA'
  legend(leg2, 'location', 'northeast');
  text(0.03, 0.94, 'ANA', 'units', 'normalized', ...
       'fontsize', 30, 'fontweight', 'bold');

% Only for bottom pane
  xlabel('Latitude (^oN)'); ylabel('Altitude (km)');

  set(gca, 'Position', [0.1100, 0.4373-0.25, 0.7750+0.09, 0.4700+0.28]);
  hgexport(gcf, ['figs/m2cc_noaa_curtain_ana_', xx(nn).tag, '.png'], style);
end
