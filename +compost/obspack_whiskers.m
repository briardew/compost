% Bug testing
%fb = input('Backgound filename: ', 's');
%fa = input('Analysis filename: ', 's');
fb = 'covid_free/20211019a/covid_free__obspack_co2_mip_v3.mat';
fa = 'm2cc-v1_ana/20211202a/m2cc_ana__obspack_co2_mip_v3.mat';

disp(['Comparing ', fb, ' to ']);
disp(['          ', fa, ' ...']);

load(fb);
va = load(fa);

cell_gasana = va.cell_gasmod;
clear va;

COLOR1 = [0   135 255]/255;
COLOR2 = [255 2   51 ]/255;
COLOR3 = [153 102 255]/255;

snames = {'co2_abp_surface-flask_1_representative',	...
          'co2_alt_surface-flask_1_representative',	...
          'co2_ams_surface-flask_1_representative',	...
          'co2_asc_surface-flask_1_representative',	...
          'co2_avi_surface-flask_1_representative',	...
          'co2_azr_surface-flask_1_representative',	...
          'co2_bme_surface-flask_1_representative',	...
          'co2_bmw_surface-flask_1_representative',	...
          'co2_brw_surface-flask_1_representative',	...
          'co2_cba_surface-flask_1_representative',	...
          'co2_cgo_surface-flask_1_representative',	...
          'co2_chr_surface-flask_1_representative',	...
          'co2_crz_surface-flask_1_representative',	...
          'co2_gmi_surface-flask_1_representative',	...
          'co2_hba_surface-flask_1_representative',	...
          'co2_ice_surface-flask_1_representative',	...
          'co2_key_surface-flask_1_representative',	...
          'co2_kum_surface-flask_1_representative',	...
          'co2_mbc_surface-flask_1_representative',	...
          'co2_mhd_surface-flask_1_representative',	...
          'co2_mid_surface-flask_1_representative',	...
          'co2_psa_surface-flask_1_representative',	...
          'co2_rpb_surface-flask_1_representative',	...
          'co2_shm_surface-flask_1_representative',	...
          'co2_smo_surface-flask_1_representative',	...
          'co2_spo_surface-flask_1_representative',	...
          'co2_stm_surface-flask_1_representative',	...
          'co2_syo_surface-flask_1_representative',	...
          'co2_zep_surface-flask_1_representative'};

omf1s = [];
lat1s = [];
num1s = [];

omf2s = [];
lat2s = [];
num2s = [];

ncomp = 0;
scomp = {};

for nn = 1:numel(snames)
  sname = snames{nn};

% Find site number
  for ic = 1:NSITES
    kk = strfind(cell_fobs{ic}, sname);
    if (~isempty(kk)), break; end
  end

  if (isempty(kk))
    error(['Could not find matching site name: ', sname, ' ...']);
  end

% *** Might need to sub-sample or something here ***
  iokd  = find(cell_qcflags{ic}(1,:) == '.' & ...
               cell_qcflags{ic}(2,:) == '.');
  iokx  = find(cell_qcflags{ic}(1,:) == 'x' & ...
               cell_qcflags{ic}(2,:) == 'x');
  iok   = union(iokd, iokx);
  iuse1 = intersect(find(~isnan(cell_gasmod{ic})), iok);
  iuse2 = intersect(find(~isnan(cell_gasana{ic})), iok);
  iuse  = intersect(iuse1, iuse2);

  if (~isempty(iuse))
    ncomp = ncomp + 1;
    scomp = {scomp{:}, snames{nn}};
  end

% A. Apply additive mass fix
% --------------------------
  fix1  = 0.063*(cell_dnobs{ic} - datenum(2015,01,01))/365.25;
  fix2  = 0*cell_dnobs{ic};
  omf1  = cell_gasobs{ic}(iuse) - (cell_gasmod{ic}(iuse) - fix1(iuse));
  omf2  = cell_gasobs{ic}(iuse) - (cell_gasana{ic}(iuse) - fix2(iuse));

  omf1s = [omf1s; omf1];
  lat1s = [lat1s; mean(cell_lat{ic}(iuse))*ones(size(iuse))];
  num1s = [num1s; ncomp*ones(size(iuse))];

  omf2s = [omf2s; omf2];
  lat2s = [lat2s; mean(cell_lat{ic}(iuse))*ones(size(iuse))];
  num2s = [num2s; ncomp*ones(size(iuse))];
end

[sval,sidx] = sort(lat1s);
slab = {};
nlab = unique(num1s(sidx), 'stable');
for ii = 1:numel(nlab)
  nn   = nlab(ii);
  npts = nnz(num1s == nn);
  slab = {slab{:},  upper(scomp{nn}(5:7))};
% slab = {slab{:}, [upper(scomp{nn}(5:7)), sprintf(' (%u)',npts)]};
end

figure;
hold off;
hp = plot([0 0], [-91 91], 'k');
hold on;
% The actual plots
hh = boxplot(omf1s, lat1s, 'orientation', 'horizontal', ...
             'position', lat1s+0.8, 'colors', COLOR1,   ...
             'plotstyle', 'compact', 'jitter', 0,       ...
             'outliersize', 2);
set(hh(end,:), 'visible', 'off');
hh = boxplot(omf2s, lat1s, 'orientation', 'horizontal', ...
             'position', lat1s-0.8, 'colors', COLOR2,   ...
             'plotstyle', 'compact', 'jitter', 0,       ...
             'outliersize', 2);
set(hh(end,:), 'visible', 'off');
% An extra for labels
hh = boxplot(omf2s, lat1s, 'orientation', 'horizontal', ...
             'position', lat1s, 'colors', 'k',          ...
             'plotstyle', 'compact', 'jitter', 0,       ...
             'outliersize', 2, 'labels', slab);
set(hh, 'visible', 'off');

xl0 = xlim;
xl1 = floor(max(abs(xl0)));
xl1 = 6; % overriding by hand
xlim(xl1*[-1 1]); ylim([-91 91]);

% Hack to make it look like a grid
grid on;
hp = plot(xl1*[-1 1], 0*[1 1], 'k:');
uistack(hp, 'bottom');
for ll = -90:10:90
  hp = plot(xl1*[-1 1], ll*[1 1], 'color', [0.9 0.9 0.9]);
  uistack(hp, 'bottom');
end

% Hacks to get the legend right
hp = plot([0 0], [-91 91], 'color', COLOR2);
uistack(hp, 'bottom');
hp = plot([0 0], [-91 91], 'color', COLOR1);
uistack(hp, 'bottom');
hold off;

legend('Simulated', 'Analysis', 'location', 'southwest');
xlabel('CO2 (dry-air ppmv)');
title('Marine boundary layer OMF statistics');

% Zonal band table
errs = @(ii) [mean(omf1s(ii)) mean(omf2s(ii)) ...
              sqrt(mean(omf1s(ii).^2)) sqrt(mean(omf2s(ii).^2))];

sout = ['      BIAS/S  BIAS/A  RMSE/S  RMSE/A'; ...
        '------------------------------------'];
ibs  = find( 45.0 < lat1s);
sout = [sout; 'NHHL ', sprintf('%7.4f %7.4f %7.4f %7.4f', errs(ibs))];
ibs  = find( 22.5 < lat1s & lat1s <= 45.0);
sout = [sout; 'NHML ', sprintf('%7.4f %7.4f %7.4f %7.4f', errs(ibs))];
ibs  = find(  0.0 < lat1s & lat1s <= 22.5);
sout = [sout; 'NHTR ', sprintf('%7.4f %7.4f %7.4f %7.4f', errs(ibs))];
ibs  = find(-22.5 < lat1s & lat1s <= 0);
sout = [sout; 'SHTR ', sprintf('%7.4f %7.4f %7.4f %7.4f', errs(ibs))];
ibs  = find(-45.0 < lat1s & lat1s <= -22.5);
sout = [sout; 'SHML ', sprintf('%7.4f %7.4f %7.4f %7.4f', errs(ibs))];
ibs  = find(lat1s <= -45.0);
sout = [sout; 'SHHL ', sprintf('%7.4f %7.4f %7.4f %7.4f', errs(ibs))];

disp(sout);
