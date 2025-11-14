% Bug testing
fb = '../covid_free/covid_free.obspack_co2_mip_v3.mat';
fa = 'm2cc_ana.obspack_co2_mip_v3.mat';
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

% Move to plot settings
COLOR1 = [0   135 255]/255;
COLOR2 = [255 2   51 ]/255;
COLOR3 = [153 102 255]/255;

% Get gas name
fobs = cell_fobs{1};
ii   = strfind(fobs, '_');
gas  = fobs(1:ii-1);

slist = {[gas, '_abp_surface-flask_1_representative'],	...
         [gas, '_alt_surface-flask_1_representative'],	...
         [gas, '_ams_surface-flask_1_representative'],	...
         [gas, '_asc_surface-flask_1_representative'],	...
         [gas, '_avi_surface-flask_1_representative'],	...
         [gas, '_azr_surface-flask_1_representative'],	...
         [gas, '_bme_surface-flask_1_representative'],	...
         [gas, '_bmw_surface-flask_1_representative'],	...
         [gas, '_brw_surface-flask_1_representative'],	...
         [gas, '_cba_surface-flask_1_representative'],	...
         [gas, '_cgo_surface-flask_1_representative'],	...
         [gas, '_chr_surface-flask_1_representative'],	...
         [gas, '_crz_surface-flask_1_representative'],	...
         [gas, '_gmi_surface-flask_1_representative'],	...
         [gas, '_hba_surface-flask_1_representative'],	...
         [gas, '_ice_surface-flask_1_representative'],	...
         [gas, '_key_surface-flask_1_representative'],	...
         [gas, '_kum_surface-flask_1_representative'],	...
         [gas, '_mbc_surface-flask_1_representative'],	...
         [gas, '_mhd_surface-flask_1_representative'],	...
         [gas, '_mid_surface-flask_1_representative'],	...
         [gas, '_psa_surface-flask_1_representative'],	...
         [gas, '_rpb_surface-flask_1_representative'],	...
         [gas, '_shm_surface-flask_1_representative'],	...
         [gas, '_smo_surface-flask_1_representative'],	...
         [gas, '_spo_surface-flask_1_representative'],	...
         [gas, '_stm_surface-flask_1_representative'],	...
         [gas, '_syo_surface-flask_1_representative'],	...
         [gas, '_zep_surface-flask_1_representative']};
NLIST = numel(slist);

% Derive y-axis ordering
% ----------------------
ylist = zeros(NLIST, 1);
for nn = 1:NLIST
% Find site number
  sname = slist{nn};
  for ic = 1:NSITES
    kk = strfind(cell_fobs{ic}, sname);
    if (~isempty(kk)), break; end
  end

  if (isempty(kk))
    error(['Could not find matching site name: ', sname, ' ...']);
  end

% Get site latitude
  iokd  = find(cell_qcflags{ic}(1,:) == '.' & ...
               cell_qcflags{ic}(2,:) == '.');
  iokx  = find(cell_qcflags{ic}(1,:) == 'x' & ...
               cell_qcflags{ic}(2,:) == 'x');
  iok   = union(iokd, iokx);

  ylist(nn) = mean(cell_lat{ic}(iok));
end

% Order by latitude
[sval, groups] = sort(ylist);
slab = {};
for nn = 1:NLIST
  slab = {slab{:}, upper(slist{groups(nn)}(5:7))};
end

ll0 = 7*find(-23.43643 <= sval, 1, 'first');
ll1 = 7*find(sval <= 23.43643,  1, 'last');

% Build variables for box plots
% -----------------------------
nums  = [];
lats  = [];
mons  = [];

omf1s = [];
omf2s = [];

for nn = 1:NLIST
% Find site number
  sname = slist{groups(nn)};
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

  dvecs = datevec(cell_dnobs{ic});
  nums  = [nums;  nn*ones(size(iuse))];
  lats  = [lats;  mean(cell_lat{ic}(iuse))*ones(size(iuse))];
  mons  = [mons;  dvecs(iuse,2)];

% Apply additive mass fix to simulations
% --------------------------------------
  fix1  = 0.063*(cell_dnobs{ic} - datenum(2015,01,01))/365.25;
  fix2  = 0*cell_dnobs{ic};
  omf1  = cell_gasobs{ic}(iuse) - (cell_gasmod{ic}(iuse) - fix1(iuse));
  omf2  = cell_gasobs{ic}(iuse) - (cell_gasana{ic}(iuse) - fix2(iuse));

  omf1s = [omf1s; omf1];
  omf2s = [omf2s; omf2];
end

% Function for computing error tables
errs = @(ii) [mean(omf1s(ii)) mean(omf2s(ii)) ...
              sqrt(mean(omf1s(ii).^2)) sqrt(mean(omf2s(ii).^2))];

% Build 5 plots: ALL, DJF, MAM, JJA, SON
% --------------------------------------
for nn = 1:5
  figure(nn);
  clf; hold on;

% Lazy coding
  if (nn == 1), tag = 'all'; iplt = [1:numel(mons)]; end
  if (nn == 2), tag = 'djf'; iplt = find(mons == 12 | mons <=  2); end
  if (nn == 3), tag = 'mam'; iplt = find(2 < mons   & mons <=  5); end
  if (nn == 4), tag = 'jja'; iplt = find(5 < mons   & mons <=  8); end
  if (nn == 5), tag = 'son'; iplt = find(8 < mons   & mons <= 11); end

% The actual plots
  hh = boxplot(omf1s(iplt), nums(iplt), 'orientation', 'horizontal', ...
               'position', 7*nums(iplt)+1, 'colors', COLOR1, ...
               'plotstyle', 'compact', 'whisker', 0);
  set(hh(end,:), 'visible', 'off');
  hh = boxplot(omf2s(iplt), nums(iplt), 'orientation', 'horizontal', ...
               'position', 7*nums(iplt)-1, 'colors', COLOR2, ...
               'plotstyle', 'compact', 'whisker', 0);
  set(hh(end,:), 'visible', 'off');
  hold off;

% Build labels for groups with data
  set(gca, 'YTick', [1:max(nums(:))]*7);
  set(gca, 'YTickLabel', slab');

  xlim([-3 3]); ylim([0.5 max(nums(:))+0.5]*7);
  xl0 = xlim;
  xl1 = floor(max(abs(xl0)));

% Highlight tropics
  hold on;
  hp = plot(xl1*[-1 1], ll0*[1 1], '--k');
  uistack(hp, 'bottom');
  hp = plot(xl1*[-1 1], ll1*[1 1], '--k');
  uistack(hp, 'bottom');

% Hack to make plot look like it's on a grid
  hold on;
  grid on;
  for ll = 7:7:7*max(nums(:))
    hp = plot(xl1*[-1 1], ll*[1 1], 'color', [0.9 0.9 0.9]);
    uistack(hp, 'bottom');
  end
  hp = plot([0 0], ylim, 'k');
  uistack(hp, 'bottom');
  xlim(xl0);
  hold off;

% Hack to get the legend right
  hold on;
  hp = plot([0 0], ylim, 'color', COLOR2);
  uistack(hp, 'bottom');
  hp = plot([0 0], ylim, 'color', COLOR1);
  uistack(hp, 'bottom');
  hold off;

  legend('Simulated', 'Analysis', 'location', 'southeast');
  xlabel(['In situ ', upper(gas), ' (dry-air ', units, ')']);
  title('NOAA ObsPack OMF statistics');

% Beef up the fonts, indicate season, and widen
  set(gca, 'fontsize', 14);
  text(0.03, 0.92, upper(tag), 'units', 'normalized', ...
       'fontsize', 16, 'fontweight', 'bold');
  set(gca, 'Position', [0.1400 0.1100 0.7950 0.8150]);

  hgexport(gcf, ['figs/', expid, '_noaa_', gas, '_whiskers_', tag, '.eps']);

% A pared down figure for panels
  title(''); if nn > 1, legend off; xlabel(''); end
  hgexport(gcf, ['figs/', expid, '_noaa_', gas, '_whiskpan_', tag, '.eps']);

% Zonal-mean statistics, can copy into bottom of EPS file as
% a comment
  sout = ['(', upper(tag(1:3)), ')       BAND   BIAS/S  BIAS/A  RMSE/S  RMSE/A'; ...
          '-',           '---', '---------------------------------------------'];
  ibs  = intersect(find( 45 < lats), iplt);
  sout = [sout; ' 45 < Lat        ', sprintf(' %7.4f', errs(ibs))];
  ibs  = intersect(find( 23 < lats & lats <=  45), iplt);
  sout = [sout; ' 23 < Lat <=  45 ', sprintf(' %7.4f', errs(ibs))];
  ibs  = intersect(find(  0 < lats & lats <=  23), iplt);
  sout = [sout; '  0 < Lat <=  23 ', sprintf(' %7.4f', errs(ibs))];
  ibs  = intersect(find(-23 < lats & lats <=   0), iplt);
  sout = [sout; '-23 < Lat <=   0 ', sprintf(' %7.4f', errs(ibs))];
  ibs  = intersect(find(-45 < lats & lats <= -23), iplt);
  sout = [sout; '-45 < Lat <= -23 ', sprintf(' %7.4f', errs(ibs))];
  ibs  = intersect(find(             lats <= -45), iplt);
  sout = [sout; '      Lat <= -45 ', sprintf(' %7.4f', errs(ibs))];

  disp(' ');
  disp(sout);

  fid = fopen(['figs/', expid, '_noaa_', gas, '_whiskers_', tag, '.tsv'], 'w');
  fprintf(fid, '%s', sout);
  fclose(fid);
end
