% Bug testing
fb = '../covid_free/covid_free__tccon_co2_v21.mat';
fa = 'm2cc_ana__tccon_co2_v21.mat';
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

cell_xgasana = va.cell_xgasmod;
clear va;

% Move to plot settings
COLOR1 = [0   135 255]/255;
COLOR2 = [255 2   51 ]/255;
COLOR3 = [153 102 255]/255;

% Hack because it wasn't defined (bug has been fixed)
cell_fobs = {'co2_ae20120522_20181031.converted.nc', ...
             'co2_an20150202_20180418.converted.nc', ...
             'co2_bi20090301_20181001.converted.nc', ...
             'co2_br20100122_20201105.converted.nc', ...
             'co2_bu20170303_20200331.converted.nc', ...
             'co2_ci20120920_20201229.converted.nc', ...
             'co2_db20050828_20200430.converted.nc', ...
             'co2_df20130720_20201231.converted.nc', ...
             'co2_et20161007_20200906.converted.nc', ...
             'co2_eu20100724_20200706.converted.nc', ...
             'co2_fc20130316_20131004.converted.nc', ...
             'co2_gm20070716_20201127.converted.nc', ...
             'co2_hf20150918_20161231.converted.nc', ...
             'co2_if20120823_20121201.converted.nc', ...
             'co2_iz20070518_20210526.converted.nc', ...
             'co2_jc20070731_20080622.converted.nc', ...
             'co2_jf20110519_20180514.converted.nc', ...
             'co2_js20110728_20201229.converted.nc', ...
             'co2_ka20100419_20201130.converted.nc', ...
             'co2_lh20040629_20101209.converted.nc', ...
             'co2_ll20100202_20181031.converted.nc', ...
             'co2_lr20181003_20201231.converted.nc', ...
             'co2_ma20141001_20150624.converted.nc', ...
             'co2_ni20190831_20201129.converted.nc', ...
             'co2_oc20080706_20201228.converted.nc', ...
             'co2_or20090829_20201023.converted.nc', ...
             'co2_pa20040602_20201229.converted.nc', ...
             'co2_pr20140923_20200622.converted.nc', ...
             'co2_ra20110916_20200718.converted.nc', ...
             'co2_rj20131116_20190930.converted.nc', ...
             'co2_so20090516_20201020.converted.nc', ...
             'co2_sp20140406_20200925.converted.nc', ...
             'co2_tk20110804_20190930.converted.nc', ...
             'co2_wg20080626_20200630.converted.nc', ...
             'co2_zs20150424_20201127.converted.nc'};   

% Get gas name
fobs = cell_fobs{1};
ii   = strfind(fobs, '_');
gas  = fobs(1:ii-1);

% Lauder hack
cell_lat{    22} = [cell_lat{    20}; cell_lat{    21}; cell_lat{    22}];
cell_dnobs{  22} = [cell_dnobs{  20}; cell_dnobs{  21}; cell_dnobs{  22}];
cell_xgasobs{22} = [cell_xgasobs{20}; cell_xgasobs{21}; cell_xgasobs{22}];
cell_xgasmod{22} = [cell_xgasmod{20}; cell_xgasmod{21}; cell_xgasmod{22}];
cell_xgasana{22} = [cell_xgasana{20}; cell_xgasana{21}; cell_xgasana{22}];
cell_xgaserr{22} = [cell_xgaserr{20}; cell_xgaserr{21}; cell_xgaserr{22}];

[sval, sidx] = sort(cell_dnobs{22});

cell_lat{    22} = cell_lat{    22}(sidx);
cell_dnobs{  22} = cell_dnobs{  22}(sidx);
cell_xgasobs{22} = cell_xgasobs{22}(sidx);
cell_xgasmod{22} = cell_xgasmod{22}(sidx);
cell_xgasana{22} = cell_xgasana{22}(sidx);
cell_xgaserr{22} = cell_xgaserr{22}(sidx);

cell_lat{    20} = []; cell_lat{    21} = [];
cell_dnobs{  20} = []; cell_dnobs{  21} = [];
cell_xgasobs{20} = []; cell_xgasobs{21} = [];
cell_xgasmod{20} = []; cell_xgasmod{21} = [];
cell_xgasana{20} = []; cell_xgasana{21} = [];
cell_xgaserr{20} = []; cell_xgaserr{21} = [];

% JPL/CalTech hack
cell_lat{    06} = [cell_lat{    16}; cell_lat{    17}; cell_lat{    06}];
cell_dnobs{  06} = [cell_dnobs{  16}; cell_dnobs{  17}; cell_dnobs{  06}];
cell_xgasobs{06} = [cell_xgasobs{16}; cell_xgasobs{17}; cell_xgasobs{06}];
cell_xgasmod{06} = [cell_xgasmod{16}; cell_xgasmod{17}; cell_xgasmod{06}];
cell_xgasana{06} = [cell_xgasana{16}; cell_xgasana{17}; cell_xgasana{06}];
cell_xgaserr{06} = [cell_xgaserr{16}; cell_xgaserr{17}; cell_xgaserr{06}];

[sval, sidx] = sort(cell_dnobs{06});

cell_lat{    06} = cell_lat{    06}(sidx);
cell_dnobs{  06} = cell_dnobs{  06}(sidx);
cell_xgasobs{06} = cell_xgasobs{06}(sidx);
cell_xgasmod{06} = cell_xgasmod{06}(sidx);
cell_xgasana{06} = cell_xgasana{06}(sidx);
cell_xgaserr{06} = cell_xgaserr{06}(sidx);

cell_lat{    16} = []; cell_lat{    17} = [];
cell_dnobs{  16} = []; cell_dnobs{  17} = [];
cell_xgasobs{16} = []; cell_xgasobs{17} = [];
cell_xgasmod{16} = []; cell_xgasmod{17} = [];
cell_xgasana{16} = []; cell_xgasana{17} = [];
cell_xgaserr{16} = []; cell_xgaserr{17} = [];

% Four Corners/Indianapolis hack
cell_lat{    11} = []; cell_lat{    14} = [];
cell_dnobs{  11} = []; cell_dnobs{  14} = [];
cell_xgasobs{11} = []; cell_xgasobs{14} = [];
cell_xgasmod{11} = []; cell_xgasmod{14} = [];
cell_xgasana{11} = []; cell_xgasana{14} = [];
cell_xgaserr{11} = []; cell_xgaserr{14} = [];

slist = {'Ascension', ...
         'Anmey.', ...
         'Bialystok', ...
         'Bremen', ...
         'Burgos', ...
         'JPL', ...
         'Darwin', ...
         'Edwards', ...
         'E. Trout L.', ...
         'Eureka', ...
         'Four C.', ...
         'Garmisch', ...
         'Hefei', ...
         'Indy', ...
         ['Iza', char(241), 'a'], ...
         'JPL', ...
         'JPL', ...
         'Saga', ...
         'Karlsruhe', ...
         'Lauder', ...
         'Lauder', ...
         'Lauder', ...
         'Manaus', ...
         'Nicosia', ...
         'Lamont', ...
         ['Orl', char(233), 'ans'], ...
         'Park Falls', ...
         'Paris', ...
         ['R', char(233), 'union'], ...
         'Rikubetsu', ...
         ['Sodankyl', char(228)], ...
         ['Ny ', char(197), 'lesund'], ...
         'Tsukuba', ...
         'Wollon.', ...
         'Zugspitze'};

NLIST = numel(slist);
% Derive y-axis ordering
% ----------------------
ylist = zeros(NLIST, 1);
for nn = 1:NLIST
  ylist(nn) = mean(cell_lat{nn});
end

% Order by latitude
[sval, groups] = sort(ylist);
slab = {};
for nn = 1:NLIST
  slab = {slab{:}, slist{groups(nn)}};
end

ll0 = 7*find(-23.43643 <= sval, 1, 'first');
ll1 = 7*find(sval <= 23.43643,  1, 'last');

% Brutal hack to cut out the hack above (simplify this)
slab = {slab{1:end-6}};

% Build variables for box plots
% -----------------------------
nums  = [];
lats  = [];
mons  = [];

omf1s = [];
omf2s = [];

for nn = 1:NSITES
  ic = groups(nn);

  iok   = find(abs(cell_xgasobs{ic} - cell_xgasmod{ic}) < 7.*mean(cell_xgaserr{ic}));
  iuse1 = intersect(find(~isnan(cell_xgasmod{ic})), iok);
  iuse2 = intersect(find(~isnan(cell_xgasana{ic})), iok);
  iuse  = intersect(iuse1, iuse2);

  dvecs = datevec(cell_dnobs{ic});
  nums  = [nums;  nn*ones(size(iuse))];
  lats  = [lats;  mean(cell_lat{ic}(iuse))*ones(size(iuse))];
  mons  = [mons;  dvecs(iuse,2)];

% Apply additive mass fix to simulations
% --------------------------------------
  fix1  = 0.063*(cell_dnobs{ic} - datenum(2015,01,01))/365.25;
  fix2  = 0*cell_dnobs{ic};
  omf1  = cell_xgasobs{ic}(iuse) - (cell_xgasmod{ic}(iuse) - fix1(iuse));
  omf2  = cell_xgasobs{ic}(iuse) - (cell_xgasana{ic}(iuse) - fix2(iuse));

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
% set(gca, 'YTickLabel', {slab{2:end}}');
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
  xlabel(['Column ', upper(gas), ' (dry-air ', units, ')']);
  title('TCCON OMF statistics');

% Beef up the fonts, indicate season, and widen
  set(gca, 'fontsize', 14);
  text(0.03, 0.92, upper(tag), 'units', 'normalized', ...
       'fontsize', 16, 'fontweight', 'bold');
  set(gca, 'Position', [0.1400 0.1100 0.7950 0.8150]);

  hgexport(gcf, ['figs/', expid, '_tccon_', gas, '_whiskers_', tag, '.eps']);

% A pared down figure for panels
  title(''); if nn > 1, legend off; xlabel(''); end
  hgexport(gcf, ['figs/', expid, '_tccon_', gas, '_whiskpan_', tag, '.eps']);

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

  fid = fopen(['figs/', expid, '_tccon_', gas, '_whiskers_', tag, '.tsv'], 'w');
  fprintf(fid, '%s', sout);
  fclose(fid);
end
