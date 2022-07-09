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
cell_names = {'Ascension', ...
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

% Get gas name
fobs = cell_fobs{1};
ii   = strfind(fobs, '_');
gas  = fobs(1:ii-1);

% Lauder hack
cell_dnobs{  22} = [cell_dnobs{  20}; cell_dnobs{  21}; cell_dnobs{  22}];
cell_xgasobs{22} = [cell_xgasobs{20}; cell_xgasobs{21}; cell_xgasobs{22}];
cell_xgasmod{22} = [cell_xgasmod{20}; cell_xgasmod{21}; cell_xgasmod{22}];
cell_xgasana{22} = [cell_xgasana{20}; cell_xgasana{21}; cell_xgasana{22}];
cell_xgaserr{22} = [cell_xgaserr{20}; cell_xgaserr{21}; cell_xgaserr{22}];

[sval, sidx] = sort(cell_dnobs{22});

cell_dnobs{  22} = cell_dnobs{  22}(sidx);
cell_xgasobs{22} = cell_xgasobs{22}(sidx);
cell_xgasmod{22} = cell_xgasmod{22}(sidx);
cell_xgasana{22} = cell_xgasana{22}(sidx);
cell_xgaserr{22} = cell_xgaserr{22}(sidx);

cell_dnobs{  20} = []; cell_dnobs{  21} = [];
cell_xgasobs{20} = []; cell_xgasobs{21} = [];
cell_xgasmod{20} = []; cell_xgasmod{21} = [];
cell_xgasana{20} = []; cell_xgasana{21} = [];
cell_xgaserr{20} = []; cell_xgaserr{21} = [];

% JPL/CalTech hack
cell_dnobs{  06} = [cell_dnobs{  16}; cell_dnobs{  17}; cell_dnobs{  06}];
cell_xgasobs{06} = [cell_xgasobs{16}; cell_xgasobs{17}; cell_xgasobs{06}];
cell_xgasmod{06} = [cell_xgasmod{16}; cell_xgasmod{17}; cell_xgasmod{06}];
cell_xgasana{06} = [cell_xgasana{16}; cell_xgasana{17}; cell_xgasana{06}];
cell_xgaserr{06} = [cell_xgaserr{16}; cell_xgaserr{17}; cell_xgaserr{06}];

[sval, sidx] = sort(cell_dnobs{06});

cell_dnobs{  06} = cell_dnobs{  06}(sidx);
cell_xgasobs{06} = cell_xgasobs{06}(sidx);
cell_xgasmod{06} = cell_xgasmod{06}(sidx);
cell_xgasana{06} = cell_xgasana{06}(sidx);
cell_xgaserr{06} = cell_xgaserr{06}(sidx);

cell_dnobs{  16} = []; cell_dnobs{  17} = [];
cell_xgasobs{16} = []; cell_xgasobs{17} = [];
cell_xgasmod{16} = []; cell_xgasmod{17} = [];
cell_xgasana{16} = []; cell_xgasana{17} = [];
cell_xgaserr{16} = []; cell_xgaserr{17} = [];

for ic = 1:NSITES
  it  = strfind(cell_fobs{ic}, '_');
  tag = cell_fobs{ic}(1:it(1)+2);

  ioku = find(abs(cell_xgasobs{ic} - cell_xgasmod{ic}) < 7.*mean(cell_xgaserr{ic}));
  iokn = find(~isnan(cell_xgasmod{ic}) & ~isnan(cell_xgasana{ic}));
  iok  = intersect(ioku, iokn);
  ino  = setxor(iok, [1:numel(cell_dnobs{ic})]);

  if isempty(iok), continue; end

  drift = 0.063*(cell_dnobs{ic} - datenum(2015,01,01))/365.25;

% Old way
% dnums = cell_dnobs{ic}(iok);
% obs = cell_xgasobs{ic}(iok);
% sim = cell_xgasmod{ic}(iok) - drift(iok);
% ana = cell_xgasana{ic}(iok);

% New way
  dvecs = datevec(cell_dnobs{ic});
  davgs = dvecs(:,1)*100^3 + dvecs(:,2)*100^2 + dvecs(:,3)*100^1 ...
        + dvecs(:,4);
  [cc, id, it] = unique(davgs(iok), 'stable');

  dnums = accumarray(it,   cell_dnobs{ic}(iok), [], @mean);
  obs   = accumarray(it, cell_xgasobs{ic}(iok), [], @mean);
  sim   = accumarray(it, cell_xgasmod{ic}(iok) - drift(iok), [], @mean);
  ana   = accumarray(it, cell_xgasana{ic}(iok), [], @mean);

  errsa = [mean(obs - ana), std(obs - ana), sqrt(mean((obs - ana).^2))];
  errsb = [mean(obs - sim), std(obs - sim), sqrt(mean((obs - sim).^2))];

  lega = ['Ana (', num2str(errsa(1),'%.2f'), ', ', ...
                   num2str(errsa(3),'%.2f'), ')'];
  legb = ['Sim (', num2str(errsb(1),'%.2f'), ', ', ...
                   num2str(errsb(3),'%.2f'), ')'];

% Plot primary stuff
  plot(dnums, obs, 'k.');
  hold on;
  plot(dnums, sim, '.', 'color', COLOR1);
  plot(dnums, ana, '.', 'color', COLOR2);
  hold off;
  xlim([datenum(2015,01,01) datenum(2020,01,01)]);
  datetick('x','keeplimits'); grid on

% Hack to get the legend right
  yl0 = ylim;
  hold on;
  hp = plot([0 0], [2 3]*abs(yl0(end)), '.', 'color', COLOR2, 'markersize', 18);
  uistack(hp, 'bottom');
  hp = plot([0 0], [2 3]*abs(yl0(end)), '.', 'color', COLOR1, 'markersize', 18);
  uistack(hp, 'bottom');
  hp = plot([0 0], [2 3]*abs(yl0(end)), '.', 'color', 'k',    'markersize', 18);
  uistack(hp, 'bottom');
  hold off;
  ylim(yl0);

  legend('Obs', legb, lega, 'location', 'southeast', 'autoupdate', 'off');
  title(cell_names{ic});
  ylabel(['Column ', upper(gas), ' (dry-air ', units, ')']);

% Plot flagged data in background
% hold on;
% yl = ylim;
% hno = plot(cell_dnobs{ic}(ino), cell_xgasobs{ic}(ino), '.', 'color', [0.8 0.8 0.8]);
% ylim(yl);
% uistack(hno, 'bottom');
% hold off;

  set(gca, 'fontsize', 14);
  drawnow;
  hgexport(gcf, ['figs/tccon/', expid, '_tccon_', tag, '.eps']);
end
