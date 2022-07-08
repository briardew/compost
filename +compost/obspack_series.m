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

% Move to plot settings
COLOR1 = [0   135 255]/255;
COLOR2 = [255 2   51 ]/255;
COLOR3 = [153 102 255]/255;

% Get gas name
fobs = cell_fobs{1};
ii   = strfind(fobs, '_');
gas  = fobs(1:ii-1);

for ic = 1:NSITES
  tag = cell_fobs{ic}(1:end-3);

  iokd = find(cell_qcflags{ic}(1,:) == '.' & cell_qcflags{ic}(2,:) == '.');
  iokx = find(cell_qcflags{ic}(1,:) == 'x' & cell_qcflags{ic}(2,:) == 'x');
  iokn = find(~isnan(cell_gasmod{ic}) & ~isnan(cell_gasana{ic}));
  iok  = intersect(union(iokd, iokx), iokn);
  ino  = setxor(iok, [1:numel(cell_dnobs{ic})]);

  if isempty(iok), continue; end

  drift = 0.063*(cell_dnobs{ic}(iok) - datenum(2015,01,01))/365.25;

  dnums = cell_dnobs{ic}(iok);
  obs = cell_gasobs{ic}(iok);
  sim = cell_gasmod{ic}(iok) - drift;
  ana = cell_gasana{ic}(iok);

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
  title(upper(strrep(tag(5:end), '_', ' ')));
  ylabel(['In situ ', upper(gas), ' (dry-air ', units, ')']);

% Plot flagged data in background
  hold on;
  yl = ylim;
  hno = plot(cell_dnobs{ic}(ino), cell_gasobs{ic}(ino), '.', 'color', [0.8 0.8 0.8]);
  ylim(yl);
  uistack(hno, 'bottom');
  hold off;

  set(gca, 'fontsize', 14);
  drawnow;
  hgexport(gcf, ['figs/noaa/', expid, '_noaa_', tag, '.eps']);
end
