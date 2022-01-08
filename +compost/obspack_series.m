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
  plot(dnums, obs, 'ko', dnums, sim, 'b.', dnums, ana, 'r.');
  xlim([datenum(2015,01,01) datenum(2020,01,01)]);
  datetick('x','keeplimits'); grid on
  legend('Obs', legb, lega, 'location', 'southeast', 'autoupdate', 'off');
  title(upper(strrep(tag(5:end), '_', ' ')));
  ylabel('CO2 (dry-air ppmv)');

% Plot flagged data in background
  hold on;
  yl = ylim;
  hno = plot(cell_dnobs{ic}(ino), cell_gasobs{ic}(ino), '.', 'color', [0.8 0.8 0.8]);
  ylim(yl);
  uistack(hno, 'bottom');
  hold off;

  drawnow;
  hgexport(gcf, ['figs/',tag,'.eps']);
end
