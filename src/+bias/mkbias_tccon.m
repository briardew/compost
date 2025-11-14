%MAKE_BIAS  Compare and plot model bias at NOAA ObsPack locations
%
%   I don't know

% Author: Brad Weir
%
% Changelog:
% 2018/04/06	Site is based on string now
%
% TODO:
%==============================================================================%

%==============================================================================%
dnout = [datenum(2009,01,01):1/24:datenum(2018,01,01)]';

biases = zeros(numel(dnout), NSITES);

% db, lh, ra, wg
%bcsites = {'co2_db', 'co2_ll', 'co2_ra', 'co2_wg'};
bcsites = {'co2_bi', 'co2_br', 'co2_gm', 'co2_js', 'co2_ka', 'co2_oc', 'co2_or', 'co2_pa', 'co2_rj'};

figure;
hold on;
for ii = 1:numel(bcsites)
  name = bcsites{ii};

% Find site number
% ----------------
  for ic = 1:NSITES
    kk = strfind(cell_fobs{ic}, name);
    if (~isempty(kk)), break; end
  end

  if (isempty(kk))
    error(['Could not find matching site name: ', name, ' ...']);
  end

% Determine quality flags
% -----------------------
  isok  = find(cell_xgaserr{ic} < 0.5);
  isok  = intersect(find(~isnan(cell_xgasalt{ic})), isok);

% Average repeat obs
% ------------------
  [dnfix,iin,iout] = unique(cell_dnobs{ic}(isok), 'stable');

  omf   = cell_xgasobs{ic}(isok) - cell_xgasalt{ic}(isok);
  omfix = accumarray(iout, omf, [], @mean); 

  bias = -fit.fourier(dnfix, omfix, [1:numel(dnfix)]');
  biases(:,ii) = interp1(dnfix, bias, dnout);
  codes{ii}    = upper(name(5:6));

  plot(dnfix, bias, '.');
end
hold off;
legend(codes);
grid on; datetick('x', 'keeplimits');
